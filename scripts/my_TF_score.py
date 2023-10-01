import pandas as pd
import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope
from scipy.special import softmax
from tensorflow import keras
from tensorflow.keras.utils import get_custom_objects
from scipy.spatial.distance import jensenshannon
import tensorflow as tf
import kipoiseq
from tqdm import tqdm
import os

def make_dir(dir_path):
    """ Make directory if doesn't exist."""
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    return dir_path


def one_hot_encode(sequence):
    """Convert sequence to one-hot."""
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)


#from https://github.com/kundajelab/basepair/blob/cda0875571066343cdf90aed031f7c51714d991a/basepair/losses.py#L87
def multinomial_nll(true_counts, logits):
  """Compute the multinomial negative log-likelihood
  Args:
   true_counts: observed count values
   logits: predicted logit values
  """
  counts_per_example = tf.reduce_sum(true_counts, axis=-1)
  dist = tfp.distributions.Multinomial(total_count=counts_per_example,
                     logits=logits)
  return (-tf.reduce_sum(dist.log_prob(true_counts)) /
      tf.cast(tf.shape(true_counts)[0], dtype=tf.float32))


#from chrombpnet repo 
def load_model_wrapper(args, summary=False):
    import tensorflow as tf
    # read .h5 model
    custom_objects = {"multinomial_nll": multinomial_nll, "tf": tf}
    get_custom_objects().update(custom_objects)
    model = load_model(args)
    print("got the model")
    if summary:
        model.summary()
    return model


def createMotifSequences(df, fasta_ref, rc, L_peak=None, L_input_seq=2114):
    # Returns one-hot encoded list of motif sequences and lengths (both reference and ablated)
    lengths = []
    
    
    N_seqs = df.shape[0]
    
    ref_seqs = np.empty((N_seqs, L_input_seq, 4))
    alt_seqs = np.empty((N_seqs, L_input_seq, 4))


    for index, row in tqdm(df.iterrows()):
        if L_peak:
            print('Setting peak length manually')
            length = L_peak
        else:
            length = int(row[2] - row[1]) # peak length
        
        if length % 2 == 1:
            length = length + 1
        pos = (length / 2) + row[1] # peak center
        chrom = row[0] # chromosome
        start = pos - (L_input_seq / 2) # half window left from center
        end = pos + (L_input_seq / 2) # half window right from center
        
        # Store length for later normalizing - note that we need an even length so this is rounded
        lengths.append(length)

        ref_seq = fasta_ref.fetch(chrom, start, end).upper() # get reference sequence
        
        if rc:
            # Take the reverse complement, one-hot encode and store
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'} 
            ref_seq = "".join(complement.get(base, base) for base in reversed(ref_seq))
        ref_seqs[index] = one_hot_encode(ref_seq) # onehot it
        
        # Ablate the motif by filling in N for the length of the motif <= I think they use the term motif as peak
        left = ref_seq[:(1057 - int(length / 2))]  # seq center now is 1057 
        right = ref_seq[(1057 + int(length / 2)):] # define left and right edges within the sequence where the peak/motif is not
        Nstring = "N" * length # fill with Ns
        alt_seq = left + Nstring + right # WT === NNN === WT
        alt_seqs[index] = one_hot_encode(alt_seq)
        


    return ref_seqs, alt_seqs, lengths

def divide_chunks(l, n):
      
    # looping till length l
    for i in range(0, len(l), n): 
        yield l[i:i + n]

def tf_score(model, ref_list, alt_list, rev_alt_list, rev_ref_list):
    score = []
    jsd_score = []

    rev_score = []
    rev_jsd_score = []
    
    n = 5


    batches_ref = list(divide_chunks(ref_list, n))
    batches_alt = list(divide_chunks(alt_list, n))
    batches_rev_alt = list(divide_chunks(rev_alt_list, n))
    batches_rev_ref = list(divide_chunks(rev_ref_list, n))
    print(batches_ref[1].shape)
    for idx in tqdm(range(0, len(batches_ref))):
        with tf.device('/device:GPU:0'):
            refpred_prob, refpred_cts = model.predict_on_batch(batches_ref[idx])
            altpred_prob, altpred_cts = model.predict_on_batch(batches_alt[idx])
            # RC
            rev_refpred_prob, rev_refpred_cts = model.predict_on_batch(batches_rev_ref[idx])
            rev_altpred_prob, rev_altpred_cts = model.predict_on_batch(batches_rev_alt[idx])
        # profile normalization
        ref_prob_preds = softmax(refpred_prob)
        alt_prob_preds = softmax(altpred_prob)
        
        ref_logcount_preds=np.squeeze(np.array(refpred_cts))
        alt_logcount_preds=np.squeeze(np.array(altpred_cts))
        log_counts_diff = alt_logcount_preds - ref_logcount_preds # delta counts
        # signed jensen shannon distance
        probs_jsd_diff = np.array([jensenshannon(x,y) for x,y in zip(alt_prob_preds, ref_prob_preds)])*np.sign(log_counts_diff)
        
        jsd_score.extend(probs_jsd_diff) 
        score.extend(altpred_cts - refpred_cts)
        
        # Repeat for reverse complement 
        ref_prob_preds = softmax(rev_refpred_prob)
        alt_prob_preds = softmax(rev_altpred_prob)

        ref_logcount_preds=np.squeeze(np.array(refpred_cts))
        alt_logcount_preds=np.squeeze(np.array(altpred_cts))
        log_counts_diff = alt_logcount_preds - ref_logcount_preds
        probs_jsd_diff = np.array([jensenshannon(x,y) for x,y in zip(alt_prob_preds, ref_prob_preds)])*np.sign(log_counts_diff)
    
        rev_jsd_score.extend(probs_jsd_diff) 
        rev_score.extend(altpred_cts - refpred_cts)
    
    score = np.array(score).flatten()
    jsd_score = np.array(jsd_score).flatten()
    rev_score = np.array(rev_score).flatten()
    rev_jsd_score = np.array(rev_jsd_score).flatten()

    return score, jsd_score, rev_score, rev_jsd_score

def runModel(df, fasta_ref, model_string, output):
    model = load_model_wrapper(model_string)
    ref_list, alt_list, rev_alt_list, rev_ref_list, lengths = createMotifSequences(df, fasta_ref)
    print("Got Sequences")
    score, jsd_score, rev_score, rev_jsd_score = tf_score(model,ref_list,alt_list,rev_alt_list, rev_ref_list)
    
    df['score'] = score
    df['jsd_score'] = jsd_score 
    df['rev_score']=rev_score
    df['rev_jsd_score']=rev_jsd_score
    df['length'] = lengths
    df.to_csv(output, sep='\t')
    return None
