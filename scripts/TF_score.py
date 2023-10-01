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
from tqdm import tqdm

def dna_to_one_hot(seqs):
    """
    Converts a list of DNA ("ACGT") sequences to one-hot encodings, where the
    position of 1s is ordered alphabetically by "ACGT". `seqs` must be a list
    of N strings, where every string is the same length L. Returns an N x L x 4
    NumPy array of one-hot encodings, in the same order as the input sequences.
    All bases will be converted to upper-case prior to performing the encoding.
    Any bases that are not "ACGT" will be given an encoding of all 0s.

    Adapted from: Written by Alex Tseng https://gist.github.com/amtseng/010dd522daaabc92b014f075a34a0a0b
    """
    seq_len = len(seqs[0])
    assert np.all(np.array([len(s) for s in seqs]) == seq_len)

    # Join all sequences together into one long string, all uppercase
    seq_concat = "".join(seqs).upper() + "ACGT"
    # Add one example of each base, so np.unique doesn't miss indices later

    one_hot_map = np.identity(5)[:, :-1].astype(np.int8)

    # Convert string into array of ASCII character codes;
    base_vals = np.frombuffer(bytearray(seq_concat, "utf8"), dtype=np.int8)

    # Anything that's not an A, C, G, or T gets assigned a higher code
    base_vals[~np.isin(base_vals, np.array([65, 67, 71, 84]))] = 85

    # Convert the codes into indices in [0, 4], in ascending order by code
    _, base_inds = np.unique(base_vals, return_inverse=True)

    # Get the one-hot encoding for those indices, and reshape back to separate
    return one_hot_map[base_inds[:-4]].reshape((len(seqs), seq_len, 4))

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
def load_model_wrapper(args):
    import tensorflow as tf
    # read .h5 model
    custom_objects = {"multinomial_nll": multinomial_nll, "tf": tf}
    get_custom_objects().update(custom_objects)
    model = load_model(args)
    print("got the model")
    model.summary()
    return model


def createMotifSequences(df, fasta_ref):
    # Returns one-hot encoded list of motif sequences and lengths (both reference and ablated)
    lengths = []
    ref_list = []
    alt_list = []
    rev_ref_list = []
    rev_alt_list = []

    for index, row in tqdm(df.iterrows()):
        length = int(row[2] - row[1])
        if length % 2 == 1:
            length = length + 1
        pos = (length / 2) + row[1]
        chrom = row[0]
        start = pos - (2114 / 2)
        end = pos + (2114 / 2)
        
        # Store length for later normalizing - note that we need an even length so this is rounded
        lengths.append(length)

        refseq = fasta_ref.fetch(chrom, start, end).upper()
        seq = fasta_ref.fetch(chrom, start, end).upper()

        # Ablate the motif by filling in N for the length of the motif
        left = seq[:(1057 - int(length / 2))]
        right = seq[(1057 + int(length / 2)):]
        Nstring = "N" * length
        altseq = left + Nstring + right

        ref_seqs_one = dna_to_one_hot(refseq)
        alt_seqs_one = dna_to_one_hot(altseq)

        alt_z = np.reshape(alt_seqs_one, (1, 2114, 4))
        ref_z = np.reshape(ref_seqs_one, (1, 2114, 4))
        
        alt_list.extend(alt_z)
        ref_list.extend(ref_z)

        # Take the reverse complement, one-hot encode and store
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
        Rev_refseq = "".join(complement.get(base, base) for base in reversed(refseq))
        Rev_altseq = "".join(complement.get(base, base) for base in reversed(altseq))
        ref_seqs_one = dna_to_one_hot(Rev_refseq)
        alt_seqs_one = dna_to_one_hot(Rev_altseq)

        alt_z = np.reshape(alt_seqs_one, (1, 2114, 4))
        ref_z = np.reshape(ref_seqs_one, (1, 2114, 4))
        
        rev_alt_list.extend(alt_z)
        rev_ref_list.extend(ref_z)
    
    return ref_list, alt_list, rev_alt_list, rev_ref_list, lengths

def divide_chunks(l, n):
      
    # looping till length l
    for i in range(0, len(l), n): 
        yield l[i:i + n]

def tf_score(model, ref_list, alt_list, rev_alt_list, rev_ref_list):
    score = []
    jsd_score = []

    rev_score = []
    rev_jsd_score = []
    
    n = 50

    ref_list = np.array(ref_list)
    alt_list = np.array(alt_list)
    rev_alt_list = np.array(rev_alt_list)
    rev_ref_list = np.array(rev_ref_list)

    batches_ref = list(divide_chunks(ref_list, n))
    print(batches_ref[0].shape)
    batches_alt = list(divide_chunks(alt_list, n))
    batches_rev_alt = list(divide_chunks(rev_alt_list, n))
    batches_rev_ref = list(divide_chunks(rev_ref_list, n))

    for idx in tqdm(range(0, len(batches_ref))):
        with tf.device('/device:GPU:0'):
            # print("Here", idx)
            refpred_prob, refpred_cts = model.predict_on_batch(batches_ref[idx])
            altpred_prob, altpred_cts = model.predict_on_batch(batches_alt[idx])
            
            rev_refpred_prob, rev_refpred_cts = model.predict_on_batch(batches_rev_ref[idx])
            rev_altpred_prob, rev_altpred_cts = model.predict_on_batch(batches_rev_alt[idx])

        ref_prob_preds = softmax(refpred_prob)
        alt_prob_preds = softmax(altpred_prob)

        ref_logcount_preds=np.squeeze(np.array(refpred_cts))
        alt_logcount_preds=np.squeeze(np.array(altpred_cts))
        log_counts_diff = alt_logcount_preds - ref_logcount_preds
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
