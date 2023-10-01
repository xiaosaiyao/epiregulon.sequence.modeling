import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.utils import get_custom_objects
from kipoiseq import Interval
import pyfaidx
import numpy as np
import pickle
import kipoiseq
from tqdm import tqdm
import dinuc_shuffle
import itertools
import pybedtools
import os
import upsetplot
import matplotlib_venn
import matplotlib.pyplot as plt
from scipy.special import softmax
import logomaker


def get_scores_bpnet(model, df, perturb=0.25):
    scores = []
    for i, row in tqdm(df.iterrows()):
        wt, mut = get_seq(row['seqnames'], row['start'], row['end'], 
                     final_length=2114, perturb=perturb)
        _, wt_pred = normalize_pred(model, wt, False)
        _, mut_pred = normalize_pred(model, mut, False)
        scores.append(((wt_pred-mut_pred)/wt_pred).numpy().flatten()[0])
    return scores


def get_scores(model, df, perturb=0.25):
    scores = []
    for i, row in tqdm(df.iterrows()):
        wt, mut = get_seq(row['seqnames'], row['start'], row['end'], 
                     final_length=2048, perturb=perturb)
        wt_pred = model.predict(wt[np.newaxis])
        mut_pred = model.predict(mut[np.newaxis])
        scores.append(((wt_pred-mut_pred)/wt_pred).mean())
    return scores

def load_model(model_path):
    return tf.keras.models.load_model(model_path, 
                                   custom_objects={"multinomial_nll": multinomial_nll,
                                                  "GELU": GELU}
                                  )

def get_seq(chrom, start, end, final_length=None, perturb='None'):
    fasta_extractor = FastaStringExtractor()
    motif_length = end - start
    center = final_length // 2
    if final_length:
        exp_motif_range = kipoiseq.Interval(chrom, start, end).resize(final_length)
    else:
        exp_motif_range = kipoiseq.Interval(chrom, start, end)
    wt_seq = fasta_extractor.extract(exp_motif_range)
    wt_onehot = one_hot_encode(wt_seq)
    if perturb!='None':
        assert type(perturb) == int or type(perturb) == float, 'Bad perturbation'
        test_seqs = wt_onehot.copy()
        test_seqs[center-motif_length//2: center+motif_length//2, :] = perturb

        return wt_onehot, test_seqs
    else:
        return wt_onehot

class GELU(tf.keras.layers.Layer):
    def __init__(self, name=None, **kwargs):
        super(GELU, self).__init__(**kwargs)

    def call(self, x):
        # return tf.keras.activations.sigmoid(1.702 * x) * x
        return tf.keras.activations.sigmoid(tf.constant(1.702) * x) * x


# @tf.function
def saliency_map(X, model, class_index=None, func=tf.math.reduce_mean):
    """fast function to generate saliency maps"""
    if not tf.is_tensor(X):
        X = tf.Variable(X)

    with tf.GradientTape() as tape:
        tape.watch(X)
        profile, count = model(X)
        outputs = softmax(profile) * count

        
    return tape.gradient(outputs, X)



def plot_attribution_map(saliency_df, ax=None, figsize=(10,1)):
    """plot an attribution map using logomaker"""

    logomaker.Logo(saliency_df, figsize=figsize, ax=ax, fade_below=.5,)
    if ax is None:
        ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    plt.xticks([])
    plt.yticks([])
    
    
def normalize_pred(model, seq, scale_preds=True):
    if len(seq.shape) == 2:
        seq = np.expand_dims(seq, axis=0)
    profile, count = model(seq)
    if scale_preds:
        return softmax(profile) * count
    else:
        return softmax(profile), count


def bw_to_pd(bw):
    return pd.read_table(bw.fn, header=None)

def intersect_beds(b1, b2, v=False, wa=False, wb=False):
    if type(b1) == pd.core.frame.DataFrame: b1 = pybedtools.BedTool.from_dataframe(b1)
    if type(b2) == pd.core.frame.DataFrame: b2 = pybedtools.BedTool.from_dataframe(b2)

    intersect_result = b1.intersect(b2, v=v, wa=wa, wb=wb)
    df = bw_to_pd(intersect_result)

    return df

def plot_bed_intersect(beds, labels):
    matplotlib_venn.venn2(subsets = (
                          len(beds[0]), len(beds[1]), len(beds[0].intersect(beds[1]))), 
                          set_labels=labels)

def write_bed(df, path):
    df.to_csv(path, sep='\t', header=None, index=None)
    
def read_bed(path):
    return pd.read_csv(path, header=None, sep='\t')

def unique_combinations(elements, N):
    return list(itertools.combinations(elements, N))


# def cut_motif(seq, start, end, padding=''):
    
class FastaStringExtractor:
    
    def __init__(self, fasta_file='/gstore/project/lineage/shush/genomes/hg38.fa'):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()
    
def shuffle_str(seq):
    list_seq = list(seq)
    np.random.shuffle(list_seq)
    return "".join(list_seq)



def one_hot_encode(sequence):
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

def dict_save(d, path):
    with open(path, 'wb') as f:
        pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
        
        
def dict_load(dict_path):
    with open(dict_path, 'rb') as f:
        return pickle.load(f)
    
def adjust_range(df, new_L):
    start, end = df[1], df[2]
    center = start + (end - start) // 2

    exp_start = center - new_L//2
    exp_end = center + new_L//2
    df[1] = exp_start
    df[2] = exp_end
    return df

def load_cbp(model_path):
    return tf.keras.models.load_model(model_path, custom_objects={"multinomial_nll": multinomial_nll})



def predict_rows_of(df, model_path, fasta_extractor, L_input, save_prefix, N_subset=None):
    
    save_path = f'{save_prefix}_WT_model_{model_path.split("/")[-1].split(".")[0]}.csv'
    if os.path.isfile(save_path):
        print(f'File exists, skipping {save_path}')
        df = pd.read_csv(save_path)
    else:
        model = tf.keras.models.load_model(model_path, custom_objects={"multinomial_nll": multinomial_nll})
        df = df.reset_index(drop=True)
        if N_subset:
            np.random.seed(42)
            print('Random seed set')
            df = df.iloc[np.random.choice(list(df.index), N_subset, replace=False)].reset_index(drop=True)
        nrows = df.shape[0]
        print(f'N rows = {nrows}')

        results = []
        for i, one_row in tqdm(df.iterrows()):
            chrom, start, end = one_row[:3]
            if int(one_row[2] - one_row[1]) % 2 == 1:
                end = end + 1
            wt_seq = fasta_extractor.extract(kipoiseq.Interval(chrom, start, end).resize(L_input))
            wt_onehot = one_hot_encode(wt_seq)
            wt_pred = model.predict_on_batch(np.expand_dims(wt_onehot, axis=0))

            results.append(wt_pred[1][0][0])
        df['wt_count'] = results

        df.to_csv(save_path)
    
    return  df


def perturb_seqs(df, model_path, fasta_extractor, L_input, save_prefix, N=10,
                 offset=0):
    
    
    # result_path = f'{save_prefix}_{N}_model_{model_path.split("/")[-1].split(".")[0]}_shuffle_motif.csv'

    # if os.path.isfile(result_path):
    #     print('Loading from existing results')
    #     df = pd.read_csv(result_path)
    # else:
    model = tf.keras.models.load_model(model_path, custom_objects={"multinomial_nll": multinomial_nll})
    df = df.reset_index(drop=True)
    print(df.head())

    negative_control = []
    wt_control = []
    test = []

    for i, one_row in tqdm(df.iterrows()):
        chrom, start, end = one_row[:3]
        wt_seq = fasta_extractor.extract(kipoiseq.Interval(chrom, start, end).resize(L_input))
        motif_half_length = (end - start) // 2
        wt_onehot = one_hot_encode(wt_seq)
        whole_seq_shuffle = dinuc_shuffle.dinuc_shuffle(wt_seq, num_shufs=N) 
        whole_seq_dshuffle_onehot = np.array([one_hot_encode(w) for w in whole_seq_shuffle])
        wt_control.append(model.predict_on_batch(np.expand_dims(wt_onehot, axis=0))[1][0][0])
        negative_control.append(model.predict_on_batch(whole_seq_dshuffle_onehot)[1].mean()) # mean across shuffles

        test_seq = np.array([wt_onehot.copy() for _ in range(10)])


        test_seq[:,(L_input//2-motif_half_length + offset) : (L_input//2+motif_half_length+ offset),:] = whole_seq_dshuffle_onehot[:, (L_input//2-motif_half_length) : (L_input//2+motif_half_length),:]


        test.append(model.predict_on_batch(test_seq)[1].mean())

    df['WT'] = wt_control
    df['negative_control'] = negative_control
    df['test'] = test
        # df.to_csv(result_path)
    
    return df



def cut_motifs(chrom, start, end, motifs, model, padding, fasta_extractor, L_input):
    sum_L_motifs = 0
    for motif_start, motif_end in motifs:
        sum_L_motifs += motif_end - motif_start

    wt_seq = fasta_extractor.extract(kipoiseq.Interval(chrom, start, end).resize(L_input))
    wt_onehot = one_hot_encode(wt_seq)

    if padding == 'real':
        wt_seq_with_flanks = fasta_extractor.extract(kipoiseq.Interval(chrom, start-half_motif, end+half_motif).resize(L_input))
        mutant = wt_seq_with_flanks[:r_start] + wt_seq_with_flanks[r_end:]
    elif padding == 'edge':
        mutant = wt_seq[:half_motif] + wt_seq[:r_start] + wt_seq[r_end:] + wt_seq[-half_motif:]
   
    wt_pred = model.predict_on_batch(np.expand_dims(wt_onehot, axis=0))
    mutant_pred = model.predict_on_batch(np.expand_dims(mutant_onehot, axis=0))
    return wt_pred, mutant_pred


# bw = pyBigWig.open("/gstore/project/lineage/shush/convert_bw_unstranded.bw")
# vals = bw.values('chr1', )
# bw.close()