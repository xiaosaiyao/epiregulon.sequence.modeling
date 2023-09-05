from kipoiseq import Interval
import numpy as np
import tensorflow as tf
import kipoiseq
from scipy.special import softmax
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from tqdm import tqdm
import pyfaidx

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



def load_cbp(model_path):
    return tf.keras.models.load_model(model_path, custom_objects={"multinomial_nll": multinomial_nll})

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
    

def one_hot_encode(sequence):
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)


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


def normalize_pred(model, seq, scale_preds=True):
    if len(seq.shape) == 2:
        seq = np.expand_dims(seq, axis=0)
    profile, count = model(seq)
    if scale_preds:
        return np.squeeze(softmax(profile) * count)
    else:
        return np.squeeze(softmax(profile)), np.squeeze(count)

def get_scores_bpnet(model, df, perturb=0.25):
    scores = []
    for i, row in tqdm(df.iterrows()):
        wt, mut = get_seq(row['seqnames'], row['start'], row['end'], 
                     final_length=2114, perturb=perturb)
        _, wt_pred = normalize_pred(model, wt, False)
        _, mut_pred = normalize_pred(model, mut, False)
        scores.append(((wt_pred-mut_pred)/wt_pred).flatten()[0])
    return scores

# def get_scores_bpnet(model, df, perturb=0.25):
#     scores_count = []
#     scores_scaled_profile_pcc = []
#     scores_scaled_profile_jsd = []

#     for i, row in tqdm(df.iterrows()):
#         wt, mut = get_seq(row['seqnames'], row['start'], row['end'], 
#                      final_length=2114, perturb=perturb)
#         wt_profile, wt_pred = normalize_pred(model, wt, False)
#         mut_profile, mut_pred = normalize_pred(model, mut, False)
#         scores_count.append(((wt_pred-mut_pred)/wt_pred).flatten()[0])
        
#         scaled_wt, scaled_mut = wt_profile*wt_pred, mut_profile*mut_pred
        
#         scores_scaled_profile_pcc.append(pearsonr(scaled_wt, scaled_mut)[0])
#         scores_scaled_profile_jsd.append(jensenshannon(scaled_wt, scaled_mut)*np.sign(mut_pred - wt_pred))
#     return scores_count, scores_scaled_profile_pcc, scores_scaled_profile_jsd