from kipoiseq import Interval
import numpy as np
import tensorflow as tf
import kipoiseq
from scipy.special import softmax
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from tqdm import tqdm
import pyfaidx



class GELU(tf.keras.layers.Layer):
    def __init__(self, name=None, **kwargs):
        super(GELU, self).__init__(**kwargs)

    def call(self, x):
        # return tf.keras.activations.sigmoid(1.702 * x) * x
        return tf.keras.activations.sigmoid(tf.constant(1.702) * x) * x


def load_model(model_path):
    return tf.keras.models.load_model(model_path, 
                                   custom_objects={"GELU": GELU}
                                  )

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
        # pyfaidx wants a 1-based interval which should be the same as those from R
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
    
    if final_length:
        exp_motif_range = kipoiseq.Interval(chrom, start, end).resize(final_length)
    else:
        exp_motif_range = kipoiseq.Interval(chrom, start, end)
        final_length = motif_length
    wt_seq = fasta_extractor.extract(exp_motif_range)
    wt_onehot = one_hot_encode(wt_seq)
    if perturb!='None':
        center = final_length // 2
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


def get_preds(model, df):
    preds = []
    for i, row in tqdm(df.iterrows()):
        wt = get_seq(row['seqnames'], row['start'], row['end'], 
                     final_length=2048)
        wt_pred = model.predict(wt[np.newaxis])
        preds.append(np.mean(wt_pred))
    return preds

def get_scores(model, df, perturb=0.25, normalize_by_wt=True):
    scores = []
    for i, row in tqdm(df.iterrows()):
        wt, mut = get_seq(row['seqnames'], row['start'], row['end'], 
                     final_length=2048, perturb=perturb)
        wt_pred = model.predict(wt[np.newaxis])
        mut_pred = model.predict(mut[np.newaxis])
        delta = wt_pred-mut_pred
        if normalize_by_wt:
            scores.append(((delta)/wt_pred).mean())
        else:
            scores.append((delta).mean())
    return scores
