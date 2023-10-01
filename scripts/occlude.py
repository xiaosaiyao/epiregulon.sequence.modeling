import util
import dinuc_shuffle
import pandas as pd
import pysam
import kipoiseq
import numpy as np
from scipy.spatial import distance

L_input = 2114
method = 'dshuffle'
padding = 'edge'

df = pd.read_csv("motifs_in_peaks/C1_Motifs.bed", sep='\t', header=None)
print(df.shape)
print(df.head())
# model = load_model_wrapper('doubletRemoved_models/C1_chrombpnet_nobias.h5')

one_row = next(df.iterrows())[1]
chrom, start, end = one_row[:3]

if int(one_row[2] - one_row[1]) % 2 == 1:
    end = end + 1
L_motif = end - start
half_motif = L_motif // 2

fasta_extractor = util.FastaStringExtractor('../genomes/hg38.fa')
wt_seq = fasta_extractor.extract(kipoiseq.Interval(chrom, start, end).resize(L_input))
r_start = L_input // 2 - half_motif
r_end = L_input // 2 + half_motif

print(wt_seq[r_start:r_end])

if method == 'cut':
    if padding == 'real':
        wt_seq_with_flanks = fasta_extractor.extract(kipoiseq.Interval(chrom, start-half_motif, end+half_motif).resize(L_input))
        mutant = wt_seq_with_flanks[:start] + wt_seq_with_flanks[end:]
    elif padding == 'edge':
        mutant = wt_seq[:half_motif] + wt_seq[:r_start] + wt_seq[r_end:] + wt_seq[-half_motif:]
elif method == 'rshuffle':
    mutant = wt_seq[:r_start] + util.shuffle_str(wt_seq[r_start:r_end]) + wt_seq[r_end:]
elif method == 'dshuffle':
    for i in range(10):
        mutant = wt_seq[:r_start] + dinuc_shuffle.dinuc_shuffle(wt_seq[r_start:r_end]) + wt_seq[r_end:]

# print((len(mutant)))
# print(len(wt_seq))

