# TF_RE_mapping

Steps for reproducing the ROC curves:

1. save_motif_positions.R = Save motif positions as csv files (for use in the python script); output = 3 csvs (AR_689.csv, Gata6.Zf_109.csv, Nkx2.1.Homeobox_182.csv)  
2. score_all_tfs.ipynb = Score motif importance using Deep neural network (in this case Basenji) and save the scores; output = 3 csvs (basenji_AR_689.csv, basenji_Gata6.Zf_109.csv, basenji_Nkx2.1.Homeobox_182.csv)
3. reprogram_seq_chip_and_motif.R = compare ROC values using (i) public ChIP, (ii) motif annotation, (iii) DNN thresholds as methods for TF-RE overlapping for 2 motifs - GATA6 and NKX2-1 for 2 ways of clustering - by hash assignment and by 'cluster'
4. AR_chip_and_motif.R = compare ROC values using (i) public ChIP (ii) cell line matched ChIP (iii) motif annotations (iv) public ChIP that overlaps motif (v) cell line ChIP that overlaps motif for AR motif clustered by drug treatment - Enza vs DMSO

