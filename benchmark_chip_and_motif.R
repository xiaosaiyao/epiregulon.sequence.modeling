rm(list=ls())
library(epiregulon)
library(GenomicRanges)
library(compEpiTools)
source('./utils.R')
library(glue)
library("ggplot2") 
library(dplyr)
library(VennDiagram)

mae <- scMultiome::reprogramSeq()
# peak matrix
PeakMatrix <- mae[["PeakMatrix"]]
# expression matrix
GeneExpressionMatrix <- mae[["GeneExpressionMatrix"]]
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name
assay(GeneExpressionMatrix) <- as(log2(assay(GeneExpressionMatrix)+1), "CsparseMatrix")
reducedDimMatrix <- reducedDim(mae[['TileMatrix500']], "LSI_ATAC")
reducedDim(GeneExpressionMatrix, "UMAP_Combined") <- reducedDim(mae[['TileMatrix500']], "UMAP_Combined")

set.seed(1010)
p2g <- calculateP2G(peakMatrix = PeakMatrix,
                    expMatrix = GeneExpressionMatrix,
                    reducedDim = reducedDimMatrix)

tf <- 'NKX2-1'
if (tf=='NKX2-1'){
  pos_cluster <- 'C3'
  motif_name <-  'Nkx2.1.Homeobox_182'
} else if (tf=='GATA6'){
  pos_cluster <- 'C1'
  motif_name <-  'Gata6.Zf_109' 
  }
 

# get ChIP-seq data
grl_chip <- getTFMotifInfo(genome = "hg38")[tf]
# get motif positions
motif_positions_all <- readRDS("/gstore/project/ar_ligands/NE/reprogram_seq/multiome_arrayed/OUTPUT/doubletremoved/Annotations/Motif-Positions-In-Peaks.rds")
motif_positions<- motif_positions_all[[motif_name]]
write.csv(motif_positions, glue('{motif_name}.csv'))

evaluateRegulon <- function(grl, p2g, PeakMatrix, GeneExpressionMatrix, tf,
                            positive_class, add_weights, label,
                            normalize.expression=FALSE){
  overlap <- addTFMotifInfo(grl = grl, p2g = p2g, peakMatrix = PeakMatrix)
  #Number of TF-RE links
  print(glue('Number of TF-RE links in {label}: {length(overlap$idxATAC)}'))
  regulon <- getRegulon(p2g = p2g, overlap = overlap, aggregate = FALSE)
  if (add_weights){
    regulon <- addWeights(regulon = regulon,
                           expMatrix  = GeneExpressionMatrix,
                           exp_assay  = "logcounts",
                           peakMatrix = PeakMatrix,
                           peak_assay = "counts",
                           clusters = GeneExpressionMatrix$Clusters,
                           block_factor = NULL,
                           tf_re.merge = TRUE,
                           method = "corr")
  }else {
    print('SETTING WEIGHTS TO 1')
    regulon$weight <- 1
  }
  result <- calculate_ROC(regulon, tf, positive_class, GeneExpressionMatrix, 
                          normalize.expression=normalize.expression)

  print(glue('AUROC of {label}: {calculate_AUC(result$accuracy$FPR, result$accuracy$TPR)}'))
  return(list("regulon" = regulon, "roc" = result))
}


o <- findOverlaps(grl_chip[[tf]], motif_positions)
grl_chip_in_motifs <- grl_chip[[tf]][unique(queryHits(o)),]

##############add CBP threshold
dataset_all <- read.csv(glue("/gstore/project/lineage/shush/CBP_Epiregulon/analysis/bpnet_scores_{motif_name}_{pos_cluster}_all.csv"))

##############

gr_set = list(
              motif=motif_positions,
              chip_in_motifs=grl_chip_in_motifs 
              )
grls <- list('ChIP'=grl_chip)
for (l in names(gr_set)){
  g= GRangesList()
  g[[tf]] = gr_set[[l]]
  grls[[l]] = g
}

# for (t_i in quantile(dataset_all$bpnet_0.25_score, probs=c(0.005, 0.05, 0.5, 0.75, 0.95))){
#   dataset = dataset_all[dataset_all$bpnet_0.25_score>t_i,]
#   print(length(dataset$seqnames))
#   cbp_filtered <- GRanges(seqnames = Rle(dataset$seqnames),
#                           ranges = IRanges(dataset$start, dataset$end),
#                           strand = Rle("*"))
#   g= GRangesList()
#   g[[tf]] = cbp_filtered
#   grls[[glue('CBP {signif(t_i, 3)}')]] = g
# }



dfs_weighted = list()
dfs_1 = list()
regulons = list()

for (n in names(grls)){
  print(n)
  results_w_weight <- evaluateRegulon(grl=grls[[n]], p2g, PeakMatrix,
                                      GeneExpressionMatrix, tf, pos_cluster,
                                         TRUE, n)
  roc_w = signif(calculate_AUC(results_w_weight$roc$accuracy$FPR, results_w_weight$roc$accuracy$TPR), 3)
  
  results_no_weight <- evaluateRegulon(grl=grls[[n]], p2g, PeakMatrix, 
                                     GeneExpressionMatrix, tf, pos_cluster,
                                     FALSE, n)
  roc_no = signif(calculate_AUC(results_no_weight$roc$accuracy$FPR, results_no_weight$roc$accuracy$TPR), 3)
  
  regulons[[n]] = results_w_weight$regulon
  n_rows = length(results_w_weight$roc$accuracy$FPR)
  dfs_weighted[[glue('w{n}')]] <- data.frame('FPR'=results_w_weight$roc$accuracy$FPR, 
                                     'TPR'=results_w_weight$roc$accuracy$TPR,
                                     'set'=rep(glue('w{n} {roc_w}'), n_rows))
  
  dfs_1[[glue('{n}')]] <- data.frame('FPR'=results_no_weight$roc$accuracy$FPR, 
                                            'TPR'=results_no_weight$roc$accuracy$TPR,
                                            'set'=rep(glue('{n} {roc_no}'), n_rows))
}

dfs_weighted = bind_rows(dfs_weighted)
dfs_1 = bind_rows(dfs_1)


p=ggplot(data = dfs_weighted, aes(x = FPR, y = TPR, color=set)) + geom_line() +
  ggtitle(tf) + theme(legend.position = c(0.75,0.55))

p
p=ggplot(data = dfs_1, aes(x = FPR, y = TPR, color=set)) + geom_line() + ggtitle(tf) + theme(legend.position = c(0.75,0.55))
p      

# 
# # Venn diagram of regulon overlaps
# # RE overlap between ChIP and motif (and CBP)
# col_to_venn = 'idxATAC'
# venn_diagram_info = list()
# for (i in names(regulons)){
#   if (!(grepl('_', i, fixed=TRUE))){
#     venn_diagram_info[[i]] = regulons[[i]][[col_to_venn]]
#   }
#   
# }
# 
# temp <- venn.diagram(venn_diagram_info,
#                      # fill = c("red", "blue"), alpha = c(0.5, 0.5), 
#                      cex = 2, main=col_to_venn,
#                     filename = NULL)
# dev.off(dev.list()["RStudioGD"])
# grid.draw(temp)

