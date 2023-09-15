rm(list=ls())
library(epiregulon)
library(GenomicRanges)
library(compEpiTools)

library(glue)
library("ggplot2") 
library(dplyr)
library(VennDiagram)
library(rstudioapi) 


#### ARGUMENTS
quantile_probability_thresholds = c(0.005, 0.05, 0.5, 0.75, 0.95)
model_name = "basenji"
plot_dir = glue("{base_dir}/ROC_plots")
dir.create(plot_dir)

base_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
source(glue('{base_dir}/utils.R'))

# motif positions from the ArchR project
motif_positions_all <- readRDS("/gstore/project/ar_ligands/NE/reprogram_seq/multiome_arrayed/OUTPUT/doubletremoved/Annotations/Motif-Positions-In-Peaks.rds")


# function to get ROC values for a given TF granges
evaluateRegulon <- function(grl, p2g, PeakMatrix, GeneExpressionMatrix, tf,
                            positive_class, add_weights, label, label_column,
                            scale_expression=FALSE){
  # TF RE overlap and regulon construction
  overlap <- addTFMotifInfo(grl = grl, p2g = p2g, peakMatrix = PeakMatrix)
  regulon <- getRegulon(p2g = p2g, overlap = overlap, aggregate = FALSE)
  if (add_weights){ # add weights 
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
  # regulon
  result <- calculate_ROC(regulon, tf, positive_class, GeneExpressionMatrix,
                          label_column=label_column)
  
  return(list("regulon" = regulon, "roc" = result))
}

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
# set tf of interest 

for (cluster_by in c('Clusters', 'hash_assignment')){  
  
  for (tf in c('GATA6', 'NKX2-1')){ 
    # set motif name and positive cluster name 
    if (tf=='NKX2-1'){
      if (cluster_by == 'Clusters'){
        pos_cluster <- 'C3'
      } else if (cluster_by == 'hash_assignment'){
        pos_cluster = 'HTO8_NKX2.1_UTR'
      }
      
      motif_name <-  'Nkx2.1.Homeobox_182'
    } else if (tf=='GATA6'){
      if (cluster_by == 'Clusters'){
        pos_cluster <- 'C1'
      } else if (cluster_by == 'hash_assignment'){
        pos_cluster = 'GATA6'
      }
      motif_name <-  'Gata6.Zf_109' 
      }
     
    
    # get ChIP-seq data
    grl_chip <- getTFMotifInfo(genome = "hg38")[tf] # public ChIP
    # get motif positions
    motif_positions<- motif_positions_all[[motif_name]]

    # find ChIP ranges overlapping with a motif
    o <- findOverlaps(grl_chip[[tf]], motif_positions)
    grl_chip_in_motifs <- grl_chip[[tf]][unique(queryHits(o)),]
    
    #load DNN model scores
    dataset_all <- read.csv(glue("{base_dir}/scores/{model_name}_{motif_name}.csv"))
    
    
    # assemble ChIP, motif and ChIP in motif granges
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
    # assemble DNN thresholding based grl-s
    thresholds = quantile(dataset_all$count, probs=quantile_probability_thresholds)
    for (t_i in thresholds){
      dataset = dataset_all[dataset_all$count>t_i,]
      print(length(dataset$seqnames))
      cbp_filtered <- GRanges(seqnames = Rle(dataset$seqnames),
                              ranges = IRanges(dataset$start, dataset$end),
                              strand = Rle("*"))
      g= GRangesList()
      g[[tf]] = cbp_filtered
      grls[[glue('{model_name} {signif(t_i, 3)}')]] = g
    }
    
    
    dfs_weighted = list() # to save weighted accuracy results
    dfs_1 = list() # to save results with weights set to 1
    regulons = list() # to save the regulons
    
    for (n in names(grls)){
      print(n)
      # results with weights
      results_w_weight <- evaluateRegulon(grl=grls[[n]], p2g, PeakMatrix,
                                          GeneExpressionMatrix, tf, pos_cluster,
                                             TRUE, n, label_column=cluster_by)
      # ROC values for weighted regulon
      roc_w = signif(calculate_AUC(results_w_weight$roc$accuracy$FPR, results_w_weight$roc$accuracy$TPR), 3)
      # results with weights set to 1
      results_no_weight <- evaluateRegulon(grl=grls[[n]], p2g, PeakMatrix, 
                                         GeneExpressionMatrix, tf, pos_cluster,
                                         FALSE, n, label_column=cluster_by)
      # ROC values for regulon with weights set to 1
      roc_no = signif(calculate_AUC(results_no_weight$roc$accuracy$FPR, results_no_weight$roc$accuracy$TPR), 3)
      
      regulons[[n]] = results_w_weight$regulon # keep the regulon to plot Venn diagrams
      n_rows = length(results_w_weight$roc$accuracy$FPR)
      n_TG = length(unique(results_w_weight$regulon$target))
      # save weighted regulon
      dfs_weighted[[glue('w{n}')]] <- data.frame('FPR'=results_w_weight$roc$accuracy$FPR, 
                                         'TPR'=results_w_weight$roc$accuracy$TPR,
                                         'set'=rep(glue('w{n} {roc_w}'), n_rows))
      # save regulon with weights = 1
      dfs_1[[glue('{n}')]] <- data.frame('FPR'=results_no_weight$roc$accuracy$FPR, 
                                                'TPR'=results_no_weight$roc$accuracy$TPR,
                                                'set'=rep(glue('{n} {roc_no}'), n_rows))
    }
    # combine results
    dfs_weighted = bind_rows(dfs_weighted) 
    dfs_1 = bind_rows(dfs_1)
    
    
    dfs_shuffled_all = shuffle_weights(regulons$ChIP, tf, pos_cluster, cluster_by)
    


    pdf(glue("{plot_dir}/reprogram_seq_ROC_plots_{tf}_{cluster_by}.pdf"))
    
    # plot 
    
    p_weighted=ggplot(data = dfs_weighted, aes(x = FPR, y = TPR, color=set)) + geom_line() +
      ggtitle(tf) + theme(legend.position = c(0.75,0.55)) + theme_classic()
    p_1=ggplot(data = dfs_1, aes(x = FPR, y = TPR, color=set)) + geom_line() +
      ggtitle(tf) + theme(legend.position = c(0.75,0.55)) + theme_classic()
    print(p_weighted) 
    print(p_1)     
    
    non_cbp = Filter(function(x) !any(grepl(model_name, x)), unique(dfs_weighted$set))
    dfs_weighted_non_cbp = subset(dfs_weighted, set %in% non_cbp)
    non_cbp = Filter(function(x) !any(grepl(model_name, x)), unique(dfs_1$set))
    dfs_1_non_cbp = subset(dfs_1, set %in% non_cbp)
    
    p_weighted=ggplot(data = dfs_weighted_non_cbp, aes(x = FPR, y = TPR, color=set)) + geom_line() +
      ggtitle(tf) + theme(legend.position = c(0.75,0.55)) + theme_classic()
    p_1=ggplot(data = dfs_1_non_cbp, aes(x = FPR, y = TPR, color=set)) + geom_line() +
      ggtitle(tf) + theme(legend.position = c(0.75,0.55)) + theme_classic()
    print(p_weighted) 
    print(p_1)   
    
    shuffled_p = ggplot(data = dfs_shuffled_all, aes(x = FPR, y = TPR, group=set, alpha=alpha)) + geom_line(show.legend = FALSE) +
                        ggtitle(glue('{tf} ChIP')) + theme_classic() 
    print(shuffled_p)
    dev.off() 
  }}





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

