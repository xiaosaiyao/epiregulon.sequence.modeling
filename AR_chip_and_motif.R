library(scran)
library(scater)
library(epiregulon)
library(ArchR)
library(BiocParallel)
library(eulerr)
library(ggpubr)
library(dplyr)
library(glue)
source('./utils.R')

#checkdata
data_dir = '/gstore/project/ar_ligands/AR/scRNAseq/nonpipeline/'
archR_project_path <- "/gstore/project/ar_ligands/AR/scRNAseq/nonpipeline/OUTPUT/ArchRProject"
proj.all <- loadArchRProject(path = archR_project_path, showLogo = TRUE)

proj.all <- proj.all[!proj.all$hash_assignment2 %in% c(paste0("SAM24425416HTO-",14:16),
                                                       paste0("SAM24425417HTO-",1:10),
                                                       paste0("SAM24428812HTO-",1:6)),]




extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
# grl[["LNCaP"]][["AR"]] <- rtracklayer::import(glue("{data_dir}/data/SRX2017926.05.bed"), extraCols = extraCols_narrowPeak) #GSM2277148
# grl[["VCaP"]][["AR"]] <- rtracklayer::import(glue("{data_dir}/data/SRX471847.05.bed"), extraCols = extraCols_narrowPeak) #GSM1328945

motif_name = 'AR_689'
base_dir = '~/tf_re_mapping/'
motif_positions_all = readRDS('/gstore/project/ar_ligands/AR/scRNAseq/nonpipeline/OUTPUT/ArchRProject/Annotations/Motif-Positions-In-Peaks.rds')
motif_positions = motif_positions_all[[motif_name]]
motif_positions
write.csv(motif_positions, glue('{base_dir}/{motif_name}.csv')) # save for running CBP predictions

grl_public <- getTFMotifInfo(genome = "hg38")

# grl[["LNCaP"]][["FOXA1"]] <- rtracklayer::import(glue("{data_dir}/data/SRX1280790.05.bed"), extraCols = extraCols_narrowPeak) #GSM1891830
# grl[["VCaP"]][["FOXA1"]] <- rtracklayer::import(glue("{data_dir}/data/SRX2642364.05.bed"), extraCols = extraCols_narrowPeak)  #GSM2537225
# grl[["22Rv1"]][["FOXA1"]] <- rtracklayer::import(glue("{data_dir}/data/SRX5124035.05.bed"), extraCols = extraCols_narrowPeak)  #GSM3507235
# grl[["MDA"]][["FOXA1"]] <- rtracklayer::import("/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/diffbind/diffbind.allFOXA1.DBA_DESEQ2.bed")
# 
# grl[["MDA"]][["SMARCA4"]] <- rtracklayer::import("/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/diffbind/diffbind.allSMARCA4.DBA_DESEQ2.bed")


celllines <- c("LNCaP", "VCaP")



###############################
# load gene expression matrix
GeneExpressionMatrix <- getMatrixFromProject(
    ArchRProj = proj.all,
    useMatrix = "GeneExpressionMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
)



GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix, rename="normalizedCounts")
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name
GeneExpressionMatrix$TEST_ARTICLE <- factor(as.character(GeneExpressionMatrix$TEST_ARTICLE),
                                            levels = c("DMSO", "Enza","ARV110", "A9690", "7883", "4983") )




# Add embeddings
reducedDim(GeneExpressionMatrix, "UMAP_ATAC") <- getEmbedding(ArchRProj = proj.all,
                                                              embedding = "UMAP_ATAC",
                                                              returnDF = TRUE)[colnames(GeneExpressionMatrix), ]



#load peakMatrix
peakMatrix <- getMatrixFromProject(
    ArchRProj = proj.all,
    useMatrix = "PeakMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
)

peakMatrix <- as(peakMatrix, "SingleCellExperiment")
peakMatrix <- peakMatrix[, colnames(GeneExpressionMatrix)]
names(assays(peakMatrix)) <- "counts"

cell_bed_mapping = list(LNCaP='SRX2017926.05', 
                        VCaP='SRX471847.05')

reg_dir = '~/tf_re_mapping/regulons_fin/'
dir.create(reg_dir)
model_name = 'basenji'




for (cell in celllines){
    
    dataset_all <- read.csv(glue("{base_dir}/{model_name}_scores_{motif_name}_all.csv"))
  
    message(cell)
    # pick cell line specific ChIP
    bed_prefix=cell_bed_mapping[[cell]]
    grl_chip = rtracklayer::import(glue("{data_dir}/data/{bed_prefix}.bed"), extraCols = extraCols_narrowPeak)
    # find ChIP ranges overlapping with a motif
    o <- findOverlaps(grl_chip, motif_positions)
    grl_chip_in_motifs <- grl_chip[unique(queryHits(o)),]
    
    o = findOverlaps(grl_public$AR, motif_positions)
    grl_public_chip_in_motifs = grl_public$AR[unique(queryHits(o)),]
    # threshold DNN scores
    thresholds = quantile(dataset_all$count, probs=c(0.005, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.75))


    tf_mappings = list(public_chip=grl_public$AR, # public AR ChIP
                     chip=grl_chip, # cell line ChIP
                     motif=motif_positions, # motif ChIP
                     chip_in_motif=grl_chip_in_motifs, # ChIP overlapping motif
                     public_chip_in_motifs=grl_public_chip_in_motifs)

    for (t_i in thresholds){
      dataset = dataset_all[dataset_all$count>t_i,]
      
      print(length(dataset$seqnames))
      cbp_filtered <- GRanges(seqnames = Rle(dataset$seqnames),
                              ranges = IRanges(dataset$start, dataset$end),
                              strand = Rle("*"))
      tf_mappings[[glue('{model_name}-{signif(t_i, 3)}')]] = cbp_filtered
    }
    
    
    selected <- which(proj.all$Cell == cell )
    proj <- proj.all[selected,]
    GeneExpressionMatrix.select <- GeneExpressionMatrix[, selected]
    peakMatrix.select <- peakMatrix[, selected]


    # peak2gene links

    set.seed(1010)

    # find cell line specific peaks
    cell_peaks <- lapply(unique(proj$hash_assignment2),
                         function(cluster) {readRDS(file.path(archR_project_path, "PeakCalls", paste0(make.names(cluster), "-reproduciblePeaks.gr.rds")))})
    cell_peaks <- GenomicRanges::reduce(do.call(c, GRangesList(cell_peaks)))
    peaks_overlaps <- findOverlaps(cell_peaks, rowRanges(peakMatrix.select ))
    peakMatrix.cell <- peakMatrix.select[unique(subjectHits(peaks_overlaps)),]

    # find peak to gene links
    p2g <- calculateP2G(peakMatrix = peakMatrix.cell,
                        expMatrix = GeneExpressionMatrix.select ,
                        reducedDim = reducedDim(GeneExpressionMatrix.select ),
                        peak_assay = "counts",
                        exp_assay = "normalizedCounts",
                        cor_cutoff = 0.5
    )

    dfs_weighted = list()
    dfs_noweight = list()
    
    # Construct regulons
    for (label in names(tf_mappings)){
      print(label)
      grl <- list()
      grl[[cell]][["AR"]] = tf_mappings[[label]]

      overlap <- addTFMotifInfo(grl = GRangesList(grl[[cell]]),
                              p2g = p2g,
                              peakMatrix = peakMatrix.cell)



      regulon_df_full <- getRegulon(p2g, overlap, aggregate = FALSE)




      # prune network
      pruned.regulon <- pruneRegulon(regulon = regulon_df_full,
                                   expMatrix = GeneExpressionMatrix.select ,
                                   exp_assay = "normalizedCounts",
                                   peakMatrix = peakMatrix.cell ,
                                   peak_assay = "counts",
                                   prune_value = "pval",
                                   clusters = GeneExpressionMatrix.select$TEST_ARTICLE
      )
    
      # Add weights
      regulon.w <- addWeights(regulon = pruned.regulon,
                            expMatrix = GeneExpressionMatrix.select,
                            exp_assay = "normalizedCounts",
                            peakMatrix = peakMatrix.cell,
                            peak_assay = "counts",
                            method = "wilcoxon",
                            clusters = GeneExpressionMatrix.select$TEST_ARTICLE
      )

      # calculate activity
      score.combine <- calculateActivity(expMatrix = GeneExpressionMatrix.select,
                                       regulon = regulon.w,
                                       method = "weightedMean",
                                       exp_assay = "normalizedCounts",
                                       mode = "weight",
                                       FUN = "mean")
      res = calculate_ROC(regulon.w, 'AR', 'DMSO', GeneExpressionMatrix.select,
                          label_column='TEST_ARTICLE', negative_group_label='Enza',
                          scale_expression=FALSE, exp_assay = "normalizedCounts")
      roc_w = signif(calculate_AUC(res$accuracy$FPR, res$accuracy$TPR), 3)
      n_TG = length((unique(pruned.regulon$target)))
      dfs_weighted[[label]] <- data.frame('FPR'=res$accuracy$FPR, 
                                          'TPR'=res$accuracy$TPR,
                                          'set'=rep(glue('w{label} {n_TG} {roc_w}'), length(res$accuracy$FPR)))
      
      write.csv(regulon.w, glue('{reg_dir}/{cell}_{label}.csv'), quote=FALSE)
      regulon.noweights = regulon.w
      regulon.noweights$weight = 1
      res = calculate_ROC(regulon.noweights, 'AR', 'DMSO', GeneExpressionMatrix.select,
                          label_column='TEST_ARTICLE', negative_group_label='Enza',
                          scale_expression=FALSE, exp_assay = "normalizedCounts")
      roc_n = signif(calculate_AUC(res$accuracy$FPR, res$accuracy$TPR), 3)
      
      dfs_noweight[[label]] <- data.frame('FPR'=res$accuracy$FPR, 
                                          'TPR'=res$accuracy$TPR,
                                          'set'=rep(glue('{label} {n_TG} {roc_n}'), length(res$accuracy$FPR)))
      
    dfs_weighted_all = bind_rows(dfs_weighted) 
    plot_dir = glue("~/tf_re_mapping/plots_fin/")
    dir.create(plot_dir)
    
    pdf(glue("{plot_dir}/AR_ROC_plots_{cell}.pdf"))
    
    title = glue('AR {cell} with weights')
    p_weighted=ggplot(data = dfs_weighted_all, aes(x = FPR, y = TPR, color=set)) + geom_line() +
      ggtitle(title) + theme(legend.position = c(0.75,0.55)) + theme_classic()
    print(p_weighted)
    
    
    dfs_1_all = bind_rows(dfs_noweight) 
    title = glue('AR {cell} without weights')
    p_noweight=ggplot(data = dfs_1_all, aes(x = FPR, y = TPR, color=set)) + geom_line() +
      ggtitle(title) + theme(legend.position = c(0.75,0.55)) + theme_classic()
    print(p_noweight)
    dev.off() 
  }
}



# plot Venn diagrams
library(VennDiagram)
cell = 'LNCaP'
file_prefix = glue("{reg_dir}/{cell}_")

col_to_venn = 'idxATAC'
venn_diagram_info = list()
regulons = list()
for (reg_type in c('chip', 'public_chip')){
  reg = read.csv(glue('{file_prefix}{reg_type}.csv'))
  venn_diagram_info[[reg_type]] = reg[[col_to_venn]]
  
}



# for (i in names(regulons)){
#   if (!(grepl('_', i, fixed=TRUE))){
#     venn_diagram_info[[i]] = regulons[[i]][[col_to_venn]]
#   }
#   
# }

temp <- venn.diagram(venn_diagram_info,
                     # fill = c("red", "blue"), alpha = c(0.5, 0.5),
                     cex = 2, main=glue('{col_to_venn} {cell}'),
                    filename = NULL)
dev.off(dev.list()["RStudioGD"])
grid.draw(temp)



