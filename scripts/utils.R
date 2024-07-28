
calculate_ROC <- function(regulon, tf, positive_group_label, GeneExpressionMatrix,
                          label_column = "Clusters", negative_group_label= NULL,
                          plot_results = FALSE, control=FALSE,
                          exp_assay = "normalizedCounts"){
    # regulon <- regulon[regulon$tf == tf,]
    if (control){
      print('replacing weights with 1')
      regulon$weight = 1
    }
    activity.scores <- calculateActivity(regulon = regulon,
                                         expMatrix = GeneExpressionMatrix,
                                         exp_assay = exp_assay)
    activity.scores <- activity.scores[,colnames(GeneExpressionMatrix)]
    activity_values <- as.vector(activity.scores)

    # use clusters to define positive and negative group
    positive_group <- grep(positive_group_label, GeneExpressionMatrix[[label_column]])
    if(is.null(negative_group_label))
        negative_group <- setdiff(seq_len(ncol(GeneExpressionMatrix)), positive_group)
    else
        negative_group <- grep(negative_group_label, GeneExpressionMatrix[[label_column]])
    accuracy_stats <- calculate_accuracy_metrics(activity_values, positive_group, negative_group)
    if(plot_results) {
      plot(accuracy_stats$TPR~accuracy_stats$FPR, type = "l", xlab = "False positive rate", ylab = "True positive rate")
      hist(activity_values)
    }
    calculate_AUC(accuracy_stats$FPR,accuracy_stats$TPR)
    return(list("accuracy"=accuracy_stats, "activity"=activity_values))
}

#' @export
calculate_accuracy_metrics <- function(values, positive_elements_ind, negative_elements_ind, n_steps = 1e3){
    values <- as.vector(values)
    values <- values[unique(c(positive_elements_ind, negative_elements_ind))]
    positive_elements_ind <- which(unique(c(positive_elements_ind, negative_elements_ind)) %in% positive_elements_ind)
    max_val <- max(values)
    min_val <- min(values)
    max_val <- max_val+(max_val-min_val)/n_steps #add to account for ">=" used for threshold
    steps <- seq(min_val, max_val, length.out = n_steps)
    is_positive <- rep(FALSE, length(values))
    is_positive[positive_elements_ind] <- TRUE
    res_list <- list()
    all_combinations <- data.frame(threshold_reached=as.logical(c(1,1,0,0)), category = as.logical(c(1,0,1,0)))
    for(i in seq_along(steps)){
        threshold_reached <- values >= steps[i]
        confusion_matrix <- data.frame(threshold_reached =threshold_reached, category = is_positive)
        confusion_matrix <- rbind(all_combinations, confusion_matrix)
        tab <- table(confusion_matrix)
        # account for adding one observation for combination
        tab <- tab - 1
        res_list[[i]] <- c(TP = tab[2,2], FP = tab[2, 1], TN = tab[1, 1],
                           FN = tab[1, 2], cutoff = steps[i])
    }
    res <- as.data.frame(do.call(rbind, res_list))
    TPR <- res$TP/(res$TP + res$FN)
    FPR <- res$FP/(res$FP + res$TN)
    list(TPR = TPR, FPR = FPR, cutoff = res$cutoff, confusion_matrix_data = res)
}

randomizeWeights <- function(regulon, tf, positive_group_label,
                              n=100, weight_column = "weight",
                              label_column = "Clusters"){
    res <- c()
    mae <- scMultiome::reprogramSeq()
    GeneExpressionMatrix <- mae[["GeneExpressionMatrix"]]
    rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name
    for(i in seq_len(n)){
        message(sprintf("Randomization %d out of %d\n",i, n))
        regulon[[weight_column]] <- sample(regulon[[weight_column]])
        res[i] <- calculate_ROC(regulon, tf, positive_group_label, GeneExpressionMatrix,
                                label_column = label_column, negative_group_label= NULL,
                                plot_results = FALSE)

    }
    res
}

shuffle_weights = function(regulon, tf, pos_cluster,
                           cluster_by, n_shuffles = 10){
  regulon_shuffled = regulon
  dfs_shuffled = list()

  for(i in seq_len(n_shuffles)){
    message(sprintf("Randomization %d out of %d\n",i, n_shuffles))
    regulon_shuffled[['weight']] <- sample(regulon_shuffled[['weight']])
    res <- calculate_ROC(regulon_shuffled, tf, pos_cluster, GeneExpressionMatrix,
                         label_column = cluster_by, negative_group_label= NULL,
                         plot_results = FALSE)
    n_rows = length(res$accuracy$FPR)
    dfs_shuffled[[i]] <- data.frame('FPR'=res$accuracy$FPR,
                                    'TPR'=res$accuracy$TPR,
                                    'set'=rep(i, n_rows),
                                    'alpha' = rep(0.6, n_rows))
  }

  chip_res = calculate_ROC(regulon, tf, pos_cluster, GeneExpressionMatrix,
                           label_column = cluster_by, negative_group_label= NULL,
                           plot_results = FALSE)
  dfs_shuffled[['original']] <- data.frame('FPR'=chip_res$accuracy$FPR,
                                           'TPR'=chip_res$accuracy$TPR,
                                           'set'=rep(0, n_rows),
                                           'alpha' = rep(1, n_rows))

  dfs_shuffled_all = bind_rows(dfs_shuffled)
  return(dfs_shuffled_all)
}

#' @export
calculate_AUC <- function(x,y){
  y <- y[order(x)]
  x <- x[order(x)]
  non_unique_x_ind <- which(duplicated(x))
  non_unique_x_ind <- sort(unique(c(non_unique_x_ind, non_unique_x_ind-1)))
  y[non_unique_x_ind] <- sort(y[non_unique_x_ind])
  x_intervals <- diff(x)
  pair_mean_y <- (y[1:(length(y)-1)] + y[2:length(y)])/2
  sum(x_intervals*pair_mean_y)
}


getResultsFromActivity <- function(activity.matrix, add_plot=FALSE, tf,
                                   labels,
                                   positive_elements_label,
                                   negative_elements_label, ...){
  if(length(positive_elements_label) > 1)
    positive_elements_label <- paste(positive_elements_label, collapse = "|", sep ="")
  if(length(negative_elements_label) > 1)
    negative_elements_label <- paste(negative_elements_label, collapse = "|", sep ="")
  positive_elements_ind <- grep(positive_elements_label, labels)
  negative_elements_ind <- grep(negative_elements_label, labels)
  activity_values <- activity.matrix[tf, ]
  res <- calculate_accuracy_metrics(activity_values, positive_elements_ind,
                                    negative_elements_ind)
  if(add_plot)
    lines(res$TPR~res$FPR, ...)
  else
    plot(res$TPR~res$FPR, type ="l", xlab = "False postive rate",
         ylab = "True positive rate", ...)
  return(calculate_AUC(res$FPR, res$TPR))
}



get_activity_matrix <- function(method = NULL,
                                GRN = NULL,
                                n_bin =24,
                                tfs = NULL,
                                geneExprMatrix.sce = NULL){
  stopifnot(method %in% c("FigR", "Epiregulon", "cellOracle", "Pando", "Scenic"))
  library(epiregulon)
  if(method == "FigR"){
    # adjust gene names to FigR
    GRN <- GRN[,c("Motif", "DORC", "Score")]
    colnames(GRN) <- c("tf", "target", "weight")
    if(length(intersect(tfs, GRN$tf))==0) {
      warning("Tfs not found in the FigR output.")
      return(NULL)
    }
    return(calculateActivity(expMatrix = geneExprMatrix.sce,
                             regulon = GRN))
  }
  else if(method == "Epiregulon"){
    return(calculateActivity(expMatrix = geneExprMatrix.sce,
                             regulon = GRN))
  }
  else if(method == "Pando"){
    library(Signac)
    library(Seurat)
    library(Pando)
    Seurat_obj <- GRN
    test_srt <- find_modules(Seurat_obj, rsq_thresh = 0.05)
    TFmodules <- NetworkModules(test_srt)
    DefaultAssay(Seurat_obj) <- "RNA"
    tfs <- tfs[tfs %in% names(TFmodules@features$genes_pos)]
    if(length(tfs)==0) {
      warning("Tfs not found in the Pando output.")
      return(NULL)
    }
    activity.matrix <- matrix(nrow = length(tfs), ncol = dim(Seurat_obj)[2])
    for(i in seq_along(tfs)){
      tf <- tfs[i]
      # extract target genes (only upregulated)
      targets <- TFmodules@features$genes_pos[[tf]]
      # adjust gene name to Pando requirements
      tf <- tfs[i] <- gsub("-", "", tf)
      Seurat_obj <- AddModuleScore(Seurat_obj, features = list(targets), name = paste0(tf, "_activity"),
                                   nbin = n_bin)
      activity.matrix[i,] <- Seurat_obj@meta.data[[paste0(tf, "_activity1")]]
    }
    rownames(activity.matrix) <- tfs
    colnames(activity.matrix) <- colnames(Seurat_obj)
    return(activity.matrix)
  }
  else if(method == "cellOracle"){
    clusters <- colData(geneExprMatrix.sce)$cluster_cellOracle
    regulon <- GRN[,c("target","source", "coef_mean", "cluster")]
    colnames(regulon) <- c("target", "tf", "weight", "cluster")
    tfs <- tfs[tfs %in% regulon$tf]
    if(length(tfs)==0) {
      warning("Tfs not found in the cellOracle output.")
      return(NULL)
    }
    regulon <- regulon[regulon$tf %in% tfs,]
    regulon <- split(regulon, regulon$cluster)
    unique_clusters <- as.character(unique(clusters))
    cluster_cells <- split(colnames(geneExprMatrix.sce), clusters)
    activity.matrix <- matrix(ncol = ncol(assay(geneExprMatrix.sce)), nrow  = length(tfs), dimnames = list(tfs))
    cell_ind <- 1
    matrix_columns <- c()
    for (cluster_id in unique_clusters){
      selected_cells <- cluster_cells[[cluster_id]]
      activity_part <- calculateActivity(expMatrix = geneExprMatrix.sce[,selected_cells],
                                         regulon = regulon[[cluster_id]])
      if(is.null(activity_part))
        activity_part <- matrix(0, ncol = length(selected_cells),
                                nrow  = length(tfs),
                                dimnames = list(tfs))
      # add values for tfs which have 0 activity in the cluster
      if (nrow(activity_part) < nrow(activity.matrix)){
        missing_tfs <- setdiff(tfs, rownames(activity_part))
        activity_part <- rbind(activity_part, matrix(0, ncol = ncol(activity_part), nrow = length(missing_tfs), dimnames = list(missing_tfs)))
      }
      activity_part <- activity_part[tfs,,drop = FALSE]
      activity.matrix[,cell_ind:(cell_ind+length(selected_cells)-1)] <- as.matrix(activity_part)
      cell_ind <- cell_ind + length(selected_cells)
      matrix_columns <- c(matrix_columns, selected_cells)
    }
    # adjust row order to sce object
    colnames(activity.matrix) <- matrix_columns
    activity.matrix <- activity.matrix[Matrix::rowSums(activity.matrix)!=0,,drop= FALSE]
    return(activity.matrix)
  }
}
