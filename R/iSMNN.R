#' @include iSMNN_utility.R
#'
NULL


#' @title iSMNN
#'
#' @description This function iSMNN is designed to perform iterative supervised batch effect correction for scRNA-seq data by refining mutual nearest neighbors (MNNs) within corresponding clusters (or cell types) on the top of corrected data.
#' It takes as input raw expression matrices from two or more batches and a list of the unified cluster labels (output from unifiedClusterLabelling of SMNN package).
#' It outputs a Seurat object that contains the the batch-corrected expression matrix for batches
#' @usage iSMNN(object.list = merge.list, batch.cluster.labels = batch.cluster.labels,
#'              matched.clusters = c("Endothelial cells", "Macrophage", "Fibroblast"),
#'              strategy = "Short.run", iterations = 5, dims = 1:20, npcs = 30)
#'
#' @param object.list A list of \code{\link{Seurat}} objects between which to find anchors for downstream integration.
#' @param assay A vector of assay names specifying which assay to use when constructing anchors. If NULL, the current default assay for each object is used.
#' @param batch.cluster.labels is a list of vectors specifying the cluster labels of each cell from each batch. Cells not belonging to any clusters should be set to 0.
#' @param matched.clusters specifies the cell clusters matched between two or more batches.
#' @param strategy specifies the iteration option chosen for batch effect correction that in the first option "\code{Short.run}", iSMNN runs for a fixed number of iterations (default = 5) and takes the output with the lowest F statistic as the optimal correction results;  
#' in the second option "\code{Long.run}", after the first local minimum is observed, an additional number of iterations (default = 3) is run to allow leveraging possible further decrease of F statistic after the first local minimal value.
#' @param iterations defines the number of iterations to execute.
#' @param reference A vector specifying the object/s to be used as a reference during integration. If NULL (default),
#' all pairwise anchors are found (no reference/s). If not NULL, the corresponding objects in \code{object.list}
#' will be used as references. When using a set of specified references, anchors are first found between each query and each reference
#' The references are then integrated through pairwise integration. Each query is then mapped to the integrated reference
#' @param anchor.features Can be either:
#' \itemize{
#'   \item{A numeric value. This will call \code{\link{SelectIntegrationFeatures}} to select the provided number of features to be used in anchor finding}
#'   \item{A vector of features to be used as input to the anchor finding process}
#' }
#' @param scale Whether or not to scale the features provided. Only set to FALSE
#' if you have previously scaled the features you want to use for each object in
#' the object.list
#' @param reduction Dimensional reduction to perform when finding anchors. Can
#' be one of:
#' \itemize{
#'   \item{cca: Canonical correlation analysis}
#'   \item{rpca: Reciprocal PCA}
#' }
#' @param l2.norm Perform L2 normalization on the CCA cell embeddings after dimensional reduction
#' @param dims Which dimensions to use from the CCA to specify the neighbor search space
#' @param k.anchor How many neighbors (k) to use when picking anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the neighborhood search space in the anchor filtering
#' @param nn.method Method for nearest neighbor finding. Options include: rann, annoy
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param k.weight Number of neighbors to consider when weighting. Default is \code{k.weight = 100}
#' @param sd.weigth defines the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell. Default is \code{sd.weigth = 1}
#' @param verbose Print progress bars and output
#'
#' @return iSMNN returns a Seurat object that contains the the batch-corrected expression matrix for batches
#'
#' @author Yuchen Yang <yyuchen@email.unc.edu>, Gang Li <franklee@live.unc.edu>, Li Qian <li_qian@med.unc.edu>, Yun Li <yunli@med.unc.edu>
#' @references Yuchen Yang, Gang Li, Li Qian, Yun Li. iSMNN 2020
#'
#' @examples
#' # Load the example data data_SMNN
#' data("data_iSMNN")
#'
#' # Provide the marker genes for cluster matching
#' markers <- c("Col1a1", "Pdgfra", "Ptprc", "Pecam1")
#'
#' # Specify the cluster labels for each marker gene
#' cluster.info <- c("fibroblast", "fibroblast", "macrophage", "endothelial cells")
#'
#' # Harmonize cluster labels across batches
#' library(SMNN)
#' batch.cluster.labels <- unifiedClusterLabelling(batches = list(data_SMNN$batch1.mat, data_iSMNN$batch2.mat), features.use = markers,
#'                                                 cluster.labels = cluster.info, min.perc = 0.3)
#' names(batch.cluster.labels[[1]]) <- colnames(data_iSMNN$batch1.mat)
#' names(batch.cluster.labels[[2]]) <- colnames(data_iSMNN$batch2.mat)
#'
#' # Construct the input object for batches using Seurat
#' library(Seurat)
#' merge <- CreateSeuratObject(counts = cbind(data_iSMNN$batch1.mat, data_iSMNN$batch2.mat), min.cells = 0, min.features = 0)
#' batch_id <- c(rep("batch1", ncol(data_iSMNN$batch1.mat)), rep("batch2", ncol(data_iSMNN$batch2.mat)))
#' names(batch_id) <- colnames(merge)
#' merge <- AddMetaData(object = merge, metadata = batch_id, col.name = "batch_id")
#' merge.list <- SplitObject(merge, split.by = "batch_id")
#'
#' merge.list <- lapply(X = merge.list, FUN = function(x) {
#'   x <- NormalizeData(x)
#'   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#' })
#'
#' # Correct batch effect
#' corrected.results <- iSMNN(object.list = merge.list, batch.cluster.labels = batch.cluster.labels,
#'                            matched.clusters = c("endothelial cells", "macrophage", "fibroblast"),
#'                            strategy = "Short.run", iterations = 5, dims = 1:20, npcs = 30, k.filter = 30)
#'
#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @import Seurat
#' @import SMNN
#'
#' @export
iSMNN <- function(
  object.list = NULL,
  assay = NULL,
  batch.cluster.labels = NULL,
  matched.clusters = NULL,
  strategy = c("Short.run", "Long.run"),
  iterations = NULL,
  reference = NULL,
  anchor.features = 2000,
  scale = TRUE,
  reduction = c("cca", "rpca"),
  l2.norm = TRUE,
  dims = 1:30,
  npcs = 30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "rann",
  eps = 0,
  k.weight = 100,
  sd.weight = 1,
  verbose = TRUE
) {
  strategy <- match.arg(arg = strategy)
  reduction <- match.arg(arg = reduction)
  if (reduction == "rpca") {
    reduction <- "pca"
  }
  message("1st round of batch effect correcting ...")
  object.anchors_1round <- iSMNN_FindSMNNs(object.list = object.list, batch.cluster.labels = batch.cluster.labels, matched.clusters = matched.clusters, assay = assay,
                                           reference = reference, anchor.features = anchor.features, scale = scale, normalization.method = "LogNormalize",
                                           reduction = reduction, l2.norm = l2.norm, dims = dims, k.anchor = k.anchor, k.filter = k.filter,
                                           k.score = k.score, max.features = max.features, nn.method = nn.method, eps = eps, verbose = verbose)

  k.weight_input = k.weight
  if (nrow(object.anchors_1round@anchors)/2 <= k.weight){
    k.weight = nrow(object.anchors_1round@anchors)/2-1
  }
  message("Integrating data ...")
  object.combined_1round <- IntegrateData(anchorset = object.anchors_1round, dims = dims, k.weight = k.weight, sd.weight = sd.weight, verbose = FALSE)
  DefaultAssay(object.combined_1round) <- "integrated"

  # Run the standard workflow for visualization and clustering
  object.combined_1round <- ScaleData(object.combined_1round, verbose = FALSE)
  object.combined_1round <- RunPCA(object.combined_1round, npcs = npcs, verbose = FALSE)

  message("UMAP embedding ...")
  object.combined_1round <- RunUMAP(object.combined_1round, reduction = "pca", dims = 1:10, verbose = FALSE)
  object.combined_1round <- AddMetaData(object = object.combined_1round, metadata = unlist(batch.cluster.labels), col.name = "cell.anno")

  ### F measure for 1st round results
  data = data.frame(object.combined_1round@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,])
  data$Celltype = factor(unlist(batch.cluster.labels)[object.combined_1round@meta.data$cell.anno %in% matched.clusters])
  batch = c()
  for(i in 1:length(unique(object.combined_1round@meta.data$batch_id[object.combined_1round@meta.data$cell.anno %in% matched.clusters]))){
    batch = c(batch, rep(i, table(object.combined_1round@meta.data$batch_id[object.combined_1round@meta.data$cell.anno %in% matched.clusters])[unique(object.combined_1round@meta.data$batch_id[object.combined_1round@meta.data$cell.anno %in% matched.clusters])[i]]))
  }
  data$batch = factor(batch)

  if (length(levels(data$Celltype))>1){
    manova_1round = manova(object.combined_1round@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,] ~ Celltype+batch, data = data)
  } else{
    manova_1round = manova(object.combined_1round@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,] ~ batch, data = data)
  }
  temp_1round = summary(manova_1round)$stats[1:2,3]
  Fbatch_1round = as.data.frame(temp_1round)["batch",]

  # iterations
  if (iterations > 1){
     if (strategy == "Short.run"){
        object.combined_desired = object.combined_1round
        Fbatch_desired = Fbatch_1round
        for (k in 2:iterations){
          message(k, " round of batch effect correcting ...")
          corrected.object = SplitObject(object.combined_desired, split.by = "batch_id")
          object.anchors_after_correction <- iSMNN_FindSMNNs_aftercorrection(object.list = corrected.object, assay = rep("integrated", length(merge.list)),
                                                                             batch.cluster.labels = batch.cluster.labels, matched.clusters = matched.clusters, reference = NULL,
                                                                             anchor.features = anchor.features, scale = scale, normalization.method = "LogNormalize",
                                                                             reduction = reduction, l2.norm = l2.norm, dims = dims, k.anchor = k.anchor,
                                                                             k.filter = k.filter, k.score = k.score, max.features = max.features,
                                                                             nn.method = nn.method, eps = eps, verbose = verbose)
          object.anchors_k <- object.anchors_1round
          object.anchors_k@anchors <- object.anchors_after_correction@anchors

          k.weight = k.weight_input
          if (nrow(object.anchors_after_correction@anchors)/2 <= k.weight){
             k.weight = nrow(object.anchors_after_correction@anchors)/2-1
          }
          message("Integrating data...")
          object.combined_k <- IntegrateData(anchorset = object.anchors_k, dims = dims, k.weight = k.weight, sd.weight = sd.weight, verbose = FALSE)
          DefaultAssay(object.combined_k) <- "integrated"

          # Run the standard workflow for visualization and clustering
          object.combined_k <- ScaleData(object.combined_k, verbose = FALSE)
          object.combined_k <- RunPCA(object.combined_k, npcs = npcs, verbose = FALSE)

          message("UMAP embedding ...")
          object.combined_k <- RunUMAP(object.combined_k, reduction = "pca", dims = 1:10, verbose = FALSE)

          ### F measure for 1st round results
          data = data.frame(object.combined_k@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,])
          data$Celltype = factor(unlist(batch.cluster.labels)[object.combined_1round@meta.data$cell.anno %in% matched.clusters])
          data$batch = factor(batch)

          if (length(levels(data$Celltype))>1){
             manova_k = manova(object.combined_k@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,] ~ Celltype+batch, data = data)
          } else{
             manova_k = manova(object.combined_k@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,] ~ batch, data = data)
          }
          temp_k = summary(manova_k)$stats[1:2,3]
          Fbatch_k = as.data.frame(temp_k)["batch",]

          if (Fbatch_k > Fbatch_desired){
             k = k-1
             break
          } else{
             object.combined_desired = object.combined_k
             Fbatch_desired = Fbatch_k
          }
        }
        object.combined_desired <- AddMetaData(object = object.combined_desired, metadata = unlist(batch.cluster.labels), col.name = "cell.anno")

        message("Finished! Batch effect correcting was performed for ", k, " iterations.")
     } else if(strategy == "Long.run"){
        object.combined_k = object.combined_1round
        Fbatch_desired = Fbatch_1round
        Fbatch_record = c(Fbatch_1round)
        for (k in 2:step){
          message(k, " round of batch effect correcting ...")
          corrected.object = SplitObject(object.combined_k, split.by = "batch_id")
          object.anchors_after_correction <- iSMNN_FindIntegrationAnchors_aftercorrection(object.list = corrected.object, assay = rep("integrated", length(merge.list)),
                                                                                          batch.cluster.labels = batch.cluster.labels, matched.clusters = matched.clusters, reference = NULL,
                                                                                          anchor.features = anchor.features, scale = scale, normalization.method = "LogNormalize",
                                                                                          reduction = reduction, l2.norm = l2.norm, dims = dims, k.anchor = k.anchor,
                                                                                          k.filter = k.filter, k.score = k.score, max.features = max.features,
                                                                                          nn.method = nn.method, eps = eps, verbose = verbose)
         object.anchors_k <- object.anchors_1round
         object.anchors_k@anchors <- object.anchors_after_correction@anchors

         k.weight = k.weight_input
         if (nrow(object.anchors_after_correction@anchors)/2 <= k.weight){
            k.weight = nrow(object.anchors_after_correction@anchors)/2-1
         }
         message("Integrating data...")
         object.combined_k <- IntegrateData(anchorset = object.anchors_k, dims = dims, k.weight = k.weight, sd.weight = sd.weight, verbose = FALSE)
         DefaultAssay(object.combined_k) <- "integrated"

         # Run the standard workflow for visualization and clustering
         object.combined_k <- ScaleData(object.combined_k, verbose = FALSE)
         object.combined_k <- RunPCA(object.combined_k, npcs = npcs, verbose = FALSE)

         message("UMAP embedding ...")
         object.combined_k <- RunUMAP(object.combined_k, reduction = "pca", dims = 1:10, verbose = FALSE)

         ### F measure for 1st round results
         data = data.frame(object.combined_k@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,])
         data$Celltype = factor(unlist(batch.cluster.labels)[object.combined_1round@meta.data$cell.anno %in% matched.clusters])
         data$batch = factor(batch)

         if (length(levels(data$Celltype))>1){
            manova_k = manova(object.combined_k@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,] ~ Celltype+batch, data = data)
         } else{
            manova_k = manova(object.combined_k@reductions$umap@cell.embeddings[object.combined_1round@meta.data$cell.anno %in% matched.clusters,] ~ batch, data = data)
         }
         temp_k = summary(manova_k)$stats[1:2,3]
         Fbatch_k = as.data.frame(temp_k)["batch",]

         Fbatch_record = c(Fbatch_record, Fbatch_k)
         assign(paste0("object.combined_", k, "round"), object.combined_k)
       }

       desired.k <- which(Fbatch_record == min(Fbatch_record))
       object.combined_desired <- get(paste0("object.combined_", desired.k, "round"))
       object.combined_desired <- AddMetaData(object = object.combined_desired, metadata = unlist(batch.cluster.labels), col.name = "cell.anno")

       message("Finished! Batch effect correcting was performed for ", desired.k, " round.")
    }

    return(object.combined_desired)
    
  } else{
    
    return(object.combined_desired)
    
  }
}
