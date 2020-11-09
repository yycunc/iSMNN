#!/usr/bin/env Rscript

#' A list of two expression matrices for two batches. The first batch contains 400 cells of
#' three cell types, fibroblasts, macrophages and endothelial cells. And the second batches
#' has 500 cells of the same three cell types.
#'
#' @docType data
#' @usage data("data_iSMNN")
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
#' batch.cluster.labels <- unifiedClusterLabelling(data_SMNN$batch1.mat, data_iSMNN$batch2.mat, features.use = markers,
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
"data_iSMNN"
