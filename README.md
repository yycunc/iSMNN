# iSMNN
### iSMNN: Batch Effect Correction for Single-cell RNA-seq data via Iterative Supervised Mutual Nearest Neighbor Refinement

Batch effect correction is an essential step in the integrative analysis of multiple single cell RNA-seq (scRNA-seq) data. One state-of-the-art strategy for batch effect correction is via unsupervised or supervised detection of mutual nearest neighbors (MNNs). However, both two kinds of methods only detect MNNs across batches on the top of uncorrected data, where the large batch effect may affect the MNN search. The number of MNNs is supposed to be small when the batch effect across samples is large, which may lead to less accurate correction results compared to those with a larger amount of MNNs. 

Here, we present iSMNN, an iterative supervised batch effect correction method that performs multiple rounds of MNN refining and batch effect correction. Our benchmarking showed that, with the help of iterative MNN refining, iSMNN has been demonstrated advantages in removing batch effect yet maximally retaining cell type specific biological features, in particular to the scenarios where different batches largely differ.

iSMNN is maintained by Yuchen Yang [yyuchen@email.unc.edu].

## News and Updates
October 29, 2020
* Version 0.99.0 released
  + First offical release
  

## Brief introduction

The current implementation of SMNN encompasses two major steps: one optional clustering step and the other batch effect correction step. In the first step, SMNN takes the expression matrix as input, and performs clustering using Seurat v. 3.0 (Butler *et al.*, 2018). Corresponding clusters/cell types are then matched across batches based on marker genes specified by the user. This entire clustering step can be by-passed by feeding SMNN cell cluster labels. With cell cluster label information, SMNN searches mutual nearest neighbors within each cell type, and performs batch effect correction using the *SMNNcorrect* function.

In this tutorial, we will perform batch effect correction using SMNN in a toy example containing two batches. The first batch contains 400 cells from three cell types, namely fibroblasts, macrophages and endothelial cells. And the second batches has 500 cells from the same three cell types. Both two batches contain 3000 genes.


## Installation

iSMNN package can be directly installed from GitHub with:
```{r installation}
install.packages("devtools")

devtools::install_github("yycunc/iSMNN")
```


## Set up the library
```{r init, message=TRUE}
library("iSMNN")
```


## Load the input expression matrix

Once installed, one can use the following command lines to load the data used in our toy example using the following command lines: 
```{r set up for input expression data}
data("data_iSMNN")
dim(data_iSMNN$batch1.mat)
data_iSMNN$batch1.mat[1:5, 1:5]
dim(data_iSMNN$batch2.mat)
```


## The optional step for clustering

The function *unifiedClusterLabelling* from **SMNN** package is used to match/harmonize the clusters/cell type labels across multiple scRNA-seq batches. It takes as input raw expression matrices from two or more batches, a list of marker genes and their corresponding cluster labels, and outputs harmonized cluster label for every single cells across all batches. The SMNN pakcage can be installed from [here](https://github.com/yycunc/SMNN).

### Provide cell-type specific marker gene information

Two pieces of information are needed:
- Marker genes
- Corresponding cell type labels for each marker gene

```{r define the marker genes for cluster matching, warning=FALSE}
# Maker genes
markers <- c("Col1a1", "Pdgfra", "Ptprc", "Pecam1")
# Corresponding cell type labels for each marker gene
cluster.info <- c("fibroblast", "fibroblast", "macrophage", "endothelial cells")
```

### Harmonize cluster labels across batches

```{r, results='hide', fig.show="hide", message=FALSE}
library(SMNN)
batch.cluster.labels <- unifiedClusterLabelling(data_SMNN$batch1.mat, data_iSMNN$batch2.mat, features.use = markers, cluster.labels = cluster.info, min.perc = 0.3)
```

## Batch effect correction using iSMNN function

With harmonized cluster label information for single cells across batches, we implement batch effect correction using iSMNN.

### Construct the input object for batches using Seurat

Input object is first constructed following the instruction of Seurat. See the tutorial from the [Seurat website](https://satijalab.org/seurat/) for details.

```{r perform batch effect correction using iSMNN}
library(Seurat)
merge <- CreateSeuratObject(counts = cbind(data_iSMNN$batch1.mat, data_iSMNN$batch2.mat), project = "merge", min.cells = 0, min.features = 0)
batch_id <- c(rep("batch1", ncol(data_iSMNN$batch1.mat)), rep("batch2", ncol(data_iSMNN$batch2.mat))）
names(batch_id) <- colnames(merge)
merge <- AddMetaData(object = merge, metadata = batch_id, col.name = "batch_id")
merge.list <- SplitObject(merge, split.by = "batch_id")

library(cowplot)
merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
```

### Correct batch effect
```{r perform batch effect correction using iSMNN}
corrected.results <- iSMNN(object.list = merge.list, batch.cluster.labels = batch.cluster.labels, matched.clusters = c("endothelial cells", "macrophage", "fibroblast"),
                           iterations = 5, dims = 1:20, npcs = 30)
```

***iSMNN*** function will return a Seurat object that contains: (1) the batch-corrected expression matrix for each batch; and (2) information regarding mutual nearest neighbors.

In the example below, we treat batch 1 as the reference batch and batch 2 as the batch to be corrected (such that batch 2 will be corrected towards the reference batch 1). Note that the reference batch (i.e., batch 1 in our example) will only applied cosine normalization.

```{r output from SMNNcorrect}
# Output after correction for batch
## Output (1): the batch-corrected expression matrix
corrected.results$corrected[[2]][1:10,1:10]
## Output (2): mutual nearest neighbor information
corrected.results$pairs[[2]]
```

## Citation
Yang, Y., Li, G., Yifang Xie, Li Wang, Jiandong Liu, Li Qian, Li., Y. (2020) iSMNN: Batch Effect Correction for Single-cell RNA-seq data via Iterative Supervised Mutual Nearest Neighbor Refinement. *biorxiv*, 

## Credits
Some functions are borrowed from or executed according to the [Seurat v3 package](https://github.com/satijalab/seurat) (Sturat *et al.*, 2019).
