# iSMNN
### iSMNN: Batch Effect Correction for Single-cell RNA-seq data via Iterative Supervised Mutual Nearest Neighbor Refinement

Batch effect correction is an essential step in the integrative analysis of multiple single cell RNA-seq (scRNA-seq) data. One state-of-the-art strategy for batch effect correction is via unsupervised or supervised detection of mutual nearest neighbors (MNNs). However, both two kinds of methods only detect MNNs across batches on the top of uncorrected data, where the large batch effect may affect the MNN search. The number of MNNs is supposed to be small when the batch effect across samples is large, which may lead to less accurate correction results compared to those with a larger amount of MNNs. 

Here, we present iSMNN, an iterative supervised batch effect correction method that performs multiple rounds of MNN refining and batch effect correction. Our benchmarking showed that, with the help of iterative MNN refining, iSMNN has been demonstrated advantages in removing batch effect yet maximally retaining cell type specific biological features, in particular to the scenarios where different batches largely differ.

iSMNN is maintained by Yuchen Yang [yyuchen@email.unc.edu] and Yun Li [yun_li@med.unc.edu].

## News and Updates
Febrary 22, 2021
* Version 1.00 released
  + Small errors are fixed in this version

October 29, 2020
* Version 0.99.0 released
  + First offical release
  

## Brief introduction

The current implementation of iSMNN encompasses two major steps: one optional cluster harmonizing step and the other batch effect correction step. In the first step, the clusters/cell type labels are matched/harmonized across multiple scRNA-seq batches using *unifiedClusterLabelling* function from [**SMNN** package](https://github.com/yycunc/SMNN) (Yang *et al.* 2020). This entire clustering step can be by-passed by feeding iSMNN cell cluster labels. With cell cluster label information, iSMNN iteratively searches mutual nearest neighbors within each harmonized cell type, and performs batch effect correction using the *iSMNN* function.

In this tutorial, we will perform batch effect correction using iSMNN in a toy example containing two batches. The first batch contains 400 cells from three cell types, namely fibroblasts, macrophages and endothelial cells. And the second batches has 500 cells from the same three cell types. Both two batches contain 3000 genes.


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


## The optional step for cluster harmonizing

The function *unifiedClusterLabelling* from [**SMNN** package](https://github.com/yycunc/SMNN) is used to match/harmonize the clusters/cell type labels across multiple scRNA-seq batches. It takes as input raw expression matrices from two or more batches, a list of marker genes and their corresponding cluster labels, and outputs harmonized cluster label for every single cells across all batches.

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
names(batch.cluster.labels[[1]])=colnames(data_iSMNN$batch1.mat)
names(batch.cluster.labels[[2]])=colnames(data_iSMNN$batch2.mat)
```

## Batch effect correction using *iSMNN* function

With harmonized cluster label information for single cells across batches, we implement batch effect correction using iSMNN.

### Construct the input object for batches using Seurat

Input object is first constructed following the instruction of **Seurat v3 package**. See the tutorial from the [Seurat website](https://satijalab.org/seurat/) for details.

```{r perform batch effect correction using iSMNN}
library(Seurat)
merge <- CreateSeuratObject(counts = cbind(data_iSMNN$batch1.mat, data_iSMNN$batch2.mat), project = "merge", min.cells = 0, min.features = 0)
batch_id <- c(rep("batch1", ncol(data_iSMNN$batch1.mat)), rep("batch2", ncol(data_iSMNN$batch2.mat)))
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
                            strategy = "Short.run", iterations = 5, dims = 1:20, npcs = 30, k.filter = 30)
```

***iSMNN*** function will return a Seurat object that contains the batch-corrected expression matrix for batches

```{r output from SMNNcorrect}
# Output after correction for batch
corrected.results

## An object of class Seurat 
## 4624 features across 900 samples within 2 assays 
## Active assay: integrated (1624 features)
##  1 other assay present: RNA
##  2 dimensional reductions calculated: pca, umap
```

## Batch effect correction using SMNN

Our previous developed batch effect correction method SMNN (Yang *et al.*, Briefings in Bioinformatics, 2020)) can be considered as a special case of iSMNN, where only one iteration is execute. Thus, users can also implement SMNN for batch effect correction using ***iSMNN*** function via simply defining the number of iteration as 1 (iterations = 1), for example,
```{r perform batch effect correction using iSMNN}
corrected.results <- iSMNN(object.list = merge.list, batch.cluster.labels = batch.cluster.labels, matched.clusters = c("endothelial cells", "macrophage", "fibroblast"),
                            strategy = "Short.run", iterations = 1, dims = 1:20, npcs = 30, k.filter = 30)
```

## Citation
Yang, Y., Li, G., Xie, Y., Wang, L., Yang, Y., Liu, J., Qian, L., Li., Y. (2020) iSMNN: Batch Effect Correction for Single-cell RNA-seq data via Iterative Supervised Mutual Nearest Neighbor Refinement. *biorxiv*, https://www.biorxiv.org/content/10.1101/2020.11.09.375659v1.

## Credits
Some functions are borrowed from or executed according to the [Seurat v3 package](https://github.com/satijalab/seurat) (Sturat *et al.*, 2019).
