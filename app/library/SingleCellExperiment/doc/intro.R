## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)

## ----construct----------------------------------------------------------------
library(SingleCellExperiment)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(counts)
sce

## -----------------------------------------------------------------------------
sce <- SingleCellExperiment(list(counts=counts))
sce

## -----------------------------------------------------------------------------
pretend.cell.labels <- sample(letters, ncol(counts), replace=TRUE)
pretend.gene.lengths <- sample(10000, nrow(counts))

sce <- SingleCellExperiment(list(counts=counts),
    colData=DataFrame(label=pretend.cell.labels),
    rowData=DataFrame(length=pretend.gene.lengths),
    metadata=list(study="GSE111111")
)
sce

## ----coerce-------------------------------------------------------------------
se <- SummarizedExperiment(list(counts=counts))
as(se, "SingleCellExperiment")

## -----------------------------------------------------------------------------
dim(assay(sce))
colnames(colData(sce))
colnames(rowData(sce))

## ----fluidigm-----------------------------------------------------------------
library(scRNAseq)
sce <- ReprocessedAllenData("tophat_counts")
sce

## ----subset-------------------------------------------------------------------
counts <- assay(sce, "tophat_counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
assayNames(sce)

## ----pca----------------------------------------------------------------------
pca_data <- prcomp(t(logcounts(sce)), rank=50)

library(Rtsne)
set.seed(5252)
tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)

reducedDims(sce) <- list(PCA=pca_data$x, TSNE=tsne_data$Y)
sce

## -----------------------------------------------------------------------------
reducedDims(sce)
reducedDimNames(sce)
head(reducedDim(sce, "PCA")[,1:2])
head(reducedDim(sce, "TSNE")[,1:2])

## -----------------------------------------------------------------------------
dim(reducedDim(sce, "PCA"))
dim(reducedDim(sce[,1:10], "PCA"))

## -----------------------------------------------------------------------------
counts(sce) <- assay(sce, "tophat_counts")
sce
dim(counts(sce))

## -----------------------------------------------------------------------------
altExp(sce)

## -----------------------------------------------------------------------------
rowData(altExp(sce))$concentration <- runif(nrow(altExp(sce)))
rowData(altExp(sce))
rowData(sce)

## -----------------------------------------------------------------------------
is.riken <- grepl("^[0-9]", rownames(sce))
sce <- splitAltExps(sce, ifelse(is.riken, "RIKEN", "gene"))
altExpNames(sce)

## -----------------------------------------------------------------------------
swapAltExp(sce, "RIKEN", saved="original")

## -----------------------------------------------------------------------------
cell1 <- sample(ncol(sce), 100, replace=TRUE)
cell2 <- sample(ncol(sce), 100, replace=TRUE)
distance <- runif(100)

## -----------------------------------------------------------------------------
colPair(sce, "relationships") <- SelfHits(
    cell1, cell2, nnode=ncol(sce), value=distance)
colPair(sce, "relationships")
class(colPair(sce, asSparse=TRUE))

## -----------------------------------------------------------------------------
sub <- sce[,50:300]
colPair(sub) # grabs the first pairing, if no 'type' is supplied.

## -----------------------------------------------------------------------------
# Making up some size factors and storing them:
sizeFactors(sce) <- 2^rnorm(ncol(sce))
summary(sizeFactors(sce))

# Deleting the size factors:
sizeFactors(sce) <- NULL
sizeFactors(sce)

## -----------------------------------------------------------------------------
# Making up some labels and storing them:
colLabels(sce) <- sample(letters, ncol(sce), replace=TRUE)
table(colLabels(sce))

# Deleting the labels:
colLabels(sce) <- NULL
colLabels(sce)

## -----------------------------------------------------------------------------
# Packs integer or character vectors into the rowData:
rowSubset(sce, "my gene set 1") <- 1:10
which(rowSubset(sce, "my gene set 1"))

# Easy to delete:
rowSubset(sce, "my gene set 1") <- NULL
rowSubset(sce, "my gene set 1")

## -----------------------------------------------------------------------------
sessionInfo()

