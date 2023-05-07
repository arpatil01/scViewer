## ---- echo=FALSE, results="hide"----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE)
library(BiocStyle)
set.seed(10918)

## -----------------------------------------------------------------------------
library(scRNAseq)
sce <- ZeiselBrainData()

library(scuttle)
sce <- quickPerCellQC(sce, subsets=list(Mito=grep("mt-", rownames(sce))),
    sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent")) 

sce

## ---- echo=FALSE--------------------------------------------------------------
# Make the damn thing sparse for speed.
counts(sce) <- as(counts(sce), "dgCMatrix")

## -----------------------------------------------------------------------------
summary(librarySizeFactors(sce))

## -----------------------------------------------------------------------------
summary(geometricSizeFactors(sce))

## -----------------------------------------------------------------------------
summary(medianSizeFactors(sce))

## -----------------------------------------------------------------------------
sizeFactors(sce) <- librarySizeFactors(sce)

## -----------------------------------------------------------------------------
sce <- computeLibraryFactors(sce)
summary(sizeFactors(sce))

## -----------------------------------------------------------------------------
library(scran)
clusters <- quickCluster(sce)

sce <- computePooledFactors(sce, clusters=clusters)
summary(sizeFactors(sce))

## -----------------------------------------------------------------------------
sce <- computePooledFactors(sce, clusters=sce$level1class)
summary(sizeFactors(sce))

## -----------------------------------------------------------------------------
sce2 <- computeSpikeFactors(sce, "ERCC")
summary(sizeFactors(sce2))

## -----------------------------------------------------------------------------
sce <- logNormCounts(sce)
assayNames(sce)

## -----------------------------------------------------------------------------
assay(sce, "normed") <- normalizeCounts(sce, log=FALSE,
    size.factors=runif(ncol(sce)), pseudo.count=1.5)

## -----------------------------------------------------------------------------
assay(sce, "cpm") <- calculateCPM(sce)

## -----------------------------------------------------------------------------
sessionInfo()

