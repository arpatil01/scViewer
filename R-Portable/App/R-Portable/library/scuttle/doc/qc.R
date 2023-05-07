## ---- echo=FALSE, results="hide"----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE)
library(BiocStyle)
set.seed(10918)

## -----------------------------------------------------------------------------
library(scRNAseq)
sce <- ZeiselBrainData()
sce

## -----------------------------------------------------------------------------
library(scuttle)
is.mito <- grep("mt-", rownames(sce))
per.cell <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
summary(per.cell$sum)
summary(per.cell$detected)

## -----------------------------------------------------------------------------
summary(per.cell$subsets_Mito_percent)
summary(per.cell$altexps_ERCC_percent)

## -----------------------------------------------------------------------------
colData(sce) <- cbind(colData(sce), per.cell)

## -----------------------------------------------------------------------------
sce2 <- addPerCellQCMetrics(sce, subsets=list(Mito=is.mito))
colnames(colData(sce2))

## -----------------------------------------------------------------------------
low.total <- isOutlier(per.cell$sum, type="lower", log=TRUE)
summary(low.total)

## -----------------------------------------------------------------------------
attr(low.total, "threshold")

## -----------------------------------------------------------------------------
low.total.batched <- isOutlier(per.cell$sum, type="lower", log=TRUE, batch=sce$tissue)
summary(low.total.batched)

## -----------------------------------------------------------------------------
# An example with just mitochondrial filters.
qc.stats <- perCellQCFilters(per.cell, sub.fields="subsets_Mito_percent") 
colSums(as.matrix(qc.stats))

# Another example with mitochondrial + spike-in filters.
qc.stats2 <- perCellQCFilters(per.cell, 
    sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent")) 
colSums(as.matrix(qc.stats2))

## -----------------------------------------------------------------------------
filtered <- quickPerCellQC(sce, subsets=list(Mito=is.mito), sub.fields="subsets_Mito_percent")
filtered

## -----------------------------------------------------------------------------
# Pretending that the first 10 cells are empty wells, for demonstration.
per.feat <- perFeatureQCMetrics(sce, subsets=list(Empty=1:10))
summary(per.feat$mean)
summary(per.feat$detected)
summary(per.feat$subsets_Empty_ratio)

## -----------------------------------------------------------------------------
ave <- calculateAverage(sce)
summary(ave)

## -----------------------------------------------------------------------------
sessionInfo()

