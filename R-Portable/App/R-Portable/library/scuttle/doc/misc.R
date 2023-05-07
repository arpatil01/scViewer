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

## -----------------------------------------------------------------------------
library(scuttle)
agg.sce <- aggregateAcrossCells(sce, ids=sce$level1class)
head(assay(agg.sce))
colData(agg.sce)[,c("ids", "ncells")]

## -----------------------------------------------------------------------------
agg.sce <- aggregateAcrossCells(sce, 
    ids=colData(sce)[,c("level1class", "tissue")])
head(assay(agg.sce))
colData(agg.sce)[,c("level1class", "tissue", "ncells")]

## -----------------------------------------------------------------------------
sce <- logNormCounts(sce)
agg.feat <- sumCountsAcrossFeatures(sce,
    ids=list(GeneSet1=1:10, GeneSet2=11:50, GeneSet3=1:100),
    average=TRUE, exprs_values="logcounts")
agg.feat[,1:10]

## -----------------------------------------------------------------------------
agg.n <- summarizeAssayByGroup(sce, statistics="prop.detected",
    ids=colData(sce)[,c("level1class", "tissue")])
head(assay(agg.n))

## -----------------------------------------------------------------------------
# Mocking up a dataset to demonstrate:
outfile <- tempfile()
write.table(counts(sce)[1:100,], file=outfile, sep="\t", quote=FALSE)

# Reading it in as a sparse matrix:
output <- readSparseCounts(outfile)
class(output)

## -----------------------------------------------------------------------------
# Original row names are Ensembl IDs.
sce.ens <- ZeiselBrainData(ensembl=TRUE)
head(rownames(sce.ens)) 

# Replacing with guaranteed unique and non-missing symbols:
rownames(sce.ens) <- uniquifyFeatureNames(
    rownames(sce.ens), rowData(sce.ens)$originalName
)
head(rownames(sce.ens)) 

## -----------------------------------------------------------------------------
out <- makePerCellDF(sce, features="Tspan12")
colnames(out)

## -----------------------------------------------------------------------------
out2 <- makePerFeatureDF(sce, cells=c("1772063062_D05",
    "1772063061_D01", "1772060240_F02", "1772062114_F05"))
colnames(out2)

## -----------------------------------------------------------------------------
sessionInfo()

