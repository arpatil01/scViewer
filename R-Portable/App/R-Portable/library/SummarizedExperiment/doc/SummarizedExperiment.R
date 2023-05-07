## ----style, echo=FALSE, results='asis'----------------------------------------
BiocStyle::markdown()

## ---- echo=FALSE--------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(data(airway, package="airway"))

## -----------------------------------------------------------------------------
library(SummarizedExperiment)
data(airway, package="airway")
se <- airway
se

## ----assays, eval = FALSE-----------------------------------------------------
#  assays(se)$counts

## ----assays_table, echo = FALSE-----------------------------------------------
knitr::kable(assays(se)$counts[1:10,])

## ----rowRanges----------------------------------------------------------------
rowRanges(se)

## ----colData------------------------------------------------------------------
colData(se)

## ----columnSubset-------------------------------------------------------------
# subset for only those samples treated with dexamethasone
se[, se$dex == "trt"]

## ----metadata-----------------------------------------------------------------
metadata(se)

## ----metadata-formula---------------------------------------------------------
metadata(se)$formula <- counts ~ dex + albut

metadata(se)

## ----constructRSE-------------------------------------------------------------
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])

SummarizedExperiment(assays=list(counts=counts),
                     rowRanges=rowRanges, colData=colData)

## ----constructSE--------------------------------------------------------------
SummarizedExperiment(assays=list(counts=counts), colData=colData)

## ----2d-----------------------------------------------------------------------
# subset the first five transcripts and first three samples
se[1:5, 1:3]

## ----colDataExtraction--------------------------------------------------------
se[, se$cell == "N61311"]

## ----getSet-------------------------------------------------------------------
counts <- matrix(1:15, 5, 3, dimnames=list(LETTERS[1:5], LETTERS[1:3]))

dates <- SummarizedExperiment(assays=list(counts=counts),
                              rowData=DataFrame(month=month.name[1:5], day=1:5))

# Subset all January assays
dates[rowData(dates)$month == "January", ]

## ----assay_assays-------------------------------------------------------------
assays(se)

assays(se)[[1]][1:5, 1:5]

# assay defaults to the first assay if no i is given
assay(se)[1:5, 1:5]

assay(se, 1)[1:5, 1:5]

## ----overlap------------------------------------------------------------------
# Subset for only rows which are in the interval 100,000 to 110,000 of
# chromosome 1
roi <- GRanges(seqnames="1", ranges=100000:1100000)
subsetByOverlaps(se, roi)

## -----------------------------------------------------------------------------
sessionInfo()

