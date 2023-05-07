## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)

## -----------------------------------------------------------------------------
library(SingleCellExperiment)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(counts)

altExp(sce, "Spike") <- SingleCellExperiment(matrix(rpois(20, lambda = 5), ncol=10, nrow=2))
altExp(sce, "Protein") <- SingleCellExperiment(matrix(rpois(50, lambda = 100), ncol=10, nrow=5))
altExp(sce, "CRISPR") <- SingleCellExperiment(matrix(rbinom(80, p=0.1, 1), ncol=10, nrow=8))

sce

## -----------------------------------------------------------------------------
totalCount <- function(x, i=1, multiplier=1, subset.row=NULL) {
    mat <- assay(x, i)
    if (!is.null(subset.row)) {
        mat <- mat[subset.row,,drop=FALSE]
    }
    colSums(mat) * multiplier
}

## -----------------------------------------------------------------------------
totals <- applySCE(sce, FUN=totalCount)
totals

## -----------------------------------------------------------------------------
totals.manual <- list( 
    totalCount(sce),
    Spike=totalCount(altExp(sce, "Spike")),
    Protein=totalCount(altExp(sce, "Protein")),
    CRISPR=totalCount(altExp(sce, "CRISPR"))
)
stopifnot(identical(totals, totals.manual))

## -----------------------------------------------------------------------------
totals10.manual <- list( 
    totalCount(sce, multiplier=10),
    Spike=totalCount(altExp(sce, "Spike"), multiplier=10),
    Protein=totalCount(altExp(sce, "Protein"), multiplier=10),
    CRISPR=totalCount(altExp(sce, "CRISPR"), multiplier=10)
)

## -----------------------------------------------------------------------------
totals10.apply <- applySCE(sce, FUN=totalCount, multiplier=10)
stopifnot(identical(totals10.apply, totals10.manual))

## -----------------------------------------------------------------------------
totals10.lapply <- lapply(c(List(sce), altExps(sce)),
    FUN=totalCount, multiplier=10)
stopifnot(identical(totals10.apply, totals10.lapply))

## -----------------------------------------------------------------------------
totals.custom <- applySCE(sce, FUN=totalCount, multiplier=10, 
    ALT.ARGS=list(Spike=list(subset.row=2), Protein=list(subset.row=3:5)))
totals.custom

## -----------------------------------------------------------------------------
head.sce <- applySCE(sce, FUN=head, n=5)
head.sce

## -----------------------------------------------------------------------------
altExp(head.sce)
altExp(head.sce, "Protein")
altExp(head.sce, "CRISPR")

## -----------------------------------------------------------------------------
head.sce.list <- applySCE(sce, FUN=head, n=5, SIMPLIFY=FALSE) 
head.sce.list

## -----------------------------------------------------------------------------
manual.head <- head(sce, n=5)
altExp(manual.head, "Spike") <- head(altExp(sce, "Spike"), n=5)
altExp(manual.head, "Protein") <- head(altExp(sce, "Protein"), n=5)
altExp(manual.head, "CRISPR") <- head(altExp(sce, "CRISPR"), n=5)
manual.head

## -----------------------------------------------------------------------------
sessionInfo()

