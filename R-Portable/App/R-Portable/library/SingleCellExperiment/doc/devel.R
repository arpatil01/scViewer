## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)

## ----construct----------------------------------------------------------------
library(SingleCellExperiment)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(assays = list(counts = counts))
sce

## -----------------------------------------------------------------------------
# Function in package A:
AsetX <- function(sce) {
    int_colData(sce)$X <- runif(ncol(sce))
    sce
}

# Function in package B:
BsetX <- function(sce) {
    int_colData(sce)$X <- sample(LETTERS, ncol(sce), replace=TRUE)
    sce
}

## -----------------------------------------------------------------------------
sce2 <- AsetX(sce)
int_colData(sce2)$X
sce2 <- BsetX(sce2)
int_colData(sce2)$X

## -----------------------------------------------------------------------------
AsetX_better <- function(sce) {
    int_colData(sce)$A <- DataFrame(X=runif(ncol(sce)))
    sce
}

BsetX_better <- function(sce) {
    choice <- sample(LETTERS, ncol(sce), replace=TRUE)
    int_colData(sce)$B <- DataFrame(X=choice)
    sce
}

sce2 <- AsetX_better(sce)
sce2 <- BsetX_better(sce2)
int_colData(sce2)$A$X 
int_colData(sce2)$B$X 

## -----------------------------------------------------------------------------
AsetY_better <- function(sce) {
    int_elementMetadata(sce)$A <- DataFrame(Y=runif(nrow(sce)))
    sce
}

BsetY_better <- function(sce) {
    choice <- sample(LETTERS, nrow(sce), replace=TRUE)
    int_elementMetadata(sce)$B <- DataFrame(Y=choice)
    sce
}

sce2 <- AsetY_better(sce)
sce2 <- BsetY_better(sce2)
int_elementMetadata(sce2)$A$Y 
int_elementMetadata(sce2)$B$Y

## -----------------------------------------------------------------------------
AsetZ_better <- function(sce) {
    int_metadata(sce)$A <- list(Z = "Aaron")
    sce
}

BsetZ_better <- function(sce) {
    int_metadata(sce)$B <- list(Z = "Davide")
    sce
}

sce2 <- AsetZ_better(sce)
sce2 <- BsetZ_better(sce2)
int_metadata(sce2)$A$Z
int_metadata(sce2)$B$Z

## -----------------------------------------------------------------------------
sessionInfo()

