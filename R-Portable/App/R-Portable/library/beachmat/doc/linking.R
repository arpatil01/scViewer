## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----headers------------------------------------------------------------------
system.file(package="beachmat", "include")

## -----------------------------------------------------------------------------
column_sums(matrix(rnorm(1000), ncol=10))

## -----------------------------------------------------------------------------
# Integer
column_sums(matrix(rpois(1000, 10), ncol=10))

# Logical
column_sums(matrix(rbinom(1000, 1, 0.5)==1, ncol=10))

# Sparse
column_sums(Matrix::rsparsematrix(100, 10, 0.1))

## -----------------------------------------------------------------------------
column_sums_sparse(Matrix::rsparsematrix(100, 10, 0.1))

# Errors out as an ordinary matrix isn't sparse:
try(column_sums_sparse(matrix(rpois(1000, 10), ncol=10)))

## -----------------------------------------------------------------------------
column_sums_flexible(Matrix::rsparsematrix(100, 10, 0.1))

column_sums_flexible(matrix(rpois(1000, 10), ncol=10))

## -----------------------------------------------------------------------------
generate_sparse_general()

## -----------------------------------------------------------------------------
library(Matrix)
mat <- rsparsematrix(10,10,0.1)
generate_sparse_specific(mat)

## -----------------------------------------------------------------------------
library(beachmat)
blockedColSums <- function(x, ...) {
    out <- colBlockApply(x, FUN=column_sums, ...)
    unlist(out) # combining results across blocks.
}

## -----------------------------------------------------------------------------
# Ordinary matrices:
x1 <- matrix(rnorm(1000), ncol=20)
blockedColSums(x1)

# Works for Matrix classes:
x2 <- Matrix(x1)
blockedColSums(x2)

# Works for DelayedMatrix objects:
library(DelayedArray)
x3 <- DelayedArray(x1)
blockedColSums(x3)

## ---- eval=.Platform$OS.type=="unix"------------------------------------------
#  library(BiocParallel)
#  blockedColSums(x3, BPPARAM=MulticoreParam(2))

## -----------------------------------------------------------------------------
sessionInfo()

