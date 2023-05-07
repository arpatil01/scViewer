## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## -----------------------------------------------------------------------------
mat <- matrix(rnorm(10000), ncol=10)

smat1 <- scale(mat)
head(smat1)

library(DelayedArray)
smat2 <- scale(DelayedArray(mat))
head(smat2)

library(ScaledMatrix)
smat3 <- ScaledMatrix(mat, center=TRUE, scale=TRUE)
head(smat3)

## -----------------------------------------------------------------------------
library(Matrix)
mat <- rsparsematrix(20000, 10000, density=0.01)
smat <- ScaledMatrix(mat, center=TRUE, scale=TRUE)

blob <- matrix(runif(ncol(mat) * 5), ncol=5)
system.time(out <- smat %*% blob)

# The slower way with block processing.
da <- scale(DelayedArray(mat))
system.time(out2 <- da %*% blob)

## -----------------------------------------------------------------------------
library(BiocSingular)
set.seed(1000)
system.time(pcs <- runSVD(smat, k=10, BSPARAM=IrlbaParam()))

## -----------------------------------------------------------------------------
system.time(rowSums(smat))
system.time(rowSums(da))

## -----------------------------------------------------------------------------
smat[,1:5]
t(smat)
rownames(smat) <- paste0("GENE_", 1:20000)
smat

## -----------------------------------------------------------------------------
smat + 1

## -----------------------------------------------------------------------------
set.seed(1000)
mat <- matrix(rnorm(1000000), ncol=100000) 
big.mat <- mat + 1e12

# The 'correct' value, unaffected by numerical precision.
ref <- rowMeans(scale(mat))
head(ref)

# The value from scale'ing a DelayedArray.
library(DelayedArray)
smat2 <- scale(DelayedArray(big.mat))
head(rowMeans(smat2))

# The value from a ScaledMatrix.
library(ScaledMatrix)
smat3 <- ScaledMatrix(big.mat, center=TRUE, scale=TRUE)
head(rowMeans(smat3))

## -----------------------------------------------------------------------------
sessionInfo()

