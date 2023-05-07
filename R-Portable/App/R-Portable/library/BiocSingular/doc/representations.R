## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----setup, echo=FALSE, message=FALSE-----------------------------------------
library(BiocSingular)
set.seed(100)

## -----------------------------------------------------------------------------
library(Matrix)
a <- rsparsematrix(10000, 1000, density=0.01)
out <- runPCA(a, rank=5, BSPARAM=IrlbaParam(deferred=TRUE)) # deferring for speed.
recon <- LowRankMatrix(out$rotation, out$x)
recon    

## -----------------------------------------------------------------------------
summary(recon[,1])
summary(recon[2,])

## -----------------------------------------------------------------------------
sessionInfo()

