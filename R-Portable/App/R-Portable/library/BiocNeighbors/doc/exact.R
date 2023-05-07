## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
library(BiocNeighbors)

## -----------------------------------------------------------------------------
nobs <- 10000
ndim <- 20
data <- matrix(runif(nobs*ndim), ncol=ndim)

## -----------------------------------------------------------------------------
fout <- findKNN(data, k=10, BNPARAM=KmknnParam())
head(fout$index)
head(fout$distance)

## -----------------------------------------------------------------------------
fout$index[3,]

## -----------------------------------------------------------------------------
fout$distance[3,]

## -----------------------------------------------------------------------------
nquery <- 1000
ndim <- 20
query <- matrix(runif(nquery*ndim), ncol=ndim)

## -----------------------------------------------------------------------------
qout <- queryKNN(data, query, k=5, BNPARAM=KmknnParam())
head(qout$index)
head(qout$distance)

## -----------------------------------------------------------------------------
qout$index[3,]

## -----------------------------------------------------------------------------
qout$distance[3,]

## -----------------------------------------------------------------------------
findKNN(data, k=5, subset=3:5)

## -----------------------------------------------------------------------------
names(findKNN(data, k=2, get.distance=FALSE))

## -----------------------------------------------------------------------------
library(BiocParallel)
out <- findKNN(data, k=10, BPPARAM=MulticoreParam(3))

## -----------------------------------------------------------------------------
pre <- buildIndex(data, BNPARAM=KmknnParam())
out1 <- findKNN(BNINDEX=pre, k=5)
out2 <- queryKNN(BNINDEX=pre, query=query, k=2)

## -----------------------------------------------------------------------------
out.m <- findKNN(data, k=5, BNPARAM=KmknnParam(distance="Manhattan"))

## -----------------------------------------------------------------------------
sessionInfo()

