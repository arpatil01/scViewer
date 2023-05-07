## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
library(BiocNeighbors)

## -----------------------------------------------------------------------------
nobs <- 10000
ndim <- 20
data <- matrix(runif(nobs*ndim), ncol=ndim)

## -----------------------------------------------------------------------------
fout <- findNeighbors(data, threshold=1)
head(fout$index)
head(fout$distance)

## -----------------------------------------------------------------------------
fout$index[[3]]

## -----------------------------------------------------------------------------
fout$distance[[3]]

## -----------------------------------------------------------------------------
nquery <- 1000
ndim <- 20
query <- matrix(runif(nquery*ndim), ncol=ndim)

## -----------------------------------------------------------------------------
qout <- queryNeighbors(data, query, threshold=1)
length(qout$index)

## -----------------------------------------------------------------------------
pre <- buildIndex(data, BNPARAM=KmknnParam())
fout.pre <- findNeighbors(BNINDEX=pre, threshold=1)
qout.pre <- queryNeighbors(BNINDEX=pre, query=query, threshold=1)

## -----------------------------------------------------------------------------
sessionInfo()

