## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
library(BiocNeighbors)

## -----------------------------------------------------------------------------
nobs <- 10000
ndim <- 20
data <- matrix(runif(nobs*ndim), ncol=ndim)

fout <- findKNN(data, k=10, BNPARAM=AnnoyParam())
head(fout$index)
head(fout$distance)

## -----------------------------------------------------------------------------
nquery <- 1000
ndim <- 20
query <- matrix(runif(nquery*ndim), ncol=ndim)

qout <- queryKNN(data, query, k=5, BNPARAM=AnnoyParam())
head(qout$index)
head(qout$distance)

## -----------------------------------------------------------------------------
pre <- buildIndex(data, BNPARAM=AnnoyParam())
out1 <- findKNN(BNINDEX=pre, k=5)
out2 <- queryKNN(BNINDEX=pre, query=query, k=2)

## -----------------------------------------------------------------------------
AnnoyIndex_path(pre)

## -----------------------------------------------------------------------------
sessionInfo()

