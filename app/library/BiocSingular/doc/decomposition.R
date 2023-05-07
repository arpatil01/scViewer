## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## ----setup, echo=FALSE, message=FALSE-----------------------------------------
library(BiocSingular)
set.seed(100)

## -----------------------------------------------------------------------------
dummy <- matrix(rnorm(10000), ncol=25)
e.out <- runSVD(dummy, k=10, BSPARAM=ExactParam())
str(e.out)

## -----------------------------------------------------------------------------
set.seed(1000)
i.out <- runSVD(dummy, k=10, BSPARAM=IrlbaParam())

## -----------------------------------------------------------------------------
set.seed(1000)
r.out <- runSVD(dummy, k=10, BSPARAM=RandomParam())

## -----------------------------------------------------------------------------
epam <- ExactParam(fold=10)
epam
cr.out <- runSVD(dummy, k=10, BSPARAM=epam)

## -----------------------------------------------------------------------------
cs.out <- runSVD(dummy, k=10, scale=runif(ncol(dummy)), 
    center=rnorm(ncol(dummy)))

## -----------------------------------------------------------------------------
set.seed(2000)
def.out <- runSVD(dummy, k=10, scale=runif(ncol(dummy)), 
    center=rnorm(ncol(dummy)), BSPARAM=IrlbaParam(deferred=TRUE))

## -----------------------------------------------------------------------------
ipam <- IrlbaParam(tol=1e-8, extra.work=10)
ipam

## -----------------------------------------------------------------------------
set.seed(3000)
i.out <- runSVD(dummy, k=10, BSPARAM=ipam)

## -----------------------------------------------------------------------------
set.seed(4000)
library(BiocParallel)
i.out <- runSVD(dummy, k=10, BSPARAM=ipam, BPPARAM=bpparam())

## -----------------------------------------------------------------------------
pcs.out <- runPCA(dummy, rank=10, BSPARAM=ExactParam())
str(pcs.out)

## -----------------------------------------------------------------------------
head(pcs.out$x)

## -----------------------------------------------------------------------------
sessionInfo()

