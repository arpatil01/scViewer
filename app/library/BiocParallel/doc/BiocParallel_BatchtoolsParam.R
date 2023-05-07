## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----setup, echo=FALSE-----------------------------------------------------
suppressPackageStartupMessages({
    library(BiocParallel)
})

## ----intro-----------------------------------------------------------------
library(BiocParallel)

## Pi approximation
piApprox <- function(n) {
    nums <- matrix(runif(2 * n), ncol = 2)
    d <- sqrt(nums[, 1]^2 + nums[, 2]^2)
    4 * mean(d <= 1)
}

piApprox(1000)

## Apply piApprox over
param <- BatchtoolsParam()
result <- bplapply(rep(10e5, 10), piApprox, BPPARAM=param)
mean(unlist(result))

## --------------------------------------------------------------------------
registryargs <- batchtoolsRegistryargs(
    file.dir = "mytempreg",
    work.dir = getwd(),
    packages = character(0L),
    namespaces = character(0L),
    source = character(0L),
    load = character(0L)
)
param <- BatchtoolsParam(registryargs = registryargs)
param

## --------------------------------------------------------------------------
fname <- batchtoolsTemplate("slurm")
cat(readLines(fname), sep="\n")

## ----simple_sge_example, eval=FALSE----------------------------------------
#  library(BiocParallel)
#  
#  ## Pi approximation
#  piApprox <- function(n) {
#      nums <- matrix(runif(2 * n), ncol = 2)
#      d <- sqrt(nums[, 1]^2 + nums[, 2]^2)
#      4 * mean(d <= 1)
#  }
#  
#  template <- system.file(
#      package = "BiocParallel",
#      "unitTests", "test_script", "test-sge-template.tmpl"
#  )
#  param <- BatchtoolsParam(workers=5, cluster="sge", template=template)
#  
#  ## Run parallel job
#  result <- bplapply(rep(10e5, 100), piApprox, BPPARAM=param)

## ----demo_sge, eval=FALSE--------------------------------------------------
#  library(BiocParallel)
#  
#  ## Pi approximation
#  piApprox <- function(n) {
#      nums <- matrix(runif(2 * n), ncol = 2)
#      d <- sqrt(nums[, 1]^2 + nums[, 2]^2)
#      4 * mean(d <= 1)
#  }
#  
#  template <- system.file(
#      package = "BiocParallel",
#      "unitTests", "test_script", "test-sge-template.tmpl"
#  )
#  param <- BatchtoolsParam(workers=5, cluster="sge", template=template)
#  
#  ## start param
#  bpstart(param)
#  
#  ## Display param
#  param
#  
#  ## To show the registered backend
#  bpbackend(param)
#  
#  ## Register the param
#  register(param)
#  
#  ## Check the registered param
#  registered()
#  
#  ## Run parallel job
#  result <- bplapply(rep(10e5, 100), piApprox)
#  
#  bpstop(param)

## ----sessionInfo, results="asis"-------------------------------------------
toLatex(sessionInfo())

