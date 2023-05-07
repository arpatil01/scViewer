## -----------------------------------------------------------------------------
library(BiocParallel)

## ----quick_start FUN----------------------------------------------------------
FUN <- function(x) { round(sqrt(x), 4) }

## ----quick_start registry-----------------------------------------------------
registered()

## ----configure_registry, eval=FALSE-------------------------------------------
#  options(MulticoreParam=MulticoreParam(workers=4))

## ----quickstart_bplapply_default, eval=FALSE----------------------------------
#  bplapply(1:4, FUN)

## ----quickstart_snow----------------------------------------------------------
param <- SnowParam(workers = 2, type = "SOCK")
bplapply(1:4, FUN, BPPARAM = param)

## ----BiocParallelParam_SerialParam--------------------------------------------
serialParam <- SerialParam()
serialParam

## ----BiocParallelParam_MulticoreParam-----------------------------------------
multicoreParam <- MulticoreParam(workers = 8)
multicoreParam

## ----register_registered------------------------------------------------------
registered()

## ----register_bpparam---------------------------------------------------------
bpparam()

## ----register_BatchtoolsParam-------------------------------------------------
default <- registered()
register(BatchtoolsParam(workers = 10), default = TRUE)

## ----register_BatchtoolsParam2------------------------------------------------
names(registered())
bpparam()

## ----register_restore---------------------------------------------------------
for (param in rev(default))
    register(param)

## ----error-vignette, eval=FALSE-----------------------------------------------
#  browseVignettes("BiocParallel")

## ----use_cases_data-----------------------------------------------------------
library(RNAseqData.HNRNPC.bam.chr14)
fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES

## ----forking_gr, message=FALSE------------------------------------------------
library(GenomicAlignments) ## for GenomicRanges and readGAlignments()
gr <- GRanges("chr14", IRanges((1000:3999)*5000, width=1000))

## ----forking_param------------------------------------------------------------
param <- ScanBamParam(which=range(gr))

## ----forking_FUN--------------------------------------------------------------
FUN <- function(fl, param) {
    gal <- readGAlignments(fl, param = param)
    sum(countOverlaps(gr, gal))
}

## ----forking_default_multicore------------------------------------------------
MulticoreParam()

## ----db_problems, eval = FALSE------------------------------------------------
#  library(org.Hs.eg.db)
#  FUN <- function(x, ...) {
#  ...
#  mapIds(org.Hs.eg.db, ...)
#  ...
#  }
#  bplapply(X, FUN, ..., BPPARAM = MulticoreParam())

## ----cluster_FUN--------------------------------------------------------------
FUN <- function(fl, param, gr) {
    suppressPackageStartupMessages({
        library(GenomicAlignments)
    })
    gal <- readGAlignments(fl, param = param)
    sum(countOverlaps(gr, gal))
}

## ----cluster_snow_param-------------------------------------------------------
snow <- SnowParam(workers = 2, type = "SOCK")

## ----cluster_bplapply---------------------------------------------------------
bplapply(fls[1:3], FUN, BPPARAM = snow, param = param, gr = gr)

## ----db_solution_2, eval = FALSE----------------------------------------------
#  register(SnowParam()) # default evaluation
#  bpstart() # start the cluster
#  ...
#  bplapply(X, FUN1, ...)
#  ...
#  bplapply(X, FUN2, ...) # re-use workers
#  ...
#  bpstop()

## ----cluster-MPI-work, eval=FALSE---------------------------------------------
#  library(BiocParallel)
#  library(Rmpi)
#  FUN <- function(i) system("hostname", intern=TRUE)

## ----cluster-MPI, eval=FALSE--------------------------------------------------
#  param <- SnowParam(mpi.universe.size() - 1, "MPI")
#  register(param)

## ----cluster-MPI-do, eval=FALSE-----------------------------------------------
#  xx <- bplapply(1:100, FUN)
#  table(unlist(xx))
#  mpi.quit()

## ----cluster-MPI-bpstart, eval=FALSE------------------------------------------
#  param <- bpstart(SnowParam(mpi.universe.size() - 1, "MPI"))
#  register(param)
#  xx <- bplapply(1:100, FUN)
#  bpstop(param)
#  mpi.quit()

## ----slurm--------------------------------------------------------------------
tmpl <- system.file(package="batchtools", "templates", "slurm-simple.tmpl")
noquote(readLines(tmpl))

## ----cluster-batchtools, eval=FALSE-------------------------------------------
#  ## define work to be done
#  FUN <- function(i) system("hostname", intern=TRUE)
#  library(BiocParallel)
#  
#  ## register SLURM cluster instructions from the template file
#  param <- BatchtoolsParam(workers=5, cluster="slurm", template=tmpl)
#  register(param)
#  
#  ## do work
#  xx <- bplapply(1:100, FUN)
#  table(unlist(xx))

## ----devel-bplapply-----------------------------------------------------------
system.time(x <- bplapply(1:3, function(i) { Sys.sleep(i); i }))

unlist(x)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

