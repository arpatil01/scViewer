## ---- include = FALSE, setup--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment = "#>", collapse = TRUE,
                      message = FALSE)
library(BiocStyle)

## ----data_sim, message = FALSE------------------------------------------------
library(DelayedArray)

x <- do.call(cbind, lapply(1:20, function(j) {
  rpois(n = 10000, lambda = sample(20:40, 10000, replace = TRUE))
}))
colnames(x) <- paste0("S", 1:20)
x <- realize(x, "HDF5Array")
x

## ----apply--------------------------------------------------------------------
system.time(row_sds <- apply(x, 1, sd))
head(row_sds)

## ----matrixStats, error = TRUE------------------------------------------------
matrixStats::rowSds(x)

## ----realization--------------------------------------------------------------
system.time(row_sds <- matrixStats::rowSds(as.matrix(x)))
head(row_sds)

## ----DelayedMatrixStats-------------------------------------------------------
library(DelayedMatrixStats)

system.time(row_sds <- rowSds(x))
head(row_sds)

## ----API, echo = FALSE--------------------------------------------------------
matrixStats <- sort(
  c("colsum", "rowsum", grep("^(col|row)", 
                             getNamespaceExports("matrixStats"), 
                             value = TRUE)))
sparseMatrixStats <- getNamespaceExports("sparseMatrixStats")
DelayedMatrixStats <- getNamespaceExports("DelayedMatrixStats")
DelayedArray <- getNamespaceExports("DelayedArray")

api_df <- data.frame(
  Method = paste0("`", matrixStats, "()`"),
  `Block processing` = ifelse(
    matrixStats %in% DelayedMatrixStats,
    "✔",
    ifelse(matrixStats %in% c(DelayedArray, sparseMatrixStats), "☑️", "❌")),
  `_base::matrix_ optimized` = 
    ifelse(sapply(matrixStats, existsMethod, signature = "matrix_OR_array_OR_table_OR_numeric"), 
           "✔", 
           "❌"),
  `_Matrix::dgCMatrix_ optimized` = 
    ifelse(sapply(matrixStats, existsMethod, signature = "xgCMatrix") | sapply(matrixStats, existsMethod, signature = "dgCMatrix"), 
           "✔", 
           "❌"),
  `_Matrix::lgCMatrix_ optimized` = 
    ifelse(sapply(matrixStats, existsMethod, signature = "xgCMatrix") | sapply(matrixStats, existsMethod, signature = "lgCMatrix"), 
           "✔", 
           "❌"),
  `_DelayedArray::RleArray_ (_SolidRleArraySeed_) optimized` = 
    ifelse(sapply(matrixStats, existsMethod, signature = "SolidRleArraySeed"),
           "✔", 
           "❌"),
  `_DelayedArray::RleArray_  (_ChunkedRleArraySeed_) optimized` = 
    ifelse(sapply(matrixStats, existsMethod, signature = "ChunkedRleArraySeed"),
           "✔", 
           "❌"),
  `_HDF5Array::HDF5Matrix_ optimized` = 
    ifelse(sapply(matrixStats, existsMethod, signature = "HDF5ArraySeed"),
           "✔", 
           "❌"),
  `_base::data.frame_ optimized` = 
    ifelse(sapply(matrixStats, existsMethod, signature = "data.frame"),
           "✔", 
           "❌"),
  `_S4Vectors::DataFrame_ optimized` =
    ifelse(sapply(matrixStats, existsMethod, signature = "DataFrame"),
           "✔", 
           "❌"), 
  check.names = FALSE)
knitr::kable(api_df, row.names = FALSE)

## ----benchmarking, message = FALSE, echo = TRUE, error = TRUE-----------------
library(DelayedMatrixStats)
library(sparseMatrixStats)
library(microbenchmark)
library(profmem)

set.seed(666)

# -----------------------------------------------------------------------------
# Dense with values in (0, 1)
# Fast, memory-efficient column sums of DelayedMatrix with ordinary matrix seed
#

# Generate some data
dense_matrix <- matrix(runif(20000 * 600), 
                       nrow = 20000,
                       ncol = 600)

# Benchmark
dm_matrix <- DelayedArray(dense_matrix)
class(seed(dm_matrix))
dm_matrix
microbenchmark(
  block_processing = colSums2(dm_matrix, force_block_processing = TRUE),
  seed_aware = colSums2(dm_matrix),
  times = 10)
total(profmem(colSums2(dm_matrix, force_block_processing = TRUE)))
total(profmem(colSums2(dm_matrix)))

# -----------------------------------------------------------------------------
# Sparse (60% zero) with values in (0, 1)
# Fast, memory-efficient column sums of DelayedMatrix with ordinary matrix seed
#

# Generate some data
sparse_matrix <- dense_matrix
zero_idx <- sample(length(sparse_matrix), 0.6 * length(sparse_matrix))
sparse_matrix[zero_idx] <- 0

# Benchmark
dm_dgCMatrix <- DelayedArray(Matrix(sparse_matrix, sparse = TRUE))
class(seed(dm_dgCMatrix))
dm_dgCMatrix
microbenchmark(
  block_processing = colSums2(dm_dgCMatrix, force_block_processing = TRUE),
  seed_aware = colSums2(dm_dgCMatrix),
  times = 10)
total(profmem(colSums2(dm_dgCMatrix, force_block_processing = TRUE)))
total(profmem(colSums2(dm_dgCMatrix)))

# -----------------------------------------------------------------------------
# Dense with values in {0, 100} featuring runs of identical values
# Fast, memory-efficient column sums of DelayedMatrix with Rle-based seed
#

# Generate some data
runs <- rep(sample(100, 500000, replace = TRUE), rpois(500000, 100))
runs <- runs[seq_len(20000 * 600)]
runs_matrix <- matrix(runs, 
                      nrow = 20000,
                      ncol = 600)

# Benchmark
dm_rle <- RleArray(Rle(runs),
                   dim = c(20000, 600))
class(seed(dm_rle))
dm_rle
microbenchmark(
  block_processing = colSums2(dm_rle, force_block_processing = TRUE),
  seed_aware = colSums2(dm_rle),
  times = 10)
total(profmem(colSums2(dm_rle, force_block_processing = TRUE)))
total(profmem(colSums2(dm_rle)))

## ----sin----------------------------------------------------------------------
system.time(sin_dm_matrix <- sin(dm_matrix))

## ----colSums2_sin-------------------------------------------------------------
all.equal(colSums2(sin_dm_matrix), colSums(sin(as.matrix(dm_matrix))))

