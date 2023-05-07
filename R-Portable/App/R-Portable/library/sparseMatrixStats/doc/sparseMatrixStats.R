## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(sparseMatrixStats)
# Matrix defines the sparse Matrix class
# dgCMatrix that we will use
library(Matrix)
# For reproducibility
set.seed(1)

## -----------------------------------------------------------------------------
customer_ids <- seq_len(100)
item_ids <-  seq_len(30)
n_transactions <- 1000
customer <- sample(customer_ids, size = n_transactions, replace = TRUE,
                    prob = runif(100))
item <- sample(item_ids, size = n_transactions, replace = TRUE,
               prob = runif(30))

tmp <- table(paste0(customer, "-", item))
tmp2 <- strsplit(names(tmp), "-")
purchase_table <- data.frame(
  customer = as.numeric(sapply(tmp2, function(x) x[1])),
  item = as.numeric(sapply(tmp2, function(x) x[2])),
  n = as.numeric(tmp)
)

head(purchase_table, n = 10)

## -----------------------------------------------------------------------------
purchase_matrix <- sparseMatrix(purchase_table$customer, purchase_table$item, 
                x = purchase_table$n,
                dims = c(100, 30),
                dimnames = list(customer = paste0("Customer_", customer_ids),
                                item = paste0("Item_", item_ids)))
purchase_matrix[1:10, 1:15]

## -----------------------------------------------------------------------------
# How often was each item bough in total?
colSums2(purchase_matrix)

# What is the range of number of items each 
# customer bought?
head(rowRanges(purchase_matrix))

# What is the variance in the number of items
# each customer bought?
head(rowVars(purchase_matrix))

# How many items did a customer not buy at all, one time, 2 times,
# or exactly 4 times?
head(rowTabulates(purchase_matrix, values = c(0, 1, 2, 4)))

## -----------------------------------------------------------------------------
mat <- matrix(0, nrow=10, ncol=6)
mat[sample(seq_len(60), 4)] <- 1:4
# Convert dense matrix to sparse matrix
sparse_mat <- as(mat, "dgCMatrix")
sparse_mat

## -----------------------------------------------------------------------------
apply(mat, 2, var)

## -----------------------------------------------------------------------------
matrixStats::colVars(mat)

## -----------------------------------------------------------------------------
sparseMatrixStats::colVars(sparse_mat)

## -----------------------------------------------------------------------------
big_mat <- matrix(0, nrow=1e4, ncol=50)
big_mat[sample(seq_len(1e4 * 50), 5000)] <- rnorm(5000)
# Convert dense matrix to sparse matrix
big_sparse_mat <- as(big_mat, "dgCMatrix")

## -----------------------------------------------------------------------------
bench::mark(
  sparseMatrixStats=sparseMatrixStats::colVars(big_sparse_mat),
  matrixStats=matrixStats::colVars(big_mat),
  apply=apply(big_mat, 2, var)
)

## -----------------------------------------------------------------------------
sessionInfo()

