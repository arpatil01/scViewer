## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  library(devtools)
#  install_github("zdebruine/RcppML")

## -----------------------------------------------------------------------------
library(RcppML)
library(Matrix)

## -----------------------------------------------------------------------------
# construct a system of equations
X <- matrix(rnorm(2000),100,20)
btrue <- runif(20)
y <- X %*% btrue + rnorm(100)
a <- crossprod(X)
b <- crossprod(X, y)

# solve the system of equations
x <- RcppML::nnls(a, b)

# use only coordinate descent
x <- RcppML::nnls(a, b, fast_nnls = FALSE, cd_maxit = 1000, cd_tol = 1e-8)

## -----------------------------------------------------------------------------
# simulate a sparse matrix
A <- rsparsematrix(1000, 100, 0.1)

# simulate a linear factor model
w <- matrix(runif(1000 * 10), 1000, 10)

# project the model
h <- RcppML::project(A, w)

## -----------------------------------------------------------------------------
A <- rsparsematrix(100, 100, 0.1)
model <- RcppML::nmf(A, 10, verbose = F)

w <- model$w
d <- model$d
h <- model$h
model_tolerance <- tail(model$tol, 1)

## -----------------------------------------------------------------------------
A_sym <- as(crossprod(A), "dgCMatrix")

model <- RcppML::nmf(A_sym, 10, verbose = F)

## -----------------------------------------------------------------------------
RcppML::mse(A_sym, model$w, model$d, model$h)

