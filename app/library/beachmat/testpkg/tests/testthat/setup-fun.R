# Setting up functions for generating the matrices to use for testing.

library(Matrix)
library(DelayedArray)

#######################################################

integer_sFUN <- function(nr=15, nc=10, lambda=5) { 
    matrix(rpois(nr*nc, lambda=lambda), nr, nc)
}

integer_rFUN <- function(nr=15, nc=10, lambda=1, chunk.ncols=NULL) {
    x <- integer_sFUN(nr, nc, lambda=lambda)
    rle <- Rle(x)
    RleArray(rle, dim(x))
}

#######################################################

logical_sFUN <- function(nr=15, nc=10) {
    matrix(rbinom(nr*nc, 1, 0.5)==0, nr, nc)
}

logical_dFUN <- function(nr=15, nc=10) {
    Matrix(logical_sFUN(nr, nc), sparse=FALSE, doDiag=FALSE)
}

logical_csFUN <- function(nr=15, nc=10, d=0.1) {
    rsparsematrix(nrow=nr, ncol=nc, density=d) != 0
}

logical_tsFUN <- function(...) {
	as(logical_csFUN(...), 'lgTMatrix')
}

#######################################################

numeric_sFUN <- function(nr=15, nc=10) {
    matrix(rnorm(nr*nc), nr, nc)
}

numeric_dFUN <- function(nr=15, nc=10) {
    Matrix(numeric_sFUN(nr, nc), sparse=FALSE, doDiag=FALSE)
}

numeric_csFUN <- function(nr=15, nc=10, d=0.1) {
    rsparsematrix(nrow=nr, ncol=nc, density=d)
}

numeric_tsFUN <- function(...) {
	as(numeric_csFUN(...), 'dgTMatrix')
}

#######################################################

genwords <- function(n = 5000) {
    all.choices <- c(rep("", 4), LETTERS) # to get variable length strings.
    a <- do.call(paste0, replicate(5, sample(all.choices, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

character_sFUN <- function(nr=15, nc=10, lambda=5) {
    matrix(genwords(nr*nc), nr, nc)
}

character_rFUN <- function(nr=15, nc=10, lambda=1, chunk.ncols=NULL) {
    x <- matrix(sample(LETTERS[1:4], nr*nc, replace=TRUE), nr, nc)
    rle <- Rle(x)
    RleArray(rle, dim(x))
}
