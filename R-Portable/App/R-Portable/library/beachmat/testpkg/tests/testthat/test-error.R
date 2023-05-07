# This tests the error-generating machinery throughout the package.
# library(beachtest); library(testthat); source("test-error.R")

expect_fixed_error <- function(...) { expect_error(..., fixed=TRUE) }

# Simple matrix doesn't need checking; correct dimensions is enforced,
# and can't even be changed by directly modifying the attributes.

a <- matrix(100, 10, 10)

test_that("Simple matrix errors thrown", {
    expect_fixed_error(attributes(a)$dim <- c(-1L, 10L))
    attributes(a)$dim <- c(10, 10)
    expect_type(attributes(a)$dim, "integer") # type does not change.
    expect_fixed_error(attributes(a)$dim <- c(5L, 5L))
})

# Dense matrix

library(Matrix)
A <- Matrix(1:50, 5, 10)

test_that("Dense matrix errors thrown", {
    wrong <- A
    expect_fixed_error(storage.mode(wrong@Dim) <- "double")
    wrong@Dim <- c(5L, 5L)
    expect_fixed_error(check_read_class(wrong, mode="numeric", "dense"), "length of 'x' in a dgeMatrix object is inconsistent with its dimensions")
})

# Csparse matrix

set.seed(234234)
A <- rsparsematrix(10, 20, 0.5)
    
test_that("Sparse matrix errors thrown", {
    wrong <- A
    expect_fixed_error(storage.mode(wrong@Dim) <- "double")
    wrong@Dim <- c(5L, 10L)
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "length of 'p' slot in a dgCMatrix object should be equal to 'ncol+1'")
    
    wrong <- A
    expect_fixed_error(storage.mode(wrong@p) <- "double")
    wrong@p[1] <- -1L
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "first element of 'p' in a dgCMatrix object should be 0")

    wrong <- A
    wrong@p[ncol(A)+1] <- -1L
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "last element of 'p' in a dgCMatrix object should be 'length(x)'")

    wrong <- A
    wrong@p[2] <- wrong@p[2]+100L
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "'p' slot in a dgCMatrix object should be sorted")
    wrong <- A
    wrong@p <- wrong@p[1]
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "length of 'p' slot in a dgCMatrix object should be equal to 'ncol+1'")
    
    wrong <- A
    expect_fixed_error(storage.mode(wrong@i) <- "double")
    wrong@i <- rev(wrong@i)
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "'i' in each column of a dgCMatrix object should be sorted")

    wrong <- A
    wrong@i <- wrong@i[1]
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "'x' and 'i' slots in a dgCMatrix object should have the same length")

    wrong <- A
    wrong@i <- wrong@i*100L
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "'i' slot in a dgCMatrix object should contain elements in [0, nrow)")
    
    wrong <- A
    wrong@x <- wrong@x[1]
    expect_fixed_error(check_read_class(wrong, mode="numeric", "sparse"), "'x' and 'i' slots in a dgCMatrix object should have the same length") 
})

# Numeric checks

test_that("Numeric errors thrown", {
    B <- matrix(100L, 10, 10)
    expect_fixed_error(check_read_class(B, mode="numeric", "simple"), "matrix should be double")
    
    B <- Matrix(1:50, 5, 10)
    storage.mode(B@x) <- "integer"
    expect_fixed_error(check_read_class(B, mode="numeric", "dense"), "'x' slot in a dgeMatrix object should be double")
    
    set.seed(234234)
    B <- rsparsematrix(10, 20, 0.5)
    storage.mode(B@x) <- "integer"
    expect_fixed_error(check_read_class(B, mode="numeric", "sparse"), "'x' slot in a dgCMatrix object should be double")

    expect_fixed_error(check_read_class(data.frame(1:5), mode="numeric", "simple"), "data.frames")
})

# Logical checks

test_that("Logical errors thrown", {
    B <- matrix(100L, 10, 10)
    expect_fixed_error(check_read_class(B, mode="logical", "simple"), "matrix should be logical")
    
    B <- Matrix(c(TRUE, FALSE), 5, 10)
    expect_fixed_error(storage.mode(B@x) <- "integer")
    
    set.seed(234234)
    B <- rsparsematrix(10, 20, 0.5)!=0
    expect_fixed_error(storage.mode(B@x) <- "integer")
    
    B <- as(rsparsematrix(10, 20, 0.5)!=0, "lgTMatrix")
    expect_fixed_error(storage.mode(B@x) <- "integer")

    expect_fixed_error(check_read_class(data.frame(1:5), mode="logical", "simple"), "data.frames")
})

# Integer checks

test_that("Integer errors thrown", {
    B <- matrix(1, 10, 10)
    expect_fixed_error(check_read_class(B, mode="integer", "simple"), "matrix should be integer")

    expect_fixed_error(check_read_class(data.frame(B), mode="integer", "simple"), "data.frames")
})

# Character checks

test_that("Character errors thrown", {
    B <- matrix(1, 10, 10)
    expect_fixed_error(check_read_class(B, mode="character", "simple"), "matrix should be character")

    expect_fixed_error(check_read_class(data.frame(B), mode="character", "simple"), "data.frames")
})
