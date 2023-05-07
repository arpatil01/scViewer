# This tests the error-generating machinery throughout the package.
# library(beachtest); library(testthat); source("test-error.R")

check_for_error <- function(x, ...) { 
    expect_error(
        morebeachtests:::get_column(x, seq_len(ncol(x)) - 1L, 2),
        ...,
        fixed=TRUE, class = "std::runtime_error") 
}


set.seed(234234)
A <- Matrix::rsparsematrix(10, 20, 0.5)

test_that("Csparse_reader errors thrown", {
    wrong <- A
    wrong@p[1] <- -1L
    check_for_error(wrong, "first element of 'p' in a dgCMatrix object should be 0")

    wrong <- A
    wrong@p[ncol(A)+1] <- -1L
    check_for_error(wrong, "last element of 'p' in a dgCMatrix object should be 'length(x)'")

    wrong <- A
    wrong@p[2] <- -1L 
    check_for_error(wrong, "'p' slot in a dgCMatrix object should be sorted")

    wrong <- A
    wrong@p <- wrong@p[1]
    check_for_error(wrong, "length of 'p' slot in a dgCMatrix object should be equal to 'ncol+1'")
    
    wrong <- A
    wrong@i <- rev(wrong@i)
    check_for_error(wrong, "'i' in each column of a dgCMatrix object should be sorted")

    wrong <- A
    wrong@i <- wrong@i[1]
    check_for_error(wrong, "'x' and 'i' slots in a dgCMatrix object should have the same length")

    wrong <- A
    wrong@i <- wrong@i*100L
    check_for_error(wrong, "'i' slot in a dgCMatrix object should have entries in [0, nrow)")
    
    wrong <- A
    wrong@x <- wrong@x[1]
    check_for_error(wrong, "'x' and 'i' slots in a dgCMatrix object should have the same length") 
})

library(DelayedArray)
B <- as(A, "SparseArraySeed")
    
test_that("SparseArraySeed errors thrown", {
    wrong <- B
    wrong@nzindex <- wrong@nzindex[,1,drop=FALSE]
    check_for_error(wrong, "SparseArraySeed object should have two columns")

    wrong <- B
    wrong@nzdata <- wrong@nzdata[1]
    check_for_error(wrong, "lengths in a SparseArraySeed object") 

    wrong <- B
    wrong@nzindex <- wrong@nzindex * 2L
    check_for_error(wrong, "out of bounds in a SparseArraySeed object")
})

test_that("read_lin_block errors thrown", {
    check_for_error(matrix("A"), 'not a recognized matrix representation')

    x <- matrix(1)
    expect_error(morebeachtests:::get_sparse_column(x, seq_len(ncol(x)) - 1L, 2), "no 'class'")

    x <- Matrix(1)
    expect_error(morebeachtests:::get_sparse_column(x, seq_len(ncol(x)) - 1L, 2), "not a recognized sparse representation")
})

test_that("sparse writers throw appropriate errors", {
    expect_error(morebeachtests:::test_sparse_writer1(1), "entries in 'store' refer to out-of-range columns")
    expect_error(morebeachtests:::test_sparse_writer1(2), "entries in 'store' refer to out-of-range rows")
    expect_error(morebeachtests:::test_sparse_writer1(3), "entries in 'store' refer to out-of-range rows")
    expect_error(morebeachtests:::test_sparse_writer1(4), "entries in 'store' refer to out-of-range columns")
    expect_error(morebeachtests:::test_sparse_writer1(5), NA)

    expect_error(morebeachtests:::test_sparse_writer2(A, A@x[1]), "inconsistent number of non-zero entries")
    expect_identical(A, morebeachtests:::test_sparse_writer2(A, A@x))
    expect_identical(A*2, morebeachtests:::test_sparse_writer2(A, A@x*2))

    expect_error(morebeachtests:::test_sparse_writer3(), "entries in 'store' are not sorted")
})
