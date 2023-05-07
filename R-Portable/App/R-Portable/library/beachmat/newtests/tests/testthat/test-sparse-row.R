# This tests the LIN block reader.
# library(testthat); library(morebeachtests); source("setup.R"); source("test-sparse-row.R")

set.seed(10000)

test_that("sparse matrix row reads are done correctly", {
    for (mats in list(
            SPAWN(100, 20, mode=0),
            SPAWN(10, 200, mode=1),
            SPAWN(50, 50, mode=2)
        )
    ) {
        reference <- mats[[1]]
        nr <- nrow(reference)
        for (M in mats[-1]) {
            for (perm in list(
                    sample(nr),
                    seq_len(nr),
                    rev(seq_len(nr))
                )
            ) {
                for (j in c(0, 2)) {
                    out <- morebeachtests:::get_sparse_row(M, perm - 1L, j)
                    CHECK_SPARSE_IDENTITY(reference, out, mode=j)
                }
            }
        }
    }
})

test_that("sparse matrix row reads are done correctly with slices", {
   for (mats in list(
            SPAWN(100, 20, mode=0),
            SPAWN(10, 200, mode=1),
            SPAWN(50, 50, mode=2)
        )
    ) {
        reference <- mats[[1]]
        nr <- nrow(reference)
        nc <- ncol(reference)
        for (M in mats[-1]) {

            for (perm in list(
                    sample(nr),
                    seq_len(nr),
                    rev(seq_len(nr))
                )
            ) {
                x1 <- sample(nc, nr, replace=TRUE)
                x2 <- sample(nc, nr, replace=TRUE)
                starts <- pmin(x1, x2)
                ends <- pmax(x1, x2)

                ref <- SLICE_ROWS(reference, perm, starts, ends)
                for (j in c(0, 2)) {
                    mat <- morebeachtests:::get_sparse_row_slice(M, perm - 1L, starts - 1L, ends, j)
                    CHECK_SPARSE_IDENTITY(ref, mat, mode=j)
                }
            }
        }
    }
})

test_that("sparse matrix row reads are correctly bounded", {
    mats <- SPAWN(100, 20, mode=2)

    for (M in mats[-1]) {
        expect_error(morebeachtests:::get_sparse_row(M, 1000L, mode=2L),
            "row index out of range")

        expect_error(morebeachtests:::get_sparse_row_slice(M, 1000L, 0, 1, mode=2L),
            "row index out of range")

        expect_error(morebeachtests:::get_sparse_row_slice(M, 0L, -1, 1, mode=2L),
            "column start index is greater than column end index")

        expect_error(morebeachtests:::get_sparse_row_slice(M, 0L, 0, 1000, mode=2L),
            "column end index out of range")
    }
})
