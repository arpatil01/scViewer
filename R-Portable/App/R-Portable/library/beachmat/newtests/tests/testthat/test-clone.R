# Checks that the cloning machinery works correctly.
# library(beachmat); library(testthat); source('setup.R'); source('test-clone.R')

test_that("cloning works correctly", {
    mats <- SPAWN(100, 50, mode=2)
    reference <- sum(mats[[1]])
    for (M in mats) {
        expect_equal(morebeachtests:::test_clone(M), reference)
    }
})

test_that("cloning with sparse objects works correctly", {
    mats <- SPAWN(100, 50, mode=2)
    reference <- sum(mats[[1]])
    for (M in mats[-1]) {
        expect_equal(morebeachtests:::test_clone_sparse(M), reference)
    }
})
