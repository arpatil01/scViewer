# Checks that the promotion machinery works correctly.
# library(beachmat); library(testthat); source('setup.R'); source('test-promotion.R')

test_that("promotion to sparse objects works correctly", {
    mats <- SPAWN(100, 50, mode=2)
    reference <- sum(mats[[1]])
    for (M in mats[-1]) {
        expect_equal(morebeachtests:::test_promotion(M), reference)
    }

    expect_equal(morebeachtests:::test_promotion(mats[[1]]), reference + 1)
})
