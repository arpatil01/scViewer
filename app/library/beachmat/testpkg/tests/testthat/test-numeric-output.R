# This tests the ability of the API to properly output numeric matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-numeric-output.R")

sFUN <- numeric_sFUN
dFUN <- numeric_dFUN
csFUN <- numeric_csFUN
tsFUN <- numeric_tsFUN

#######################################################

set.seed(12346)
test_that("Simple numeric matrix output is okay", {
    check_write_all(sFUN, mode="numeric")
    check_write_all(sFUN, nr=5, nc=30, mode="numeric")
    check_write_all(sFUN, nr=30, nc=5, mode="numeric")

    check_write_slice(sFUN, mode="numeric")
    check_write_slice(sFUN, nr=5, nc=30, mode="numeric")
    check_write_slice(sFUN, nr=30, nc=5, mode="numeric")

    check_write_varslice(sFUN, mode="numeric")
    check_write_varslice(sFUN, nr=5, nc=30, mode="numeric")
    check_write_varslice(sFUN, nr=30, nc=5, mode="numeric")

    check_write_indexed(sFUN, mode="numeric")
    check_write_indexed(sFUN, nr=5, nc=30, mode="numeric")
    check_write_indexed(sFUN, nr=30, nc=5, mode="numeric")

    check_write_type(sFUN, mode="numeric")
    check_write_errors(sFUN, mode="numeric")

    check_write_all(sFUN, nr=0, nc=0, mode="numeric")
    check_write_all(sFUN, nr=10, nc=0, mode="numeric")
    check_write_all(sFUN, nr=0, nc=10, mode="numeric")
})

#######################################################

set.seed(23457)
test_that("sparse numeric matrix output is okay", {
    expect_s4_class(csFUN(), "dgCMatrix")

    check_write_all(csFUN, mode="numeric")
    check_write_all(csFUN, nr=5, nc=30, mode="numeric")
    check_write_all(csFUN, nr=30, nc=5, mode="numeric")

    check_write_slice(csFUN, mode="numeric")
    check_write_slice(csFUN, nr=5, nc=30, mode="numeric")
    check_write_slice(csFUN, nr=30, nc=5, mode="numeric")

    check_write_varslice(csFUN, mode="numeric")
    check_write_varslice(csFUN, nr=5, nc=30, mode="numeric")
    check_write_varslice(csFUN, nr=30, nc=5, mode="numeric")

    check_write_indexed(csFUN, mode="numeric")
    check_write_indexed(csFUN, nr=5, nc=30, mode="numeric")
    check_write_indexed(csFUN, nr=30, nc=5, mode="numeric")

    check_write_type(csFUN, mode="numeric")
    check_write_errors(csFUN, mode="numeric")

    check_write_all(csFUN, nr=0, nc=0, mode="numeric")
    check_write_all(csFUN, nr=10, nc=0, mode="numeric")
    check_write_all(csFUN, nr=0, nc=10, mode="numeric")
})

#######################################################

test_that("Numeric matrix mode choices are okay", {
    check_write_class(sFUN(), "matrix")
    check_write_class(csFUN(), "dgCMatrix")
    check_write_class(tsFUN(), "dgTMatrix")
    check_write_class(dFUN(), "dgeMatrix")
})
