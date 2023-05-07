# This tests the ability of the API to properly output logical matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-logical-output.R")

sFUN <- logical_sFUN
dFUN <- logical_dFUN
csFUN <- logical_csFUN
tsFUN <- logical_tsFUN

set.seed(12346)
test_that("Simple logical matrix output is okay", {
    check_write_all(sFUN, mode="logical")
    check_write_all(sFUN, nr=5, nc=30, mode="logical")
    check_write_all(sFUN, nr=30, nc=5, mode="logical")

    check_write_slice(sFUN, mode="logical")
    check_write_slice(sFUN, nr=5, nc=30, mode="logical")
    check_write_slice(sFUN, nr=30, nc=5, mode="logical")

    check_write_varslice(sFUN, mode="logical")
    check_write_varslice(sFUN, nr=5, nc=30, mode="logical")
    check_write_varslice(sFUN, nr=30, nc=5, mode="logical")

    check_write_indexed(sFUN, mode="logical")
    check_write_indexed(sFUN, nr=5, nc=30, mode="logical")
    check_write_indexed(sFUN, nr=30, nc=5, mode="logical")

    check_write_type(sFUN, mode="logical")
    check_write_errors(sFUN, mode="logical")

    check_write_all(sFUN, nr=0, nc=0, mode="logical")
    check_write_all(sFUN, nr=10, nc=0, mode="logical")
    check_write_all(sFUN, nr=0, nc=10, mode="logical")
})

#######################################################

set.seed(23457)
test_that("sparse logical matrix output is okay", {
    expect_s4_class(csFUN(), "lgCMatrix")

    check_write_all(csFUN, mode="logical")
    check_write_all(csFUN, nr=5, nc=30, mode="logical")
    check_write_all(csFUN, nr=30, nc=5, mode="logical")

    check_write_slice(csFUN, mode="logical")
    check_write_slice(csFUN, nr=5, nc=30, mode="logical")
    check_write_slice(csFUN, nr=30, nc=5, mode="logical")

    check_write_varslice(csFUN, mode="logical")
    check_write_varslice(csFUN, nr=5, nc=30, mode="logical")
    check_write_varslice(csFUN, nr=30, nc=5, mode="logical")

    check_write_indexed(csFUN, mode="logical")
    check_write_indexed(csFUN, nr=5, nc=30, mode="logical")
    check_write_indexed(csFUN, nr=30, nc=5, mode="logical")

    check_write_type(csFUN, mode="logical")
    check_write_errors(csFUN, mode="logical")

    check_write_all(csFUN, nr=0, nc=0, mode="logical")
    check_write_all(csFUN, nr=10, nc=0, mode="logical")
    check_write_all(csFUN, nr=0, nc=10, mode="logical")
})

#######################################################

test_that("Logical matrix mode choices are okay", {
    check_write_class(sFUN(), "matrix")
    check_write_class(dFUN(), "lgeMatrix")
    check_write_class(csFUN(), "lgCMatrix")
    check_write_class(tsFUN(), "lgTMatrix")
})
