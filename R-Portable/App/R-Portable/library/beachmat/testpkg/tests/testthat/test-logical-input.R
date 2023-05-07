# This tests the ability of the API to properly access logical matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-logical-input.R")

sFUN <- logical_sFUN
dFUN <- logical_dFUN
csFUN <- logical_csFUN
tsFUN <- logical_tsFUN

#######################################################

set.seed(12345)
test_that("Simple logical matrix input is okay", {
    check_read_all(sFUN, mode="logical")
    check_read_all(sFUN, nr=5, nc=30, mode="logical")
    check_read_all(sFUN, nr=30, nc=5, mode="logical")

    check_read_slice(sFUN, mode="logical")
    check_read_slice(sFUN, nr=5, nc=30, mode="logical")
    check_read_slice(sFUN, nr=30, nc=5, mode="logical")

    check_read_varslice(sFUN, mode="logical")
    check_read_varslice(sFUN, nr=5, nc=30, mode="logical")
    check_read_varslice(sFUN, nr=30, nc=5, mode="logical")

    check_read_const(sFUN, mode="logical")
    check_read_const(sFUN, nr=5, nc=30, mode="logical")
    check_read_const(sFUN, nr=30, nc=5, mode="logical")

    check_read_indexed(sFUN, mode="logical")
    check_read_indexed(sFUN, nr=5, nc=30, mode="logical")
    check_read_indexed(sFUN, nr=30, nc=5, mode="logical")

    check_read_multi(sFUN, mode="logical")
    check_read_multi(sFUN, nr=5, nc=30, mode="logical")
    check_read_multi(sFUN, nr=30, nc=5, mode="logical")

    check_read_type(sFUN, mode="logical")
    check_read_class(sFUN(), mode="logical", "matrix")

    check_read_errors(sFUN, mode="logical")
    check_read_all(sFUN, nr=0, nc=0, mode="logical")
    check_read_all(sFUN, nr=10, nc=0, mode="logical")
    check_read_all(sFUN, nr=0, nc=10, mode="logical")
})

#######################################################

set.seed(13579)
test_that("Dense logical matrix input is okay", {
    expect_s4_class(dFUN(), "lgeMatrix")

    check_read_all(dFUN, mode="logical")
    check_read_all(dFUN, nr=5, nc=30, mode="logical")
    check_read_all(dFUN, nr=30, nc=5, mode="logical")

    check_read_slice(dFUN, mode="logical")
    check_read_slice(dFUN, nr=5, nc=30, mode="logical")
    check_read_slice(dFUN, nr=30, nc=5, mode="logical")

    check_read_varslice(dFUN, mode="logical")
    check_read_varslice(dFUN, nr=5, nc=30, mode="logical")
    check_read_varslice(dFUN, nr=30, nc=5, mode="logical")

    check_read_const(dFUN, mode="logical")
    check_read_const(dFUN, nr=5, nc=30, mode="logical")
    check_read_const(dFUN, nr=30, nc=5, mode="logical")

    check_read_indexed(dFUN, mode="logical")
    check_read_indexed(dFUN, nr=5, nc=30, mode="logical")
    check_read_indexed(dFUN, nr=30, nc=5, mode="logical")

    check_read_multi(dFUN, mode="logical")
    check_read_multi(dFUN, nr=5, nc=30, mode="logical")
    check_read_multi(dFUN, nr=30, nc=5, mode="logical")

    check_read_type(dFUN, mode="logical")
    check_read_class(dFUN(), mode="logical", "lgeMatrix")

    check_read_errors(dFUN, mode="logical")
    check_read_all(dFUN, nr=0, nc=0, mode="logical")
    check_read_all(dFUN, nr=10, nc=0, mode="logical")
    check_read_all(dFUN, nr=0, nc=10, mode="logical")
})

#######################################################

set.seed(23456)
test_that("Sparse logical matrix input is okay", {
    expect_s4_class(csFUN(), "lgCMatrix")

    check_read_all(csFUN, mode="logical")
    check_read_all(csFUN, nr=5, nc=30, mode="logical")
    check_read_all(csFUN, nr=30, nc=5, mode="logical")

    check_read_slice(csFUN, mode="logical")
    check_read_slice(csFUN, nr=5, nc=30, mode="logical")
    check_read_slice(csFUN, nr=30, nc=5, mode="logical")

    check_read_varslice(csFUN, mode="logical")
    check_read_varslice(csFUN, nr=5, nc=30, mode="logical")
    check_read_varslice(csFUN, nr=30, nc=5, mode="logical")

    check_read_const(csFUN, mode="logical")
    check_read_const(csFUN, nr=5, nc=30, mode="logical")
    check_read_const(csFUN, nr=30, nc=5, mode="logical")

    check_read_indexed(csFUN, mode="logical")
    check_read_indexed(csFUN, nr=5, nc=30, mode="logical")
    check_read_indexed(csFUN, nr=30, nc=5, mode="logical")

    check_read_multi(csFUN, mode="logical")
    check_read_multi(csFUN, nr=5, nc=30, mode="logical")
    check_read_multi(csFUN, nr=30, nc=5, mode="logical")

    check_read_type(csFUN, mode="logical")
    check_read_class(csFUN(), mode="logical", "lgCMatrix")

    check_read_errors(csFUN, mode="logical")
    check_read_all(csFUN, nr=0, nc=0, mode="logical")
    check_read_all(csFUN, nr=10, nc=0, mode="logical")
    check_read_all(csFUN, nr=0, nc=10, mode="logical")
})

#######################################################

set.seed(23456)
test_that("lgTMatrix (i.e., unknown) input is okay", {
    expect_s4_class(tsFUN(), "lgTMatrix")

    check_read_all(tsFUN, mode="logical")
    check_read_all(tsFUN, nr=5, nc=30, mode="logical")
    check_read_all(tsFUN, nr=30, nc=5, mode="logical")

    check_read_slice(tsFUN, mode="logical")
    check_read_slice(tsFUN, nr=5, nc=30, mode="logical")
    check_read_slice(tsFUN, nr=30, nc=5, mode="logical")

    check_read_varslice(tsFUN, mode="logical")
    check_read_varslice(tsFUN, nr=5, nc=30, mode="logical")
    check_read_varslice(tsFUN, nr=30, nc=5, mode="logical")

    check_read_const(tsFUN, mode="logical")
    check_read_const(tsFUN, nr=5, nc=30, mode="logical")
    check_read_const(tsFUN, nr=30, nc=5, mode="logical")

    check_read_indexed(tsFUN, mode="logical")
    check_read_indexed(tsFUN, nr=5, nc=30, mode="logical")
    check_read_indexed(tsFUN, nr=30, nc=5, mode="logical")

    check_read_multi(tsFUN, mode="logical")
    check_read_multi(tsFUN, nr=5, nc=30, mode="logical")
    check_read_multi(tsFUN, nr=30, nc=5, mode="logical")

    check_read_type(tsFUN, mode="logical")
    check_read_class(tsFUN(), mode="logical", "")

    check_read_errors(tsFUN, mode="logical")
    check_read_all(tsFUN, nr=0, nc=0, mode="logical")
    check_read_all(tsFUN, nr=10, nc=0, mode="logical")
    check_read_all(tsFUN, nr=0, nc=10, mode="logical")
})

test_that("lgTMatrix input is okay with reduced block size", {
    old <- getAutoBlockSize()

    for (blocksize in c(100, 250, 500)) {
        setAutoBlockSize(blocksize)

        check_read_all(tsFUN, mode="logical")
        check_read_all(tsFUN, nr=5, nc=30, mode="logical")
        check_read_all(tsFUN, nr=30, nc=5, mode="logical")

        check_read_slice(tsFUN, mode="logical")
        check_read_slice(tsFUN, nr=5, nc=30, mode="logical")
        check_read_slice(tsFUN, nr=30, nc=5, mode="logical")

        check_read_varslice(tsFUN, mode="logical")
        check_read_varslice(tsFUN, nr=5, nc=30, mode="logical")
        check_read_varslice(tsFUN, nr=30, nc=5, mode="logical")

        check_read_const(tsFUN, mode="logical")
        check_read_const(tsFUN, nr=5, nc=30, mode="logical")
        check_read_const(tsFUN, nr=30, nc=5, mode="logical")

        check_read_indexed(tsFUN, mode="logical")
        check_read_indexed(tsFUN, nr=5, nc=30, mode="logical")
        check_read_indexed(tsFUN, nr=30, nc=5, mode="logical")

        check_read_multi(tsFUN, mode="logical")
        check_read_multi(tsFUN, nr=5, nc=30, mode="logical")
        check_read_multi(tsFUN, nr=30, nc=5, mode="logical")

        check_read_type(tsFUN, mode="logical")
        check_read_class(tsFUN(), mode="logical", "")

        check_read_errors(tsFUN, mode="logical")
        check_read_all(tsFUN, nr=0, nc=0, mode="logical")
        check_read_all(tsFUN, nr=10, nc=0, mode="logical")
        check_read_all(tsFUN, nr=0, nc=10, mode="logical")
    }

    setAutoBlockSize(old)
})

#######################################################

set.seed(981347)
test_that("Delayed logical matrix input is okay", {
    delfuns <- c(
        delayed_funs(sFUN, DELAYED_FUN=function(m) !m), # known seed
        delayed_funs(tsFUN, DELAYED_FUN=function(m) !m) # unknown seed
    )

    for (FUN in delfuns) {
        NR <- 10 + sample(10, 1)
        NC <- 10 + sample(10, 1)
        expect_s4_class(FUN(NR, NC), "DelayedMatrix")

        check_read_all(FUN, NR, NC, mode="logical")

        check_read_slice(FUN, NR, NC, mode="logical")

        check_read_varslice(FUN, NR, NC, mode="logical")

        check_read_const(FUN, NR, NC, mode="logical")

        check_read_indexed(FUN, NR, NC, mode="logical")

        check_read_multi(FUN, NR, NC, mode="logical")

        check_read_type(FUN, NR, NC, mode="logical")
        check_read_class(FUN(), mode="logical", "DelayedMatrix")

        check_read_errors(FUN, NR, NC, mode="logical")
        check_read_all(FUN, nr=0, nc=0, mode="logical")
        check_read_all(FUN, nr=10, nc=0, mode="logical")
        check_read_all(FUN, nr=0, nc=10, mode="logical")
    }

    # Proper type check upon coercion!
    expect_identical("integer", .Call("get_type", delfuns[[1]]() + 1L, PACKAGE="beachtest"))
})
