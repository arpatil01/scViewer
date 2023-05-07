# This tests the ability of the API to properly access integer matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-integer-input.R")

sFUN <- integer_sFUN
rFUN <- integer_rFUN

#######################################################

set.seed(12345)
test_that("Simple integer matrix input is okay", {
    check_read_all(sFUN, mode="integer")
    check_read_all(sFUN, nr=5, nc=30, mode="integer")
    check_read_all(sFUN, nr=30, nc=5, mode="integer")

    check_read_slice(sFUN, mode="integer")
    check_read_slice(sFUN, nr=5, nc=30, mode="integer")
    check_read_slice(sFUN, nr=30, nc=5, mode="integer")

    check_read_varslice(sFUN, mode="integer")
    check_read_varslice(sFUN, nr=5, nc=30, mode="integer")
    check_read_varslice(sFUN, nr=30, nc=5, mode="integer")

    check_read_const(sFUN, mode="integer")
    check_read_const(sFUN, nr=5, nc=30, mode="integer")
    check_read_const(sFUN, nr=30, nc=5, mode="integer")

    check_read_multi(sFUN, mode="integer")
    check_read_multi(sFUN, nr=5, nc=30, mode="integer")
    check_read_multi(sFUN, nr=30, nc=5, mode="integer")

    check_read_type(sFUN, mode="integer")
    check_read_class(sFUN(), mode="integer", "matrix")

    check_read_errors(sFUN, mode="integer")
    check_read_all(sFUN, nr=0, nc=0, mode="integer")
    check_read_all(sFUN, nr=10, nc=0, mode="integer")
    check_read_all(sFUN, nr=0, nc=10, mode="integer")
})

#######################################################

set.seed(23456)
test_that("RLE integer matrix input (i.e., unknown) is okay", {
    expect_s4_class(rFUN(), "RleMatrix")

    check_read_all(rFUN, mode="integer")
    check_read_all(rFUN, nr=5, nc=30, mode="integer")
    check_read_all(rFUN, nr=30, nc=5, mode="integer")

    check_read_slice(rFUN, mode="integer")
    check_read_slice(rFUN, nr=5, nc=30, mode="integer")
    check_read_slice(rFUN, nr=30, nc=5, mode="integer")

    check_read_varslice(rFUN, mode="integer")
    check_read_varslice(rFUN, nr=5, nc=30, mode="integer")
    check_read_varslice(rFUN, nr=30, nc=5, mode="integer")

    check_read_const(rFUN, mode="integer")
    check_read_const(rFUN, nr=5, nc=30, mode="integer")
    check_read_const(rFUN, nr=30, nc=5, mode="integer")

    check_read_multi(rFUN, mode="integer")
    check_read_multi(rFUN, nr=5, nc=30, mode="integer")
    check_read_multi(rFUN, nr=30, nc=5, mode="integer")

    check_read_type(rFUN, mode="integer")
    check_read_class(rFUN(), mode="integer", "")

    check_read_errors(rFUN, mode="integer")
    check_read_all(rFUN, nr=0, nc=0, mode="integer")
    check_read_all(rFUN, nr=10, nc=0, mode="integer")
    check_read_all(rFUN, nr=0, nc=10, mode="integer")
})

test_that("RLE integer matrix input is okay with reduced block size", {
    old <- getAutoBlockSize()

    for (blocksize in c(100, 250, 500)) {
        setAutoBlockSize(blocksize)

        check_read_all(rFUN, mode="integer")
        check_read_all(rFUN, nr=5, nc=30, mode="integer")
        check_read_all(rFUN, nr=30, nc=5, mode="integer")

        check_read_slice(rFUN, mode="integer")
        check_read_slice(rFUN, nr=5, nc=30, mode="integer")
        check_read_slice(rFUN, nr=30, nc=5, mode="integer")

        check_read_varslice(rFUN, mode="integer")
        check_read_varslice(rFUN, nr=5, nc=30, mode="integer")
        check_read_varslice(rFUN, nr=30, nc=5, mode="integer")

        check_read_const(rFUN, mode="integer")
        check_read_const(rFUN, nr=5, nc=30, mode="integer")
        check_read_const(rFUN, nr=30, nc=5, mode="integer")

        check_read_multi(rFUN, mode="integer")
        check_read_multi(rFUN, nr=5, nc=30, mode="integer")
        check_read_multi(rFUN, nr=30, nc=5, mode="integer")

        check_read_type(rFUN, mode="integer")
        check_read_class(rFUN(), mode="integer", "")

        check_read_errors(rFUN, mode="integer")
        check_read_all(rFUN, nr=0, nc=0, mode="integer")
        check_read_all(rFUN, nr=10, nc=0, mode="integer")
        check_read_all(rFUN, nr=0, nc=10, mode="integer")
    }

    setAutoBlockSize(old)
})

#######################################################

set.seed(981347)
test_that("Delayed integer matrix input is okay", {
    delfuns <- c(
         delayed_funs(sFUN, DELAYED_FUN=function(x) { x + sample(nrow(x)) }), # known seed
         delayed_funs(rFUN, DELAYED_FUN=function(x) { x + sample(nrow(x)) })  # unknown seed
    )

    for (FUN in delfuns) {
        NR <- 10 + sample(10, 1)
        NC <- 10 + sample(10, 1)
        expect_s4_class(FUN(), "DelayedMatrix")

        check_read_all(FUN, NR, NC, mode="integer")

        check_read_slice(FUN, NR, NC, mode="integer")

        check_read_varslice(FUN, NR, NC, mode="integer")

        check_read_const(FUN, NR, NC, mode="integer")

        check_read_multi(FUN, NR, NC, mode="integer")

        check_read_type(FUN, NR, NC, mode="integer")
        check_read_class(FUN(), mode="integer", "DelayedMatrix")

        check_read_errors(FUN, NR, NC, mode="integer")
        check_read_all(FUN, nr=0, nc=0, mode="integer")
        check_read_all(FUN, nr=10, nc=0, mode="integer")
        check_read_all(FUN, nr=0, nc=10, mode="integer")
    }

    # Proper type check upon coercion!
    expect_identical("double", .Call("get_type", delfuns[[1]]()+1, PACKAGE="beachtest"))
})
