# This tests the ability of the API to properly access character matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-character-input.R")

sFUN <- character_sFUN
rFUN <- character_rFUN

#######################################################

set.seed(12345)
test_that("Simple character matrix input is okay", {
    check_read_all(sFUN, mode="character")
    check_read_all(sFUN, nr=5, nc=30, mode="character")
    check_read_all(sFUN, nr=30, nc=5, mode="character")

    check_read_slice(sFUN, mode="character")
    check_read_slice(sFUN, nr=5, nc=30, mode="character")
    check_read_slice(sFUN, nr=30, nc=5, mode="character")

    check_read_varslice(sFUN, mode="character")
    check_read_varslice(sFUN, nr=5, nc=30, mode="character")
    check_read_varslice(sFUN, nr=30, nc=5, mode="character")

    check_read_const(sFUN, mode="character")
    check_read_const(sFUN, nr=5, nc=30, mode="character")
    check_read_const(sFUN, nr=30, nc=5, mode="character")

    check_read_multi(sFUN, mode="character")
    check_read_multi(sFUN, nr=5, nc=30, mode="character")
    check_read_multi(sFUN, nr=30, nc=5, mode="character")

    check_read_type(sFUN, mode="character")
    check_read_class(sFUN(), mode="character", "matrix")

    check_read_errors(sFUN, mode="character")
    check_read_all(sFUN, nr=0, nc=0, mode="character")
    check_read_all(sFUN, nr=10, nc=0, mode="character")
    check_read_all(sFUN, nr=0, nc=10, mode="character")
})

#######################################################
# Testing RLE matrices, treated as unknown.

set.seed(23456)
test_that("RLE character matrix input is okay", {
    expect_s4_class(rFUN(), "RleMatrix")

    check_read_all(rFUN, mode="character")
    check_read_all(rFUN, nr=5, nc=30, mode="character")
    check_read_all(rFUN, nr=30, nc=5, mode="character")

    check_read_slice(rFUN, mode="character")
    check_read_slice(rFUN, nr=5, nc=30, mode="character")
    check_read_slice(rFUN, nr=30, nc=5, mode="character")

    check_read_varslice(rFUN, mode="character")
    check_read_varslice(rFUN, nr=5, nc=30, mode="character")
    check_read_varslice(rFUN, nr=30, nc=5, mode="character")

    check_read_const(rFUN, mode="character")
    check_read_const(rFUN, nr=5, nc=30, mode="character")
    check_read_const(rFUN, nr=30, nc=5, mode="character")

    check_read_multi(rFUN, mode="character")
    check_read_multi(rFUN, nr=5, nc=30, mode="character")
    check_read_multi(rFUN, nr=30, nc=5, mode="character")

    check_read_type(rFUN, mode="character")
    check_read_class(rFUN(), mode="character", "")

    check_read_errors(rFUN, mode="character")
    check_read_all(rFUN, nr=0, nc=0, mode="character")
    check_read_all(rFUN, nr=10, nc=0, mode="character")
    check_read_all(rFUN, nr=0, nc=10, mode="character")
})

test_that("RLE character matrix input is okay with reduced block size", {
    old <- getAutoBlockSize()

    for (blocksize in c(100, 250, 500)) {
        setAutoBlockSize(blocksize)

        check_read_all(rFUN, mode="character")
        check_read_all(rFUN, nr=5, nc=30, mode="character")
        check_read_all(rFUN, nr=30, nc=5, mode="character")

        check_read_slice(rFUN, mode="character")
        check_read_slice(rFUN, nr=5, nc=30, mode="character")
        check_read_slice(rFUN, nr=30, nc=5, mode="character")

        check_read_varslice(rFUN, mode="character")
        check_read_varslice(rFUN, nr=5, nc=30, mode="character")
        check_read_varslice(rFUN, nr=30, nc=5, mode="character")

        check_read_const(rFUN, mode="character")
        check_read_const(rFUN, nr=5, nc=30, mode="character")
        check_read_const(rFUN, nr=30, nc=5, mode="character")

        check_read_multi(rFUN, mode="character")
        check_read_multi(rFUN, nr=5, nc=30, mode="character")
        check_read_multi(rFUN, nr=30, nc=5, mode="character")

        check_read_type(rFUN, mode="character")
        check_read_class(rFUN(), mode="character", "")

        check_read_errors(rFUN, mode="character")
        check_read_all(rFUN, nr=0, nc=0, mode="character")
        check_read_all(rFUN, nr=10, nc=0, mode="character")
        check_read_all(rFUN, nr=0, nc=10, mode="character")
    }

    setAutoBlockSize(old)
})

#######################################################

set.seed(981347)
test_that("Delayed character matrix input is okay", {
    delfuns <- c(
        delayed_funs(sFUN, DELAYED_FUN=DelayedArray::tolower), # known seed
        delayed_funs(rFUN, DELAYED_FUN=DelayedArray::tolower)  # unknown seed
    )

    for (FUN in delfuns) {
        NR <- 10 + sample(10, 1)
        NC <- 10 + sample(10, 1)
        expect_s4_class(FUN(), "DelayedMatrix")

        check_read_all(FUN, NR, NC, mode="character")

        check_read_slice(FUN, NR, NC, mode="character")

        check_read_varslice(FUN, NR, NC, mode="character")
   
        check_read_const(FUN, NR, NC, mode="character")
    
        check_read_multi(FUN, NR, NC, mode="character")

        check_read_type(FUN, NR, NC, mode="character")
        check_read_class(FUN(), mode="character", "DelayedMatrix")

        check_read_errors(FUN, NR, NC, mode="character")
        check_read_all(FUN, nr=0, nc=0, mode="character")
        check_read_all(FUN, nr=10, nc=0, mode="character")
        check_read_all(FUN, nr=0, nc=10, mode="character")
    }

    # Proper type check upon coercion!
    expect_identical("logical", .Call("get_type", delfuns[[1]]()=="A", PACKAGE="beachtest"))
})
