# This tests the ability of the API to properly output integer matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-integer-output.R")

sFUN <- integer_sFUN
rFUN <- integer_rFUN

#######################################################

set.seed(12346)
test_that("Simple integer matrix output is okay", {
    check_write_all(sFUN, mode="integer")
    check_write_all(sFUN, nr=5, nc=30, mode="integer")
    check_write_all(sFUN, nr=30, nc=5, mode="integer")

    check_write_slice(sFUN, mode="integer")
    check_write_slice(sFUN, nr=5, nc=30, mode="integer")
    check_write_slice(sFUN, nr=30, nc=5, mode="integer")

    check_write_varslice(sFUN, mode="integer")
    check_write_varslice(sFUN, nr=5, nc=30, mode="integer")
    check_write_varslice(sFUN, nr=30, nc=5, mode="integer")

    check_write_indexed(sFUN, mode="integer")
    check_write_indexed(sFUN, nr=5, nc=30, mode="integer")
    check_write_indexed(sFUN, nr=30, nc=5, mode="integer")

    check_write_type(sFUN, mode="integer")
    check_write_errors(sFUN, mode="integer")

    check_write_all(sFUN, nr=0, nc=0, mode="integer")
    check_write_all(sFUN, nr=10, nc=0, mode="integer")
    check_write_all(sFUN, nr=0, nc=10, mode="integer")
})

#######################################################

test_that("Integer matrix mode choices are okay", {
    check_write_class(sFUN(), "matrix")
    check_write_class(rFUN(), "RleMatrix")
    check_write_class(rFUN()+1L, "DelayedMatrix")
})
