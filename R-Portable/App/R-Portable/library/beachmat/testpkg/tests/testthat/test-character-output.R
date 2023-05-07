# This tests the ability of the API to properly output character matrices of different types.
# library(testthat); library(beachtest); source("setup-fun.R"); source("test-character-output.R")

sFUN <- character_sFUN
rFUN <- character_rFUN

#######################################################

set.seed(12346)
test_that("Simple character matrix output is okay", {
    check_write_all(sFUN, mode="character")
    check_write_all(sFUN, nr=5, nc=30, mode="character")
    check_write_all(sFUN, nr=30, nc=5, mode="character")

    check_write_slice(sFUN, mode="character")
    check_write_slice(sFUN, nr=5, nc=30, mode="character")
    check_write_slice(sFUN, nr=30, nc=5, mode="character")

    check_write_varslice(sFUN, mode="character")
    check_write_varslice(sFUN, nr=5, nc=30, mode="character")
    check_write_varslice(sFUN, nr=30, nc=5, mode="character")

    check_write_indexed(sFUN, mode="character")
    check_write_indexed(sFUN, nr=5, nc=30, mode="character")
    check_write_indexed(sFUN, nr=30, nc=5, mode="character")

    check_write_type(sFUN, mode="character")
    check_write_errors(sFUN, mode="character")

    check_write_all(sFUN, nr=0, nc=0, mode="character")
    check_write_all(sFUN, nr=10, nc=0, mode="character")
    check_write_all(sFUN, nr=0, nc=10, mode="character")
})

#######################################################

test_that("Character matrix mode choices are okay", {
    check_write_class(sFUN(), "matrix")
    check_write_class(rFUN(), "RleMatrix")
    check_write_class(DelayedArray::tolower(rFUN()), "DelayedMatrix")
})
