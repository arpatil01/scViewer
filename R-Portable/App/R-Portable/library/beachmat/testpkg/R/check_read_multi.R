#' @export
#' @importFrom testthat expect_identical
check_read_multi <- function(FUN, ..., mode) {
    check_read_all_rows(FUN(...), mode)
    check_read_all_cols(FUN(...), mode)
    check_read_slice_rows(FUN(...), mode)
    check_read_slice_cols(FUN(...), mode)
}

#' @importFrom testthat expect_identical
check_read_all_rows <- function(test.mat, mode, FUN="get_multirow_all") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    for (o in rranges) { 
        o <- sort(o)
        expect_identical(ref[o,,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o-1L, PACKAGE="beachtest"))
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_read_all_cols <- function(test.mat, mode, FUN="get_multicol_all") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- spawn_col_ordering(ncol(test.mat))
    for (o in cranges) {
        o <- sort(o)
        expect_identical(ref[,o,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o-1L, PACKAGE="beachtest"))
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_read_slice_rows <- function(test.mat, mode, FUN="get_multirow_slice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    cbounds <- spawn_col_bounds(ncol(test.mat))

    for (o in rranges) { 
        o <- sort(o)
        for (b in cbounds) {
            range <- b[1]:b[2]
            expect_identical(ref[o,range,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o-1L, b, PACKAGE="beachtest"))
        }
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_read_slice_cols <- function(test.mat, mode, FUN="get_multicol_slice") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- spawn_col_ordering(ncol(test.mat))
    rbounds <- spawn_row_bounds(nrow(test.mat))

    for (o in cranges) {
        o <- sort(o)
        for (b in rbounds) {
            range <- b[1]:b[2]
            expect_identical(ref[range,o,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o-1L, b, PACKAGE="beachtest"))
        }
    }
    return(invisible(NULL))
}




