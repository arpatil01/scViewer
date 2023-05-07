#' @export
#' @importFrom testthat expect_identical
check_read_all <- function(FUN, ..., mode) {
    check_read_all_row(FUN(...), mode)
    check_read_all_col(FUN(...), mode)
    check_read_all_single(FUN(...), mode)
}

#' @importFrom testthat expect_identical
check_read_all_row <- function(test.mat, mode, FUN="get_row_all") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    for (o in rranges) { 
        expect_identical(ref[o,,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o, PACKAGE="beachtest"))
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_read_all_col <- function(test.mat, mode, FUN="get_col_all") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- spawn_col_ordering(ncol(test.mat))
    for (o in cranges) {
        expect_identical(ref[,o,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, o, PACKAGE="beachtest"))
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_read_all_single <- function(test.mat, mode, FUN="get_single_all") {
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    cranges <- spawn_col_ordering(ncol(test.mat))
    for (ro in rranges) {
        for (co in cranges) { 
            expect_identical(ref[ro,co,drop=FALSE], .Call(paste0(FUN, "_", mode), test.mat, ro, co, PACKAGE="beachtest"))
        }
    }
    return(invisible(NULL))
}
