#' @export
#' @importFrom testthat expect_identical
check_read_varslice <- function(FUN, ..., mode) {
    check_read_varslice_row(FUN(...), mode)
    check_read_varslice_col(FUN(...), mode)
}

#' @importFrom testthat expect_identical
check_read_varslice_row <- function(test.mat, mode, FUN="get_row_varslice") {
    NCOL <- ncol(test.mat)
    rranges <- spawn_row_ordering(nrow(test.mat))

    for (o in rranges) {
        nentries <- length(o)
        bound1 <- sample(NCOL, nentries, replace=TRUE)
        bound2 <- sample(NCOL, nentries, replace=TRUE)
        cbounds <- cbind(pmin(bound1, bound2), pmax(bound1, bound2))

        ref <- get_reference_varslice(test.mat, o, cbounds, byrow=TRUE)
        expect_identical(ref, .Call(paste0(FUN, "_", mode), test.mat, o, cbounds, PACKAGE="beachtest"))
    }
}

#' @importFrom testthat expect_identical
check_read_varslice_col <- function(test.mat, mode, FUN="get_col_varslice") {
    NROW <- nrow(test.mat)
    cranges <- spawn_col_ordering(ncol(test.mat))

    for (o in cranges) {
        nentries <- length(o)
        bound1 <- sample(NROW, nentries, replace=TRUE)
        bound2 <- sample(NROW, nentries, replace=TRUE)
        rbounds <- cbind(pmin(bound1, bound2), pmax(bound1, bound2))

        ref <- get_reference_varslice(test.mat, o, rbounds, byrow=FALSE)
        expect_identical(ref, .Call(paste0(FUN, "_", mode), test.mat, o, rbounds, PACKAGE="beachtest"))
    }

    return(invisible(NULL))
}

get_reference_varslice <- function(test.mat, order, bounds, byrow=TRUE) {
    ref <- vector("list", length(order))
    for (i in seq_along(order)) {
        range <- bounds[i,1]:bounds[i,2]
        if (byrow) {
            out <- test.mat[order[i], range]
        } else {
            out <- test.mat[range, order[i]]
        }
        ref[[i]] <- as.vector(out)
    }
    return(ref)
}
