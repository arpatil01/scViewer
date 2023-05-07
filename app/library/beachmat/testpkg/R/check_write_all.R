#' @export
#' @importFrom testthat expect_identical
check_write_all <- function(FUN, ..., mode, out.class=NULL) {
    check_write_all_row(FUN(...), mode, out.class)
    check_write_all_col(FUN(...), mode, out.class)
    check_write_all_single(FUN(...), mode, out.class)
}

#' @importFrom testthat expect_identical
check_write_all_row <- function(test.mat, mode, out.class, FUN="set_row_all") {
    if (is.null(out.class)) {
        out.class <- as.character(class(test.mat))
    }
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    for (o in rranges) {
        out <- .Call(paste0(FUN, "_", mode), test.mat, o, PACKAGE="beachtest")

        REF <- ref
        REF[] <- get(mode)(1)
        REF[seq_along(o),] <- ref[o,]
        expect_matrix(REF, out[[1]], out.class)

        expect_identical(REF[o,,drop=FALSE], out[[2]])
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_write_all_col <- function(test.mat, mode, out.class, FUN="set_col_all") {
    if (is.null(out.class)) {
        out.class <- as.character(class(test.mat))
    }
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    cranges <- spawn_col_ordering(ncol(test.mat))
    for (o in cranges) {
        out <- .Call(paste0(FUN, "_", mode), test.mat, o, PACKAGE="beachtest")

        REF <- ref
        REF[] <- get(mode)(1)
        REF[,seq_along(o)] <- ref[,o]
        expect_matrix(REF, out[[1]], out.class)

        expect_identical(REF[,o,drop=FALSE], out[[2]])
    }
    return(invisible(NULL))
}

#' @importFrom testthat expect_identical
check_write_all_single <- function(test.mat, mode, out.class, FUN="set_single_all") {
    if (is.null(out.class)) {
        out.class <- as.character(class(test.mat))
    }
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL

    rranges <- spawn_row_ordering(nrow(test.mat))
    cranges <- spawn_col_ordering(ncol(test.mat))
    for (ro in rranges) {
        for (co in cranges) { 
            out <- .Call(paste0(FUN, "_", mode), test.mat, ro, co, PACKAGE="beachtest")

            REF <- ref
            REF[] <- get(mode)(1)
            REF[seq_along(ro),seq_along(co)] <- ref[ro, co]
            expect_matrix(REF, out[[1]], out.class)

            expect_identical(REF[ro,co,drop=FALSE], out[[2]])
        }
    }
    return(invisible(NULL))
}
