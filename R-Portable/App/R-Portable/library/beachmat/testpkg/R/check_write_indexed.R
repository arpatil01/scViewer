#' @export
check_write_indexed <- function(FUN, ..., mode, out.class=NULL) {
    check_write_indexed_row(FUN(...), mode, out.class)
    check_write_indexed_col(FUN(...), mode, out.class)
}

check_write_indexed_row <- function(test.mat, mode, out.class, FUN="set_row_indexed") {
    if (is.null(out.class)) {
        out.class <- as.character(class(test.mat))
    }
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    rranges <- spawn_row_ordering(nrow(test.mat))

    for (o in rranges) {
        REF <- ref
        REF[] <- get(mode)(1)

        subr <- vector("list", length(o))
        for (i in seq_along(o)) {
            n <- sample(ncol(test.mat), 1)
            curI <- sample(ncol(test.mat), n, replace=TRUE)
            curX <- sample(ref[o[i],], n, replace=TRUE)
            subr[[i]] <- list(curI, curX)

            to.use <- !duplicated(curI, fromLast=TRUE) # Last elements overwrite earlier elements.
            REF[o[i], curI[to.use]] <- curX[to.use]
        }
        
        out <- .Call(paste0(FUN, "_", mode), test.mat, o, subr)
        expect_matrix(REF, out, out.class)
    }

    return(invisible(NULL))
}

check_write_indexed_col <- function(test.mat, mode, out.class, FUN="set_col_indexed") {
    if (is.null(out.class)) {
        out.class <- as.character(class(test.mat))
    }
    ref <- as.matrix(test.mat)
    dimnames(ref) <- NULL
    cranges <- spawn_col_ordering(ncol(test.mat))

    for (o in cranges) {
        REF <- ref
        REF[] <- get(mode)(1)

        subr <- vector("list", length(o))
        for (i in seq_along(o)) {
            n <- sample(nrow(test.mat), 1)
            curI <- sample(nrow(test.mat), n, replace=TRUE)
            curX <- sample(ref[,o[i]], n, replace=TRUE)
            subr[[i]] <- list(curI, curX)

            to.use <- !duplicated(curI, fromLast=TRUE) # Last elements overwrite earlier elements.
            REF[curI[to.use], o[i]] <- curX[to.use]
        }
        
        out <- .Call(paste0(FUN, "_", mode), test.mat, o, subr)
        expect_matrix(REF, out, out.class)
    }

    return(invisible(NULL))
}
