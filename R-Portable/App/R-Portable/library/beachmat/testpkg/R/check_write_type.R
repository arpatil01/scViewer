#' @export
check_write_type <- function(FUN, ..., mode, out.class=NULL) {
    convertible <- c("logical", "numeric", "integer")
    if (!mode %in% convertible) {
        return(invisible(NULL))
    }

    alternative <- setdiff(convertible, mode)
    for (alt.mode in alternative) {
        test.mat <- FUN(...)
        test.class <- out.class
        if (is.null(test.class)) {
            test.class <- as.character(class(test.mat))
        }

        ref <- as.matrix(test.mat)
        dimnames(ref) <- NULL

        # We use a passaging strategy as some data types do not have matrix representations (e.g., sparse matrices).
        if (alt.mode!="logical") {
            storage.mode(ref) <- alt.mode # Passaging it through the alternative mode. 
        } else {
            storage.mode(ref) <- "integer" # As 'logical' is 'int' in C++.
        }
        storage.mode(ref) <- mode # Converting it back to the current mode.
    
        ordering <- seq_len(nrow(test.mat))
        out <- .Call(paste0("set_row_", mode, "_via_", alt.mode), test.mat, ordering, PACKAGE="beachtest")
        expect_matrix(ref, out, test.class)

        ordering <- seq_len(ncol(test.mat))
        out <- .Call(paste0("set_col_", mode, "_via_", alt.mode), test.mat, ordering, PACKAGE="beachtest")
        expect_matrix(ref, out, test.class)
    }
    return(invisible(NULL))
}
