#' @export
#' @importFrom testthat expect_identical
check_read_type <- function(FUN, ..., mode) {
    expect_identical(.Call("get_type", FUN(...), PACKAGE="beachtest"), if (mode=="numeric") "double" else mode)

    convertible <- c("logical", "numeric", "integer")
    if (mode %in% convertible) {
        alternative <- setdiff(convertible, mode)

        for (alt.mode in alternative) {
            if (alt.mode=="logical") {
                next # C++ numeric/integer->logical goes via int, which is not quite right.
            }

            test.mat <- FUN(...)
            ref <- as.matrix(test.mat)
            dimnames(ref) <- NULL
            storage.mode(ref) <- alt.mode
        
            ordering <- seq_len(nrow(test.mat))
            out <- .Call(paste0("get_row_", mode, "_to_", alt.mode), test.mat, ordering, PACKAGE="beachtest")
            expect_identical(ref, out)

            ordering <- seq_len(ncol(test.mat))
            out <- .Call(paste0("get_col_", mode, "_to_", alt.mode), test.mat, ordering, PACKAGE="beachtest")
            expect_identical(ref, out)
        }
    }
    return(invisible(NULL))
}
