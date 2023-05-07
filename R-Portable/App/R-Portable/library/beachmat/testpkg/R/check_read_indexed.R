#' @export
#' @export
#' @importFrom testthat expect_identical
check_read_indexed <- function(FUN, ..., mode) {
    check_read_all_col(FUN(...), mode, FUN="get_indexed_all")
    check_read_slice_col(FUN(...), mode, FUN="get_indexed_slice")
    check_read_varslice_col(FUN(...), mode, FUN="get_indexed_varslice")
}
