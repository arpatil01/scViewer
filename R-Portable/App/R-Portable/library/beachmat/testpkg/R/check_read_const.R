#' @export
#' @export
#' @importFrom testthat expect_identical
check_read_const <- function(FUN, ..., mode) {
    check_read_all_col(FUN(...), mode, FUN="get_const_all")
    check_read_slice_col(FUN(...), mode, FUN="get_const_slice")
    check_read_varslice_col(FUN(...), mode, FUN="get_const_varslice")
}
