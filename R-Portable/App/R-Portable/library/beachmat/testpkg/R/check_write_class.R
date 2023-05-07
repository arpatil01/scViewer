#' @export
#' @importFrom testthat expect_identical
check_write_class <- function(test.mat, expected) {
    expect_identical(expected, .Call("set_class_by_sexp", test.mat, PACKAGE="beachtest"))
}
