#' @export
#' @importFrom testthat expect_identical
check_read_class <- function(test.mat, mode, expected) {
    expect_identical(expected, .Call(paste0("get_class_", mode), test.mat, PACKAGE="beachtest"))
}
