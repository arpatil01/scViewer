#' @export
#' @importFrom testthat expect_true expect_error
check_read_errors <- function(FUN, ..., mode) {
    x <- FUN(...)
    cxxfun <- paste0("get_errors_", mode)

    expect_true(.Call(cxxfun, x, 0L, PACKAGE="beachtest"))
    expect_error(.Call(cxxfun, x, 1L, PACKAGE="beachtest"), "row index out of range")
    expect_error(.Call(cxxfun, x, -1L, PACKAGE="beachtest"), "column index out of range")
    expect_error(.Call(cxxfun, x, 2L, PACKAGE="beachtest"), "column start index is greater than column end index")
    expect_error(.Call(cxxfun, x, -2L, PACKAGE="beachtest"), "row start index is greater than row end index")
    expect_error(.Call(cxxfun, x, 3L, PACKAGE="beachtest"), "column end index out of range")
    expect_error(.Call(cxxfun, x, -3L, PACKAGE="beachtest"), "row end index out of range")

    cxxfun <- paste0("get_multi_errors_", mode)
    expect_error(.Call(cxxfun, x, 1L, PACKAGE="beachtest"), "row indices are not strictly increasing")
    expect_error(.Call(cxxfun, x, -1L, PACKAGE="beachtest"), "column indices are not strictly increasing")
    expect_error(.Call(cxxfun, x, 2L, PACKAGE="beachtest"), "row index out of range")
    expect_error(.Call(cxxfun, x, -2L, PACKAGE="beachtest"), "column index out of range")
    expect_error(.Call(cxxfun, x, 3L, PACKAGE="beachtest"), "column start index is greater than column end index")
    expect_error(.Call(cxxfun, x, -3L, PACKAGE="beachtest"), "row start index is greater than row end index")
    expect_error(.Call(cxxfun, x, 4L, PACKAGE="beachtest"), "column end index out of range")
    expect_error(.Call(cxxfun, x, -4L, PACKAGE="beachtest"), "row end index out of range")
    return(invisible(NULL))
}


