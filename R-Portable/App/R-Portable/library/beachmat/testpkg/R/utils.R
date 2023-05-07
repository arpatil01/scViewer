spawn_row_ordering <- function(NROW) {
    rranges <- list(forward=seq_len(NROW), random=sample(NROW, NROW), subset=sample(NROW, NROW/2L))
    rranges$reverse <- rev(rranges$forward)
    rranges
}

spawn_col_ordering <- function(NCOL) {
    cranges <- list(forward=seq_len(NCOL), random=sample(NCOL, NCOL), subset=sample(NCOL, NCOL/2L))
    cranges$reverse <- rev(cranges$forward)
    cranges
}

spawn_row_bounds <- function(NROW) {
    list(full=c(1L, NROW), left=c(1L, floor(NROW/2L)), right=c(ceiling(NROW/2L), NROW), middle=sort(sample(NROW, 2)), single=rep(sample(NROW, 1), 2))
}

spawn_col_bounds <- function(NCOL) {
    list(full=c(1L, NCOL), left=c(1L, floor(NCOL/2L)), right=c(ceiling(NCOL/2L), NCOL), middle=sort(sample(NCOL, 2)), single=rep(sample(NCOL, 1), 2))
}

#' @importFrom testthat expect_s4_class expect_identical
expect_matrix <- function(truth, observed, xclass) {
    expect_identical(xclass, as.character(class(observed)))

    ref <- as.matrix(observed)
    dimnames(ref) <- NULL
    expect_identical(ref, truth)
}
