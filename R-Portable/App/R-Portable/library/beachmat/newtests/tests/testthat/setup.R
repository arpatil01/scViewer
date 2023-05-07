library(DelayedArray)
library(Matrix)

SPAWN <- function(nr, nc, mode) {
    mat <- Matrix::rsparsematrix(nr, nc, density=0.3)

    if (mode!=1L) {
        if (mode==0L) {
            mat <- mat != 0
        }

        output <- list(
            as.matrix(mat),
            mat,
            as(mat, "SparseArraySeed")
        )
    } else {
        mat <- as.matrix(round(mat)) # necessary to ignore zeroes in SparseArraySeed coercion.

        output <- list(
            mat,
            as(mat, "SparseArraySeed")
        )

        storage.mode(output[[1]]) <- "integer"
        storage.mode(output[[2]]@nzdata) <- "integer"
    }

    # Testing a scrambled version of the SparseArraySeed.
    n <- length(output)
    sas <- output[[n]]
    shuffle <- sample(length(sas@nzdata))
    sas@nzdata <- sas@nzdata[shuffle]
    sas@nzindex <- sas@nzindex[shuffle,,drop=FALSE]
    output[[n + 1L]] <- sas

    output

}

CONVERT <- function(x, mode) {
    if (mode==0) {
        storage.mode(x) <- "integer" # as logical conversion goes via integer truncation.
        x <- x != 0L
    } else {
        storage.mode(x) <- c("integer", "double")[mode]
    }
    dimnames(x) <- NULL
    x
}

CHECK_IDENTITY <- function(ref, mat, mode) {
    ref <- CONVERT(ref, mode)
    dimnames(ref) <- NULL

    if (mode==0L) {
        mat <- !!mat # due to the fact that logicals are integers, so non-1 values behave weirdly.
    }
    expect_identical(ref, mat)
}

CHECK_SPARSE_IDENTITY <- function(ref, mat, mode) {
    ref <- CONVERT(ref, mode) 
    if (mode==0L) { 
        ref <- as(ref, "lgCMatrix")

        expect_s4_class(mat, "lgCMatrix")
        mat <- as(!!as.matrix(mat), "lgCMatrix") # for much the same reasons as above.
    } else {
        ref <- as(ref, "dgCMatrix")
    }
    expect_identical(ref, mat)
}

SLICE_COLUMNS <- function(x, order, starts, ends) {
    for (o in order) {
        y <- x[starts[o]:ends[o],o]
        x[,o] <- vector(typeof(y), 1L)
        x[starts[o]:ends[o],o] <- y
    }
    x
}

SLICE_ROWS <- function(x, order, starts, ends) {
    for (o in order) {
        y <- x[o,starts[o]:ends[o]]
        x[o,] <- vector(typeof(y), 1L)
        x[o,starts[o]:ends[o]] <- y 
    }
    x
}
