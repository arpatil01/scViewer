#' @export
#' @import methods
setClass("AaronMatrix", representation(data="matrix"))

#' @export
#' @importFrom methods new
AaronMatrix <- function(x) {
    new("AaronMatrix", data=x)
}

#' @importFrom methods show
setMethod("show", "AaronMatrix", function(object) {
    cat(sprintf("Aaron's matrix [%i x %i, %s]\n", nrow(object), ncol(object), type(object)))
})

#' @export
setMethod("as.matrix", "AaronMatrix", function(x) x@data)

#' @export
#' @method as.matrix AaronMatrix
as.matrix.AaronMatrix <- function(x) {
    x@data
}

#' @export
setMethod("dim", "AaronMatrix", function(x) dim(x@data))

#' @export
#' @importFrom methods new
setMethod("[", "AaronMatrix", function(x, i, j, ..., drop = TRUE) {
    y <- x@data
    if (!missing(i)) {
        y <- y[i,,drop=FALSE]
    } 
    if (!missing(j)) {
        y <- y[,j,drop=FALSE]
    }
    if (drop && any(dim(y)==1L)) {
        return(drop(y))
    }
    new("AaronMatrix", data=y)
})

#' @export
#' @importFrom DelayedArray type
setMethod("type", "AaronMatrix", function(x) typeof(x@data))

#' @importFrom Rcpp sourceCpp
#' @useDynLib morebeach
NULL
