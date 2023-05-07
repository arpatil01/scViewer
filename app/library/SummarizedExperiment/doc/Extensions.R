## ---- echo=FALSE, results="hide"----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE)

## ---- echo=FALSE--------------------------------------------------------------
library(SummarizedExperiment)
library(testthat)

## -----------------------------------------------------------------------------
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.CountSE <- setClass("CountSE", contains="SummarizedExperiment")

## -----------------------------------------------------------------------------
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
CountSE <- function(counts, ...) {
    se <- SummarizedExperiment(list(counts=counts), ...)
    .CountSE(se)
}

## -----------------------------------------------------------------------------
setValidity2("CountSE", function(object) {
    msg <- NULL

    if (assayNames(object)[1] != "counts") {
        msg <- c(msg, "'counts' must be first assay")
    }

    if (min(assay(object)) < 0) {
        msg <- c(msg, "negative values in 'counts'")
    }

    if (is.null(msg)) {
        TRUE
    } else msg
})

## -----------------------------------------------------------------------------
CountSE(matrix(rpois(100, lambda=1), ncol=5))

## ---- error=TRUE--------------------------------------------------------------
CountSE(matrix(rnorm(100), ncol=5))

## -----------------------------------------------------------------------------
#' @export
setGeneric("negcounts", function(x, ...) standardGeneric("negcounts"))

## -----------------------------------------------------------------------------
#' @export
#' @importFrom SummarizedExperiment assay
setMethod("negcounts", "CountSE", function(x, withDimnames=TRUE) {
    -assay(x, withDimnames=withDimnames)
})

## -----------------------------------------------------------------------------
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ExampleClass <- setClass("ExampleClass",
    slots= representation(
        rowVec="integer",
        colVec="integer",
        rowToRowMat="matrix",
        colToColMat="matrix",
        rowToColMat="matrix",
        colToRowMat="matrix"
    ),
    contains="SummarizedExperiment"
)

## -----------------------------------------------------------------------------
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
ExampleClass <- function(
    rowVec=integer(0), 
    colVec=integer(0),
    rowToRowMat=matrix(0,0,0),
    colToColMat=matrix(0,0,0),
    rowToColMat=matrix(0,0,0),
    colToRowMat=matrix(0,0,0),
    ...)
{
    se <- SummarizedExperiment(...)
    .ExampleClass(se, rowVec=rowVec, colVec=colVec,
        rowToRowMat=rowToRowMat, colToColMat=colToColMat, 
        rowToColMat=rowToColMat, colToRowMat=colToRowMat)
}

## -----------------------------------------------------------------------------
#' @export
setGeneric("rowVec", function(x, ...) standardGeneric("rowVec"))

#' @export
setGeneric("colVec", function(x, ...) standardGeneric("colVec"))

## -----------------------------------------------------------------------------
#' @export
setMethod("rowVec", "ExampleClass", function(x, withDimnames=TRUE) {
    out <- x@rowVec
    if (withDimnames) 
        names(out) <- rownames(x)
    out
})

#' @export
setMethod("colVec", "ExampleClass", function(x, withDimnames=TRUE) {
    out <- x@colVec
    if (withDimnames) 
        names(out) <- colnames(x)
    out
})

## -----------------------------------------------------------------------------
#' @export
setGeneric("rowToRowMat", function(x, ...) standardGeneric("rowToRowMat"))

#' @export
setGeneric("colToColMat", function(x, ...) standardGeneric("colToColMat"))

#' @export
setGeneric("rowToColMat", function(x, ...) standardGeneric("rowToColMat"))

#' @export
setGeneric("colToRowMat", function(x, ...) standardGeneric("colToRowMat"))

## -----------------------------------------------------------------------------
#' @export
setMethod("rowToRowMat", "ExampleClass", function(x, withDimnames=TRUE) {
    out <- x@rowToRowMat
    if (withDimnames) 
        rownames(out) <- rownames(x)
    out
})

#' @export
setMethod("colToColMat", "ExampleClass", function(x, withDimnames=TRUE) {
    out <- x@colToColMat
    if (withDimnames) 
        colnames(out) <- colnames(x)
    out
})

#' @export
setMethod("rowToColMat", "ExampleClass", function(x, withDimnames=TRUE) {
    out <- x@rowToColMat
    if (withDimnames) 
        rownames(out) <- colnames(x)
    out
})

#' @export
setMethod("colToRowMat", "ExampleClass", function(x, withDimnames=TRUE) {
    out <- x@colToRowMat
    if (withDimnames) 
        colnames(out) <- rownames(x)
    out
})

## -----------------------------------------------------------------------------
#' @export
#' @importMethodsFrom SummarizedExperiment rowData
setMethod("rowData", "ExampleClass", function(x, ...) {
    out <- callNextMethod()
    
    # Do something extra here.
    out$extra <- runif(nrow(out))

    # Returning the rowData object.
    out
})

## -----------------------------------------------------------------------------
#' @importFrom BiocGenerics NCOL NROW
setValidity2("ExampleClass", function(object) {
    NR <- NROW(object)
    NC <- NCOL(object)
    msg <- NULL

    # 1D
    if (length(rowVec(object, withDimnames=FALSE)) != NR) {
        msg <- c(msg, "'rowVec' should have length equal to the number of rows")
    }
    if (length(colVec(object, withDimnames=FALSE)) != NC) {
        msg <- c(
            msg, "'colVec' should have length equal to the number of columns"
        )
    }

    # 2D
    if (NROW(rowToRowMat(object, withDimnames=FALSE)) != NR) {
        msg <- c(
            msg, "'nrow(rowToRowMat)' should be equal to the number of rows"
        )
    }
    if (NCOL(colToColMat(object, withDimnames=FALSE)) != NC) {
        msg <- c(
            msg, "'ncol(colToColMat)' should be equal to the number of columns"
        )
    }
    if (NROW(rowToColMat(object, withDimnames=FALSE)) != NC) {
        msg <- c(
            msg, "'nrow(rowToColMat)' should be equal to the number of columns"
        )
    }
    if (NCOL(colToRowMat(object, withDimnames=FALSE)) != NR) {
        msg <- c(
            msg, "'ncol(colToRowMat)' should be equal to the number of rows"
        )
    }

    if (length(msg)) {
        msg
    } else TRUE
})

## -----------------------------------------------------------------------------
#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "ExampleClass", function(object) {
    callNextMethod()
    cat(
        "rowToRowMat has ", ncol(rowToRowMat(object)), " columns\n",
        "colToColMat has ", nrow(colToColMat(object)), " rows\n",
        "rowToColMat has ", ncol(rowToRowMat(object)), " columns\n",
        "colToRowMat has ", ncol(rowToRowMat(object)), " rows\n",
        sep=""
    )
})

## -----------------------------------------------------------------------------
#' @export
setGeneric("rowVec<-", function(x, ..., value) standardGeneric("rowVec<-"))

#' @export
setGeneric("colVec<-", function(x, ..., value) standardGeneric("colVec<-"))

## -----------------------------------------------------------------------------
#' @export
setReplaceMethod("rowVec", "ExampleClass", function(x, value) {
    x@rowVec <- value
    validObject(x)
    x
})

#' @export
setReplaceMethod("colVec", "ExampleClass", function(x, value) {
    x@colVec <- value
    validObject(x)
    x
})

## -----------------------------------------------------------------------------
#' @export
setGeneric("rowToRowMat<-", function(x, ..., value)
    standardGeneric("rowToRowMat<-")
)

#' @export
setGeneric("colToColMat<-", function(x, ..., value)
    standardGeneric("colToColMat<-")
)

#' @export
setGeneric("rowToColMat<-", function(x, ..., value) 
    standardGeneric("rowToColMat<-")
)

#' @export
setGeneric("colToRowMat<-", function(x, ..., value)
    standardGeneric("colToRowMat<-")
)

## -----------------------------------------------------------------------------
#' @export
setReplaceMethod("rowToRowMat", "ExampleClass", function(x, value) {
    x@rowToRowMat <- value
    validObject(x)
    x
})

#' @export
setReplaceMethod("colToColMat", "ExampleClass", function(x, value) {
    x@colToColMat <- value
    validObject(x)
    x
})

#' @export
setReplaceMethod("rowToColMat", "ExampleClass", function(x, value) {
    x@rowToColMat <- value
    validObject(x)
    x
})

#' @export
setReplaceMethod("colToRowMat", "ExampleClass", function(x, value) {
    x@colToRowMat <- value
    validObject(x)
    x
})

## -----------------------------------------------------------------------------
#' @export
#' @importMethodsFrom SummarizedExperiment "rowData<-"
setReplaceMethod("rowData", "ExampleClass", function(x, ..., value) {
    y <- callNextMethod() # returns a modified ExampleClass
    
    # Do something extra here.
    message("hi!\n")

    y
})

## -----------------------------------------------------------------------------
#' @export
#' @importFrom BiocGenerics normalize
setMethod("normalize", "ExampleClass", function(object) {
    # do something exciting, i.e., flip the signs
    new.row <- -rowVec(object, withDimnames=FALSE) 
    new.col <- -colVec(object, withDimnames=FALSE)
    BiocGenerics:::replaceSlots(object, rowVec=new.row, 
        colVec=new.col, check=FALSE)
})

## -----------------------------------------------------------------------------
#' @export
setMethod("[", "ExampleClass", function(x, i, j, drop=TRUE) {
    rv <- rowVec(x, withDimnames=FALSE)
    cv <- colVec(x, withDimnames=FALSE)
    rrm <- rowToRowMat(x, withDimnames=FALSE)
    ccm <- colToColMat(x, withDimnames=FALSE)
    rcm <- rowToColMat(x, withDimnames=FALSE)
    crm <- colToRowMat(x, withDimnames=FALSE)

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(x), fmt
            )
        }
        i <- as.vector(i)
        rv <- rv[i]
        rrm <- rrm[i,,drop=FALSE]
        crm <- crm[,i,drop=FALSE]
    }

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                j, colnames(x), fmt
            )
        }
        j <- as.vector(j)
        cv <- cv[j]
        ccm <- ccm[,j,drop=FALSE]
        rcm <- rcm[j,,drop=FALSE]
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, rowVec=rv, colVec=cv,
        rowToRowMat=rrm, colToColMat=ccm, 
        rowToColMat=rcm, colToRowMat=crm, check=FALSE)
})

## -----------------------------------------------------------------------------
#' @export
setReplaceMethod("[", c("ExampleClass", "ANY", "ANY", "ExampleClass"),
        function(x, i, j, ..., value) {
    rv <- rowVec(x, withDimnames=FALSE)
    cv <- colVec(x, withDimnames=FALSE)
    rrm <- rowToRowMat(x, withDimnames=FALSE)
    ccm <- colToColMat(x, withDimnames=FALSE)
    rcm <- rowToColMat(x, withDimnames=FALSE)
    crm <- colToRowMat(x, withDimnames=FALSE)

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(x), fmt
            )
        }
        i <- as.vector(i)
        rv[i] <- rowVec(value, withDimnames=FALSE)
        rrm[i,] <- rowToRowMat(value, withDimnames=FALSE)
        crm[,i] <- colToRowMat(value, withDimnames=FALSE)
    }

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                j, colnames(x), fmt
            )
        }
        j <- as.vector(j)
        cv[j] <- colVec(value, withDimnames=FALSE)
        ccm[,j] <- colToColMat(value, withDimnames=FALSE)
        rcm[j,] <- rowToColMat(value, withDimnames=FALSE)
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, rowVec=rv, colVec=cv,
        rowToRowMat=rrm, colToColMat=ccm, 
        rowToColMat=rcm, colToRowMat=crm, check=FALSE)
})

## -----------------------------------------------------------------------------
#' @export
setMethod("rbind", "ExampleClass", function(..., deparse.level=1) {
    args <- list(...)
    all.rv <- lapply(args, rowVec, withDimnames=FALSE)
    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)
    all.crm <- lapply(args, colToRowMat, withDimnames=FALSE)

    all.rv <- do.call(c, all.rv)
    all.rrm <- do.call(rbind, all.rrm)
    all.crm <- do.call(cbind, all.crm)

    # Checks for identical column state.
    ref <- args[[1]]
    ref.cv <- colVec(ref, withDimnames=FALSE)
    ref.ccm <- colToColMat(ref, withDimnames=FALSE)
    ref.rcm <- rowToColMat(ref, withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(ref.cv, colVec(x, withDimnames=FALSE)) 
            || !identical(ref.ccm, colToColMat(x, withDimnames=FALSE))
            || !identical(ref.rcm, rowToColMat(x, withDimnames=FALSE)))
        {
            stop("per-column values are not compatible")
        }
    }
 
    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, rowVec=all.rv,
        rowToRowMat=all.rrm, colToRowMat=all.crm, 
        check=FALSE)
})

## -----------------------------------------------------------------------------
#' @export
setMethod("cbind", "ExampleClass", function(..., deparse.level=1) {
    args <- list(...)
    all.cv <- lapply(args, colVec, withDimnames=FALSE)
    all.ccm <- lapply(args, colToColMat, withDimnames=FALSE)
    all.rcm <- lapply(args, rowToColMat, withDimnames=FALSE)

    all.cv <- do.call(c, all.cv)
    all.ccm <- do.call(cbind, all.ccm)
    all.rcm <- do.call(rbind, all.rcm)

    # Checks for identical column state.
    ref <- args[[1]]
    ref.rv <- rowVec(ref, withDimnames=FALSE)
    ref.rrm <- rowToRowMat(ref, withDimnames=FALSE)
    ref.crm <- colToRowMat(ref, withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(ref.rv, rowVec(x, withDimnames=FALSE)) 
            || !identical(ref.rrm, rowToRowMat(x, withDimnames=FALSE))
            || !identical(ref.crm, colToRowMat(x, withDimnames=FALSE)))
        {
            stop("per-row values are not compatible")
        }
    }

    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, colVec=all.cv,
        colToColMat=all.ccm, rowToColMat=all.rcm, 
        check=FALSE)
})

## -----------------------------------------------------------------------------
#' @exportMethods coerce
setAs("SummarizedExperiment", "ExampleClass", function(from) {
    new("ExampleClass", from, 
        rowVec=integer(nrow(from)), 
        colVec=integer(ncol(from)),
        rowToRowMat=matrix(0,nrow(from),0),
        colToColMat=matrix(0,0,ncol(from)),
        rowToColMat=matrix(0,ncol(from),0),
        colToRowMat=matrix(0,0,nrow(from)))
})

## -----------------------------------------------------------------------------
se <- SummarizedExperiment(matrix(rpois(100, lambda=1), ncol=5))
as(se, "CountSE")

## -----------------------------------------------------------------------------
RV <- 1:10
CV <- sample(50, 7)
RRM <- matrix(runif(30), nrow=10)
CCM <- matrix(rnorm(14), ncol=7)
RCM <- matrix(runif(21), nrow=7)
CRM <- matrix(rnorm(20), ncol=10)

thing <- ExampleClass(rowVec=RV, colVec=CV,
    rowToRowMat=RRM, colToColMat=CCM,
    rowToColMat=RCM, colToRowMat=CRM,
    assays=list(counts=matrix(rnorm(70), nrow=10)),
    colData=DataFrame(whee=LETTERS[1:7]),
    rowData=DataFrame(yay=letters[1:10])
)

## -----------------------------------------------------------------------------
rownames(thing) <- paste0("FEATURE_", seq_len(nrow(thing)))
colnames(thing) <- paste0("SAMPLE_", seq_len(ncol(thing)))
thing

## -----------------------------------------------------------------------------
expect_true(validObject(thing))

## -----------------------------------------------------------------------------
expect_true(validObject(.ExampleClass())) # internal
expect_true(validObject(ExampleClass())) # exported

## -----------------------------------------------------------------------------
expect_error(ExampleClass(rowVec=1), "rowVec")
expect_error(ExampleClass(colVec=1), "colVec")
expect_error(ExampleClass(rowToRowMat=rbind(1)), "rowToRowMat")
expect_error(ExampleClass(colToColMat=rbind(1)), "colToColMat")
expect_error(ExampleClass(rowToColMat=rbind(1)), "rowToColMat")
expect_error(ExampleClass(colToRowMat=rbind(1)), "colToRowMat")

## -----------------------------------------------------------------------------
se <- as(thing, "SummarizedExperiment")
conv <- as(se, "ExampleClass")
expect_true(validObject(conv))

## -----------------------------------------------------------------------------
expect_identical(names(rowVec(thing)), rownames(thing))
expect_identical(rowVec(thing, withDimnames=FALSE), RV)

expect_identical(names(colVec(thing)), colnames(thing))
expect_identical(colVec(thing, withDimnames=FALSE), CV)

## -----------------------------------------------------------------------------
expect_identical(rowToRowMat(thing, withDimnames=FALSE), RRM)
expect_identical(rownames(rowToRowMat(thing)), rownames(thing))

expect_identical(colToColMat(thing, withDimnames=FALSE), CCM)
expect_identical(colnames(colToColMat(thing)), colnames(thing))

expect_identical(rowToColMat(thing, withDimnames=FALSE), RCM)
expect_identical(rownames(rowToColMat(thing)), colnames(thing))

expect_identical(colToRowMat(thing, withDimnames=FALSE), CRM)
expect_identical(colnames(colToRowMat(thing)), rownames(thing))

## -----------------------------------------------------------------------------
expect_true("extra" %in% colnames(rowData(thing)))

## -----------------------------------------------------------------------------
rowVec(thing) <- 0:9
expect_equivalent(rowVec(thing), 0:9)

colVec(thing) <- 7:1
expect_equivalent(colVec(thing), 7:1)

## -----------------------------------------------------------------------------
old <- rowToRowMat(thing)
rowToRowMat(thing) <- -old
expect_equivalent(rowToRowMat(thing), -old)

old <- colToColMat(thing)
colToColMat(thing) <- 2 * old
expect_equivalent(colToColMat(thing), 2 * old)

old <- rowToColMat(thing)
rowToColMat(thing) <- old + 1
expect_equivalent(rowToColMat(thing), old + 1)

old <- colToRowMat(thing) 
colToRowMat(thing) <- old / 10
expect_equivalent(colToRowMat(thing), old / 10)

## -----------------------------------------------------------------------------
expect_message(rowData(thing) <- 1, "hi")

## -----------------------------------------------------------------------------
expect_error(rowVec(thing) <- 0, "rowVec")
expect_error(colVec(thing) <- 0, "colVec")
expect_error(rowToRowMat(thing) <- rbind(0), "rowToRowMat")
expect_error(colToColMat(thing) <- rbind(0), "colToColMat")
expect_error(rowToColMat(thing) <- rbind(0), "rowToColMat")
expect_error(colToRowMat(thing) <- rbind(0), "colToRowMat")

## -----------------------------------------------------------------------------
modified <- normalize(thing)
expect_equal(rowVec(modified), -rowVec(thing))
expect_equal(colVec(modified), -colVec(thing))

## -----------------------------------------------------------------------------
subbyrow <- thing[1:5,]
expect_identical(rowVec(subbyrow), rowVec(thing)[1:5])
expect_identical(rowToRowMat(subbyrow), rowToRowMat(thing)[1:5,])
expect_identical(colToRowMat(subbyrow), colToRowMat(thing)[,1:5])

# columns unaffected...
expect_identical(colVec(subbyrow), colVec(thing)) 
expect_identical(colToColMat(subbyrow), colToColMat(thing))
expect_identical(rowToColMat(subbyrow), rowToColMat(thing))

## -----------------------------------------------------------------------------
subbycol <- thing[,1:2]
expect_identical(colVec(subbycol), colVec(thing)[1:2])
expect_identical(colToColMat(subbycol), colToColMat(thing)[,1:2])
expect_identical(rowToColMat(subbycol), rowToColMat(thing)[1:2,])

# rows unaffected...
expect_identical(rowVec(subbycol), rowVec(thing)) 
expect_identical(rowToRowMat(subbycol), rowToRowMat(thing))
expect_identical(colToRowMat(subbycol), colToRowMat(thing))

## -----------------------------------------------------------------------------
norow <- thing[0,]
expect_true(validObject(norow))
expect_identical(nrow(norow), 0L)

nocol <- thing[,0]
expect_true(validObject(nocol))
expect_identical(ncol(nocol), 0L)

## -----------------------------------------------------------------------------
modified <- thing
modified[1:5,1:2] <- thing[5:1,2:1]

rperm <- c(5:1, 6:nrow(thing))
expect_identical(rowVec(modified), rowVec(thing)[rperm])
expect_identical(rowToRowMat(modified), rowToRowMat(thing)[rperm,])
expect_identical(colToRowMat(modified), colToRowMat(thing)[,rperm])

cperm <- c(2:1, 3:ncol(thing))
expect_identical(colVec(modified), colVec(thing)[cperm])
expect_identical(colToColMat(modified), colToColMat(thing)[,cperm])
expect_identical(rowToColMat(modified), rowToColMat(thing)[cperm,])

## -----------------------------------------------------------------------------
modified <- thing
modified[0,] <- thing[0,]
expect_equal(modified, thing)
modified[1,] <- thing[1,]
expect_equal(modified, thing)
modified[,0] <- thing[,0]
expect_equal(modified, thing)
modified[,1] <- thing[,1]
expect_equal(modified, thing)

## -----------------------------------------------------------------------------
expect_error(modified[1,1] <- thing[0,0], "replacement has length zero")

## -----------------------------------------------------------------------------
combined <- rbind(thing, thing)

rtwice <- rep(seq_len(nrow(thing)), 2)
expect_identical(rowVec(combined), rowVec(thing)[rtwice])
expect_identical(rowToRowMat(combined), rowToRowMat(thing)[rtwice,])
expect_identical(colToRowMat(combined), colToRowMat(thing)[,rtwice])

# Columns are unaffected:
expect_identical(colVec(combined), colVec(thing))
expect_identical(colToColMat(combined), colToColMat(thing))
expect_identical(rowToColMat(combined), rowToColMat(thing))

## -----------------------------------------------------------------------------
combined <- cbind(thing, thing)

ctwice <- rep(seq_len(ncol(thing)), 2)
expect_equivalent(colVec(combined), colVec(thing)[ctwice]) 
expect_equivalent(colToColMat(combined), colToColMat(thing)[,ctwice])
expect_equivalent(rowToColMat(combined), rowToColMat(thing)[ctwice,])

# Rows are unaffected:
expect_equivalent(rowVec(combined), rowVec(thing)) 
expect_equivalent(rowToRowMat(combined), rowToRowMat(thing))
expect_equivalent(colToRowMat(combined), colToRowMat(thing))

## -----------------------------------------------------------------------------
expect_equal(thing, rbind(thing))
expect_equal(thing, rbind(thing, thing[0,]))

expect_equal(thing, cbind(thing))
expect_equal(thing, cbind(thing, thing[,0]))

## -----------------------------------------------------------------------------
expect_error(rbind(thing, thing[,ncol(thing):1]), "not compatible")
expect_error(cbind(thing, thing[nrow(thing):1,]), "not compatible")

## -----------------------------------------------------------------------------
sessionInfo()

