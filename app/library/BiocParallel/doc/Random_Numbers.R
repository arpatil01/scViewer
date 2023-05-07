## -----------------------------------------------------------------------------
library(BiocParallel)
stopifnot(
    packageVersion("BiocParallel") > "1.27.5"
)

## -----------------------------------------------------------------------------
result1 <- bplapply(1:3, runif, BPPARAM = SerialParam(RNGseed = 100))
result1

## -----------------------------------------------------------------------------
result2 <- bplapply(1:3, runif, BPPARAM = SerialParam(RNGseed = 100))
stopifnot(
    identical(result1, result2)
)

result3 <- bplapply(1:3, runif, BPPARAM = SerialParam(RNGseed = 200))
result3

stopifnot(
    !identical(result1, result3)
)

## -----------------------------------------------------------------------------
result4 <- bplapply(1:3, runif, BPPARAM = SnowParam(RNGseed = 100))
stopifnot(
    identical(result1, result4)
)

if (!identical(.Platform$OS.type, "windows")) {
    result5 <- bplapply(1:3, runif, BPPARAM = MulticoreParam(RNGseed = 100))
    stopifnot(
        identical(result1, result5)
    )
}

## -----------------------------------------------------------------------------
result6 <- bplapply(1:3, runif, BPPARAM = SnowParam(workers = 2, RNGseed = 100))
result7 <- bplapply(1:3, runif, BPPARAM = SnowParam(workers = 3, RNGseed = 100))
result8 <- bplapply(
    1:3, runif,
    BPPARAM = SnowParam(workers = 2, tasks = 3, RNGseed = 100)
)
stopifnot(
    identical(result1, result6),
    identical(result1, result7),
    identical(result1, result8)
)

## -----------------------------------------------------------------------------
ITER_FUN_FACTORY <- function() {
    x <- 1:3
    i <- 0L
    function() {
        i <<- i + 1L
        if (i > length(x))
            return(NULL)
        x[[i]]
    }
}

## ---- collapse = TRUE---------------------------------------------------------
ITER <- ITER_FUN_FACTORY()
ITER()

ITER()

ITER()

ITER()

## -----------------------------------------------------------------------------
result9 <- bpiterate(
    ITER_FUN_FACTORY(), runif,
    BPPARAM = SerialParam(RNGseed = 100)
)
stopifnot(
    identical(result1, result9)
)

## -----------------------------------------------------------------------------
FUN1 <- function(i) {
    if (identical(i, 2L)) {
        ## error when evaluating the second element
        stop("i == 2")
    } else runif(i)
}
result10 <- bptry(bplapply(
    1:3, FUN1,
    BPPARAM = SerialParam(RNGseed = 100, stop.on.error = FALSE)
))
result10

## -----------------------------------------------------------------------------
FUN2 <- function(i) {
    if (identical(i, 2L)) {
        ## the random number stream should be in the same state as the
        ## first time through the loop, and rnorm(i) should return
        ## same result as FUN
        runif(i)
    } else {
        ## if this branch is used, then we are incorrectly updating
        ## already calculated elements -- '0' in the output would
        ## indicate this error
        0
    }
}
result11 <- bplapply(
    1:3, FUN2,
    BPREDO = result10,
    BPPARAM = SerialParam(RNGseed = 100, stop.on.error = FALSE)
)
stopifnot(
    identical(result1, result11)
)

## -----------------------------------------------------------------------------
set.seed(200)
value <- runif(1)

set.seed(200)
result12 <- bplapply(1:3, runif, BPPARAM = SerialParam(RNGseed = 100))
stopifnot(
    identical(result1, result12),
    identical(value, runif(1))
)

## -----------------------------------------------------------------------------
set.seed(100)
value <- runif(1)

set.seed(100)
result13 <- bplapply(1:3, runif, BPPARAM = SerialParam())
stopifnot(
    !identical(result1, result13),
    identical(value, runif(1))
)

## -----------------------------------------------------------------------------
param <- bpstart(SerialParam(RNGseed = 100))
result16 <- bplapply(1:3, runif, BPPARAM = param)
bpstop(param)
stopifnot(
    identical(result1, result16)
)

## -----------------------------------------------------------------------------
param <- bpstart(SerialParam(RNGseed = 100))
result16 <- bplapply(1:3, runif, BPPARAM = param)
result17 <- bplapply(1:3, runif, BPPARAM = param)
bpstop(param)
stopifnot(
    identical(result1, result16),
    !identical(result1, result17)
)

## -----------------------------------------------------------------------------
set.seed(100)
result20 <- lapply(1:3, runif)
stopifnot(
    !identical(result1, result20)
)

## ---- echo = FALSE------------------------------------------------------------
sessionInfo()

