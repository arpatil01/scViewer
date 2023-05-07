## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
N <- 1e7
piR <- function(n, rng = runif) {
    x <- rng(n)
    y <- rng(n)
    4 * sum(sqrt(x^2 + y^2) < 1.0) / n
}
set.seed(42)
system.time(cat("pi ~= ", piR(N), "\n"))

## -----------------------------------------------------------------------------
library(dqrng)
dqRNGkind("Xoroshiro128+")
dqset.seed(42)
system.time(cat("pi ~= ", piR(N, rng = dqrunif), "\n"))

## -----------------------------------------------------------------------------
system.time(runif(N))
system.time(dqrunif(N))

## -----------------------------------------------------------------------------
system.time(rexp(N))
system.time(dqrexp(N))

## -----------------------------------------------------------------------------
system.time(rnorm(N))
system.time(dqrnorm(N))

## -----------------------------------------------------------------------------
system.time(for (i in 1:100)   sample.int(N, N/100, replace = TRUE))
system.time(for (i in 1:100) dqsample.int(N, N/100, replace = TRUE))
system.time(for (i in 1:100)   sample.int(N, N/100))
system.time(for (i in 1:100) dqsample.int(N, N/100))

