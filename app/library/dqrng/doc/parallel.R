## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  library(parallel)
#  cl <- parallel::makeCluster(2)
#  res <- clusterApply(cl, 1:8, function(stream, seed, N) {
#    library(dqrng)
#    dqRNGkind("Threefry")
#    dqset.seed(seed, stream)
#    dqrnorm(N)
#  }, 42, 1e6)
#  stopCluster(cl)
#  
#  res <- matrix(unlist(res), ncol = 8)
#  symnum(x = cor(res), cutpoints = c(0.001, 0.003, 0.999),
#         symbols = c(" ", "?", "!", "1"),
#         abbr.colnames = FALSE, corr = TRUE)

