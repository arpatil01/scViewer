## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)

## -----------------------------------------------------------------------------
system.file("extensions", package="beachmat")

## -----------------------------------------------------------------------------
beachmat_AaronMatrix_integer_output <- TRUE
beachmat_AaronMatrix_character_output <- TRUE

## ---- eval=FALSE--------------------------------------------------------------
#  testpkg <- system.file("testpkg", package="beachmat")
#  devtools::install(testpkg, quick=TRUE)
#  library(beachtest)

