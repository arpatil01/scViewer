## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----verbatim, message=FALSE-----------------------------------------------
library(GenomeInfoDb)
names(genomeStyles())

## ----email-----------------------------------------------------------------
packageDescription("GenomeInfoDb")$Maintainer

