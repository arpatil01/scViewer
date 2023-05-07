## ----message=FALSE, echo=FALSE------------------------------------------------
library(Biobase)
data(sample.ExpressionSet)

## -----------------------------------------------------------------------------
sample.ExpressionSet
exprs(sample.ExpressionSet)[1,]
pData(sample.ExpressionSet)[1:2,1:3]

## -----------------------------------------------------------------------------
rbind(exprs(sample.ExpressionSet[1,]),
            sex <- t(pData(sample.ExpressionSet))[1,])

## -----------------------------------------------------------------------------
medContr <- function( y, x ) {
    ys <- split(y,x)
    median(ys[[1]]) - median(ys[[2]])
 }

## -----------------------------------------------------------------------------
apply(exprs(sample.ExpressionSet[1,,drop=F]), 1, medContr, pData(sample.ExpressionSet)[["sex"]])

## -----------------------------------------------------------------------------
medContr1 <- function(y) {
   ys <- split(y,sex)
   median(ys[[1]]) - median(ys[[2]])
}
esApply( sample.ExpressionSet, 1, medContr1)[1]

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

