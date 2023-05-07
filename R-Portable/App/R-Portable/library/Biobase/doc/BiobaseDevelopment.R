## ----echo=FALSE, message=FALSE------------------------------------------------
library(Biobase)

## -----------------------------------------------------------------------------
getClass("eSet")

## -----------------------------------------------------------------------------
getValidity(getClass("eSet"))

## ----eval=FALSE---------------------------------------------------------------
#  obj <- new("ExpressionSet",
#             phenoData = new("AnnotatedDataFrame"),
#             experimentData = new("MIAME"), annotation = character(),
#             exprs = new("matrix"))

## ----eval=FALSE---------------------------------------------------------------
#  assayDataNew("environment", elt)

## ----warning=FALSE------------------------------------------------------------
data(sample.ExpressionSet)
storageMode(sample.ExpressionSet) 
tryCatch(assayData(sample.ExpressionSet)$exprs <- log(exprs(sample.ExpressionSet)), 
    error=function(err) cat(conditionMessage(err))) 

## ----eval=FALSE---------------------------------------------------------------
#  exprs(sample.ExpressionSet) <- log(exprs(sample.ExpressionSet))`

## -----------------------------------------------------------------------------
getClass("ExpressionSet") 
getValidity(getClass("ExpressionSet")) 

## -----------------------------------------------------------------------------
setClass("SwirlSet", contains="eSet") 

## ----message=FALSE------------------------------------------------------------
setMethod("initialize", "SwirlSet",
          function(.Object,
                   R = new("matrix"),
                   G = new("matrix"),
                   Rb = new("matrix"),
                   Gb = new("matrix"),
                   ...) {
            callNextMethod(.Object,
                           R=R, G=G, Rb=Rb, Gb=Gb,
                           ...) 
            }) 

## -----------------------------------------------------------------------------
setMethod("initialize", "SwirlSet",
          function(.Object,
                   assayData=assayDataNew(
                   R=R, G=G, Rb=Rb, Gb=Gb),
                   R = new("matrix"),
                   G = new("matrix"),
                   Rb = new("matrix"),
                   Gb = new("matrix"),
                   ...) {
              if (!missing(assayData) &&
                  any(!missing(R), !missing(G), !missing(Rb), !missing(Gb))) {
                  warning("using 'assayData'; ignoring 'R', 'G', 'Rb', 'Gb'")
              }
              callNextMethod(.Object, assayData=assayData, ...)
          })

## ----eval=FALSE---------------------------------------------------------------
#  new("SwirlSet")

## ----eval=FALSE---------------------------------------------------------------
#  setMethod("initialize", "MySet",
#            function(.Object, ...) {
#                .Object <- callNextMethod(.Object, ...)
#            })

## -----------------------------------------------------------------------------
setValidity("SwirlSet", function(object) {
    assayDataValidMembers(assayData(object), c("R", "G", "Rb", "Gb"))
})

## ----eval=FALSE---------------------------------------------------------------
#  myFancyFunction <- function(obj) {
#      assayData(obj) <- fancyAssaydData # obj invalid...
#      phenoData(obj) <- justAsFancyPhenoData # but now valid
#      validObject(obj)
#      (obj)
#  }

## -----------------------------------------------------------------------------
data(sample.ExpressionSet) 
classVersion(sample.ExpressionSet) 
obj <- updateObject(sample.ExpressionSet)

## -----------------------------------------------------------------------------
isCurrent(sample.ExpressionSet)[c("eSet", "ExpressionSet")]

## -----------------------------------------------------------------------------
setClass("MySet",
         contains = "eSet",
         prototype = prototype(new("VersionedBiobase",
                               versions=c(classVersion("eSet"), MySet="1.0.0"))))
obj <- new("MySet")
classVersion(obj)

## -----------------------------------------------------------------------------
 setClass("MySet",
          contains = "eSet",
          prototype = prototype( 
              new("VersionedBiobase",
                  versions=c(classVersion("eSet"), MySet="1.0.1"))))
isCurrent(obj)

## -----------------------------------------------------------------------------
setMethod("updateObject", signature(object="MySet"),
          function(object, ..., verbose=FALSE) {
              if (verbose) message("updateObject(object = 'MySet')")
                  object <- callNextMethod()
              if (isCurrent(object) ["MySet"]) return(object)
              ## Create an updated instance.
              if (!isVersioned(object))
                  ## Radical surgery – create a new, up-to-date instance
                  new("MySet",
                      assayData = updateObject(assayData(object),
                          ..., verbose=verbose),
                      phenoData = updateObject(phenoData(object),
                          ..., verbose=verbose),
                      experimentData = updateObject(experimentData(object),
                          ..., verbose=verbose),
                      annotation = updateObject(annotation(object),
                          ..., verbose=verbose))
              else {
                  ## Make minor changes, and update version by consulting class definition
                  classVersion(object)["MySet"]<-
                  classVersion("MySet")["MySet"]
                  object
              }
          })

## -----------------------------------------------------------------------------
classVersion(updateObject(obj))

## -----------------------------------------------------------------------------
classVersion(new("AnnotatedDataFrame"))

## -----------------------------------------------------------------------------
setClass("SwirlSet", contains = "eSet",
         prototype = prototype(
             new("VersionedBiobase",
                 versions=c(classVersion("eSet"), SwirlSet="1.0.0"))))
classVersion(new("SwirlSet"))

## -----------------------------------------------------------------------------
obj <- new("SwirlSet")
classVersion(obj)["MyID"] <- "0.0.1"
classVersion(obj)

## -----------------------------------------------------------------------------
 classVersion(updateObject(obj))

## -----------------------------------------------------------------------------
sessionInfo()

