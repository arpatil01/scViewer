### R code from vignette source 'ExtendingGenomicRanges.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: options
###################################################
options(width=72)
options(showHeadLines=3)
options(showTailLines=3)


###################################################
### code chunk number 3: granges-ranges
###################################################
library(GenomicRanges)
selectMethod(ranges, "GRanges")


###################################################
### code chunk number 4: delegating-granges-ranges
###################################################
selectMethod(ranges, "DelegatingGenomicRanges")


###################################################
### code chunk number 5: gnclist-granges
###################################################
getSlots("GNCList")["granges"]


###################################################
### code chunk number 6: vranges
###################################################
GenomicRanges:::extraColumnSlotNames(VariantAnnotation:::VRanges())


