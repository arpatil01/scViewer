## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----preliminaries, echo=FALSE, message=FALSE------------------------------
library(GenomeInfoDb)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

## ----genomeStyles1---------------------------------------------------------
seqmap <- genomeStyles()
head(seqmap,n=2)

## ----name------------------------------------------------------------------
names(genomeStyles())

## ----genomeStyles2---------------------------------------------------------
head(genomeStyles("Homo_sapiens"),5)

## ----style-present---------------------------------------------------------
"UCSC" %in% names(genomeStyles("Homo_sapiens"))

## ----extractSeqlevels------------------------------------------------------
extractSeqlevels(species="Arabidopsis_thaliana", style="NCBI")

## ----extractSeqlevelsgroup-------------------------------------------------
extractSeqlevelsByGroup(species="Arabidopsis_thaliana", style="NCBI",
                         group="auto")

## ----seqlevelsStyle--------------------------------------------------------
seqlevelsStyle(paste0("chr",c(1:30)))
seqlevelsStyle(c("2L","2R","X","Xhet"))

## ----keepChr-txdb----------------------------------------------------------
newchr <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9_hap1"))
seqlevelsInGroup(newchr, group="sex")
seqlevelsInGroup(newchr, group="auto")
seqlevelsInGroup(newchr, group="circular")
seqlevelsInGroup(newchr, group="sex","Homo_sapiens","UCSC")

## ----check2----------------------------------------------------------------
seqnames <- c("chr1", "chr9", "chr2", "chr3", "chr10")
all(seqnames %in% extractSeqlevels("Homo_sapiens", "UCSC"))

## ----orderSeqlevels--------------------------------------------------------
seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
orderSeqlevels(seqnames)
seqnames[orderSeqlevels(seqnames)]

## ----rankSeqlevels---------------------------------------------------------
seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
rankSeqlevels(seqnames)

## ----find------------------------------------------------------------------
mapSeqlevels(c("chrII", "chrIII", "chrM"), "NCBI")

## ----basic-gr--------------------------------------------------------------
gr <- GRanges(paste0("ch",1:35), IRanges(1:35, width=5))
gr

## ----renameseqlevels-------------------------------------------------------
newnames <- paste0("chr",1:35)
names(newnames) <- paste0("ch",1:35)
head(newnames)
gr <- renameSeqlevels(gr,newnames)
gr

## ----dropseqlevels---------------------------------------------------------
dropSeqlevels(gr, paste0("chr",23:35), pruning.mode="coarse")

## ----keepseqlevels---------------------------------------------------------
keepSeqlevels(gr, paste0("chr",1:22), pruning.mode="coarse")

## ----keepstdchr------------------------------------------------------------
keepStandardChromosomes(gr, pruning.mode="coarse")

## ----keepstdchr-2----------------------------------------------------------
plantgr <- GRanges(c(1:5,"MT","Pltd"), IRanges(1:7,width=5))
keepStandardChromosomes(plantgr, species="Arabidopsis thaliana",
                                 pruning.mode="coarse")

## ----Seqinfo-egs-----------------------------------------------------------
## Note that all the arguments (except 'genome') must have the
## same length. 'genome' can be of length 1, whatever the lengths
## of the other arguments are.
x <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
             seqlengths=c(100, 200, NA, 15),
             isCircular=c(NA, FALSE, FALSE, TRUE),
             genome="toy")
length(x)
seqnames(x)
names(x)
seqlevels(x)
seqlengths(x)
isCircular(x)
genome(x)

x[c("chrY", "chr3", "chr1")]  # subset by names

## Rename, drop, add and/or reorder the sequence levels:
xx <- x
seqlevels(xx) <- sub("chr", "ch", seqlevels(xx))  # rename
xx
seqlevels(xx) <- rev(seqlevels(xx))  # reorder
xx
seqlevels(xx) <- c("ch1", "ch2", "chY")  # drop/add/reorder
xx
seqlevels(xx) <- c(chY="Y", ch1="1", "22")  # rename/reorder/drop/add
xx

y <- Seqinfo(seqnames=c("chr3", "chr4", "chrM"),
             seqlengths=c(300, NA, 15))
y
merge(x, y)  # rows for chr3 and chrM are merged
suppressWarnings(merge(x, y))

## Note that, strictly speaking, merging 2 Seqinfo objects is not
## a commutative operation, i.e., in general 'z1 <- merge(x, y)'
## is not identical to 'z2 <- merge(y, x)'. However 'z1' and 'z2'
## are guaranteed to contain the same information (i.e. the same
## rows, but typically not in the same order):
suppressWarnings(merge(y, x))

## This contradicts what 'x' says about circularity of chr3 and chrM:
isCircular(y)[c("chr3", "chrM")] <- c(TRUE, FALSE)
y
if (interactive()) {
  merge(x, y)  # raises an error
}

## ----quick-style-----------------------------------------------------------
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
seqlevels(txdb)
genomeStyles("Drosophila melanogaster")
mapSeqlevels(seqlevels(txdb), "NCBI")

## ----sequence, eval=FALSE--------------------------------------------------
#  sequence <- seqlevels(x)
#  
#  ## sequence is in UCSC format and we want NCBI style
#  newStyle <- mapSeqlevels(sequence,"NCBI")
#  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
#  
#  ## rename the seqlevels
#  x <- renameSeqlevels(x,newStyle)
#  
#  ## keep only the seqlevels you want (say autosomes)
#  auto <- extractSeqlevelsByGroup(species="Homo sapiens", style="NCBI",
#                                  group="auto")
#  x <- keepSeqlevels(x,auto)

## ----sessionInfo, results='asis', eval=TRUE--------------------------------
toLatex(sessionInfo())

