## ----style, echo=FALSE, results='asis'----------------------------------------
BiocStyle::markdown()

## ----BiocManager, eval=FALSE--------------------------------------------------
#  if (!require("BiocManager"))
#      install.packages("BiocManager")
#  BiocManager::install("GenomicRanges")

## ----initialize, results="hide", warning=FALSE, message=FALSE-----------------
library(GenomicRanges)

## ----example-GRanges----------------------------------------------------------
gr <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    score = 1:10,
    GC = seq(1, 0, length=10))
gr

## ----GRanges-location-accessors-----------------------------------------------
seqnames(gr)
ranges(gr)
strand(gr)

## ----granges-accessor---------------------------------------------------------
granges(gr)

## ----metadataAccess-----------------------------------------------------------
mcols(gr)
mcols(gr)$score

## ----setSeqLengths------------------------------------------------------------
seqlengths(gr) <- c(249250621, 243199373, 198022430)

## ----setSeqLengths2-----------------------------------------------------------
seqlengths(gr)

## ----names--------------------------------------------------------------------
names(gr)
length(gr)

## ----splitAppendGRanges-------------------------------------------------------
sp <- split(gr, rep(1:2, each=5))
sp

## ----combine------------------------------------------------------------------
c(sp[[1]], sp[[2]])

## ----subset1------------------------------------------------------------------
gr[2:3]

## ----subset2------------------------------------------------------------------
gr[2:3, "GC"]

## ----assign1------------------------------------------------------------------
singles <- split(gr, names(gr))
grMod <- gr
grMod[2] <- singles[[1]]
head(grMod, n=3)

## ----other--------------------------------------------------------------------
rep(singles[[2]], times = 3)
rev(gr)
head(gr,n=2)
tail(gr,n=2)
window(gr, start=2,end=4)
gr[IRanges(start=c(2,7), end=c(3,9))]

## ----IRangesStuff-------------------------------------------------------------
g <- gr[1:3]
g <- append(g, singles[[10]])
start(g)
end(g)
width(g)
range(g)

## ----flank--------------------------------------------------------------------
flank(g, 10)

## ----flank2-------------------------------------------------------------------
flank(g, 10, start=FALSE)

## ----shiftAndResize-----------------------------------------------------------
shift(g, 5)
resize(g, 30)

## ----reduce-------------------------------------------------------------------
reduce(g)

## ----gaps---------------------------------------------------------------------
gaps(g)

## ----disjoin------------------------------------------------------------------
disjoin(g)

## ----coverage-----------------------------------------------------------------
coverage(g)

## ----intervals1---------------------------------------------------------------
g2 <- head(gr, n=2)
union(g, g2)
intersect(g, g2)
setdiff(g, g2)

## ----intervals2---------------------------------------------------------------
g3 <- g[1:2]
ranges(g3[1]) <- IRanges(start=105, end=112)
punion(g2, g3)
pintersect(g2, g3)
psetdiff(g2, g3)

## ----manPage, eval=FALSE------------------------------------------------------
#  ?GRanges

## ----granges-methods, eval=FALSE----------------------------------------------
#  methods(class="GRanges")

## ----example-GRangesList------------------------------------------------------
gr1 <- GRanges(
    seqnames = "chr2",
    ranges = IRanges(103, 106),
    strand = "+",
    score = 5L, GC = 0.45)
gr2 <- GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges(c(107, 113), width = 3),
    strand = c("+", "-"),
    score = 3:4, GC = c(0.3, 0.5))
grl <- GRangesList("txA" = gr1, "txB" = gr2)
grl

## ----basicGRLAccessors--------------------------------------------------------
seqnames(grl)
ranges(grl)
strand(grl)

## ----exceptions---------------------------------------------------------------
length(grl)
names(grl)
seqlengths(grl)

## ----elementNROWS-------------------------------------------------------------
elementNROWS(grl)

## ----isEmpty------------------------------------------------------------------
isEmpty(grl)

## ----mcolsGRL-----------------------------------------------------------------
mcols(grl) <- c("Transcript A","Transcript B")
mcols(grl)

## ----mcolsGRL-unlist----------------------------------------------------------
mcols(unlist(grl))

## ----unlistGRL----------------------------------------------------------------
ul <- unlist(grl)
ul

## ----pc-grl-------------------------------------------------------------------
grl1 <- GRangesList(
    gr1 = GRanges("chr2", IRanges(3, 6)),
    gr2 = GRanges("chr1", IRanges(c(7,13), width = 3)))
grl2 <- GRangesList(
    gr1 = GRanges("chr2", IRanges(9, 12)),
    gr2 = GRanges("chr1", IRanges(c(25,38), width = 3)))

pc(grl1, grl2)

grl3 <- c(grl1, grl2)
regroup(grl3, names(grl3))

## ----intOpsGRL----------------------------------------------------------------
start(grl)
end(grl)
width(grl)

## ----List-ops-----------------------------------------------------------------
sum(width(grl))  # sum of widths of each grl element

## ----coverageGRL--------------------------------------------------------------
shift(grl, 20)
coverage(grl)

## ----subsetGRL, eval=FALSE----------------------------------------------------
#  grl[1]
#  grl[[1]]
#  grl["txA"]
#  grl$txB

## ----subsetGRL2---------------------------------------------------------------
grl[1, "score"]
grl["txB", "GC"]

## ----otherSubsetGRL-----------------------------------------------------------
rep(grl[[1]], times = 3)
rev(grl)
head(grl, n=1)
tail(grl, n=1)
window(grl, start=1, end=1)
grl[IRanges(start=2, end=2)]

## ----lapply-------------------------------------------------------------------
lapply(grl, length)
sapply(grl, length)

## ----mapply-------------------------------------------------------------------
grl2 <- shift(grl, 10)
names(grl2) <- c("shiftTxA", "shiftTxB")

mapply(c, grl, grl2)
Map(c, grl, grl2)

## ----endoapply----------------------------------------------------------------
endoapply(grl, rev)
mendoapply(c, grl, grl2)

## ----ReduceGRL----------------------------------------------------------------
Reduce(c, grl)

## ----unlist-relist------------------------------------------------------------
gr <- unlist(grl)
gr$log_score <- log(gr$score)
grl <- relist(gr, grl)
grl

## ----manPage2, eval=FALSE-----------------------------------------------------
#  ?GRangesList
#  methods(class="GRangesList")   # _partial_ list

## ----findOverlaps-------------------------------------------------------------
findOverlaps(gr, grl)

## ----countOL------------------------------------------------------------------
countOverlaps(gr, grl)

## ----subsetByOverlaps---------------------------------------------------------
subsetByOverlaps(gr,grl)

## ----select-first-------------------------------------------------------------
findOverlaps(gr, grl, select="first")
findOverlaps(grl, gr, select="first")

## ----subjectx_declaration-----------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
broads <- GenomicFeatures::genes(txdb)
x <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    score = 1:10, GC = seq(1, 0, length=10))
subject <- broads[ seqnames(broads) %in% seqlevels(gr) ]

## ----nearest------------------------------------------------------------------
nearest(x, subject)
nearest(x)

## ----precede------------------------------------------------------------------
precede(x, subject)

## ----follow-------------------------------------------------------------------
follow(x, subject)

## ----nearestKNeighbors--------------------------------------------------------
nearestKNeighbors(x, subject)
nearestKNeighbors(x, subject, k=10)

nearestKNeighbors(x)
nearestKNeighbors(x, k=10)

## ----SessionInfo--------------------------------------------------------------
sessionInfo()

