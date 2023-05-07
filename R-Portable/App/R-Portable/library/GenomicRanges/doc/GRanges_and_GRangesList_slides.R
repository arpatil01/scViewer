### R code from vignette source 'GRanges_and_GRangesList_slides.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=84)
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...)
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  par(mai=c(0.5, 0.2, 1.2, 0.2))
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main, cex.main=2.8, font.main=1)
  axis(1)
}


###################################################
### code chunk number 2: GRanges_constructor
###################################################
library(GenomicRanges)
gr1 <- GRanges(seqnames=Rle(c("ch1", "chMT"), c(2, 4)),
               ranges=IRanges(16:21, 20),
               strand=rep(c("+", "-", "*"), 2))
gr1


###################################################
### code chunk number 3: GRanges_accessors1
###################################################
length(gr1)
seqnames(gr1)
ranges(gr1)


###################################################
### code chunk number 4: GRanges_accessors2
###################################################
start(gr1)
end(gr1)
width(gr1)
strand(gr1)
strand(gr1) <- c("-", "-", "+")
strand(gr1)


###################################################
### code chunk number 5: GRanges_accessors3
###################################################
names(gr1) <- LETTERS[1:6]
gr1
names(gr1)


###################################################
### code chunk number 6: GRanges_accessors4
###################################################
mcols(gr1) <- DataFrame(score=11:16, GC=seq(1, 0, length=6))
gr1
mcols(gr1)


###################################################
### code chunk number 7: GRanges_accessors5
###################################################
seqinfo(gr1)
seqlevels(gr1)
seqlengths(gr1)
seqlengths(gr1) <- c(50000, 800)
seqlengths(gr1)


###################################################
### code chunk number 8: GRanges_Vector_ops1
###################################################
gr1[c("F", "A")]
gr1[strand(gr1) == "+"]


###################################################
### code chunk number 9: GRanges_Vector_ops2
###################################################
gr1 <- gr1[-5]
gr1


###################################################
### code chunk number 10: GRanges_Vector_ops3
###################################################
gr2 <- GRanges(seqnames="ch2",
               ranges=IRanges(start=c(2:1,2), width=6),
               score=15:13,
               GC=seq(0, 0.4, length=3))
gr12 <- c(gr1, gr2)
gr12


###################################################
### code chunk number 11: GRanges_Vector_ops4
###################################################
gr12[length(gr12)] == gr12
duplicated(gr12)
unique(gr12)


###################################################
### code chunk number 12: GRanges_sort
###################################################
sort(gr12)


###################################################
### code chunk number 13: GRanges_split
###################################################
split(gr12, seqnames(gr12))


###################################################
### code chunk number 14: ranges-ir0-plot
###################################################
library(IRanges)
ir0 <- IRanges(start=c(7, 9, 12, 14, 22:24),
               end=c(15, 11, 12, 18, 26, 27, 28))

png("ranges-ir0-plot.png", width=800, height=170)
plotRanges(ir0, xlim=c(5, 35), main="ir0", col="blue")
dev.off()


###################################################
### code chunk number 15: ranges-shift-ir0-plot
###################################################
png("ranges-shift-ir0-plot.png", width=800, height=170)
plotRanges(shift(ir0, 5), xlim=c(5, 35), main="shift(ir0, 5)", col="blue")
dev.off()


###################################################
### code chunk number 16: ranges-reduce-ir0-plot
###################################################
png("ranges-reduce-ir0-plot.png", width=800, height=170)
plotRanges(reduce(ir0), xlim=c(5, 35), main="reduce(ir0)", col="blue")
dev.off()


###################################################
### code chunk number 17: ranges-disjoin-ir0-plot
###################################################
png("ranges-disjoin-ir0-plot.png", width=800, height=170)
plotRanges(disjoin(ir0), xlim=c(5, 35), main="disjoin(ir0)", col="blue")
dev.off()


###################################################
### code chunk number 18: GRanges_range_based_ops1
###################################################
gr2
shift(gr2, 50)


###################################################
### code chunk number 19: GRanges_range_based_ops2
###################################################
gr1
resize(gr1, 12)


###################################################
### code chunk number 20: GRanges_range_based_ops3
###################################################
gr1
flank(gr1, 3)


###################################################
### code chunk number 21: GRanges_range_based_ops4
###################################################
gr3 <- shift(gr1, c(35000, rep(0, 3), 100))
width(gr3)[c(3,5)] <- 117
gr3
range(gr3)


###################################################
### code chunk number 22: GRanges_reduce
###################################################
gr3
reduce(gr3)


###################################################
### code chunk number 23: GRanges_gaps
###################################################
gr3
gaps(gr3)


###################################################
### code chunk number 24: GRanges_disjoin
###################################################
gr3
disjoin(gr3)


###################################################
### code chunk number 25: GRanges_coverage1
###################################################
cvg12 <- coverage(gr12)
cvg12


###################################################
### code chunk number 26: GRanges_coverage2
###################################################
mean(cvg12)
max(cvg12)


###################################################
### code chunk number 27: slice_coverage
###################################################
sl12 <- slice(cvg12, lower=1)
sl12
elementNROWS(sl12)
sl12$chMT
mean(sl12$chMT)
max(sl12$chMT)


###################################################
### code chunk number 28: findOverlaps1
###################################################
library(pasillaBamSubset)
untreated1_chr4()
library(GenomicAlignments)
reads <- readGAlignments(untreated1_chr4())


###################################################
### code chunk number 29: findOverlaps2
###################################################
reads <- as(reads, "GRanges")
reads[1:4]


###################################################
### code chunk number 30: findOverlaps3
###################################################
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
dm3_genes <- genes(txdb)


###################################################
### code chunk number 31: findOverlaps4
###################################################
hits <- findOverlaps(reads, dm3_genes)
head(hits)


###################################################
### code chunk number 32: GRangesList_constructor
###################################################
grl <- GRangesList(gr3, gr2)
grl


###################################################
### code chunk number 33: GRangesList_accessors1
###################################################
length(grl)


###################################################
### code chunk number 34: GRangesList_accessors2
###################################################
seqnames(grl)


###################################################
### code chunk number 35: GRangesList_accessors3
###################################################
strand(grl)


###################################################
### code chunk number 36: GRangesList_accessors4
###################################################
ranges(grl)


###################################################
### code chunk number 37: GRangesList_accessors5
###################################################
start(grl)
end(grl)
width(grl)


###################################################
### code chunk number 38: GRangesList_accessors6
###################################################
names(grl) <- c("TX1", "TX2")
grl


###################################################
### code chunk number 39: GRangesList_accessors7
###################################################
mcols(grl)$geneid <- c("GENE1", "GENE2") 
mcols(grl)
grl


###################################################
### code chunk number 40: GRangesList_accessors8
###################################################
seqinfo(grl)


###################################################
### code chunk number 41: GRangesList_Vector_ops1
###################################################
grl[c("TX2", "TX1")]


###################################################
### code chunk number 42: GRangesList_Vector_ops2
###################################################
c(grl, GRangesList(gr3))


###################################################
### code chunk number 43: GRangesList_List_ops1
###################################################
grl[[2]]
elementNROWS(grl)
unlisted <- unlist(grl, use.names=FALSE)  # same as c(grl[[1]], grl[[2]])
unlisted


###################################################
### code chunk number 44: GRangesList_List_ops2
###################################################
grl100 <- relist(shift(unlisted, 100), grl)
grl100


###################################################
### code chunk number 45: GRangesList_List_ops3
###################################################
grl100b <- endoapply(grl, shift, 100)
grl100b
mcols(grl100)
mcols(grl100b)


###################################################
### code chunk number 46: GRangesList_range_based_ops1
###################################################
grl


###################################################
### code chunk number 47: GRangesList_range_based_ops2
###################################################
shift(grl, 100)


###################################################
### code chunk number 48: GRangesList_range_based_ops3
###################################################
grl


###################################################
### code chunk number 49: GRangesList_range_based_ops4
###################################################
flank(grl, 10)


###################################################
### code chunk number 50: GRangesList_range_based_ops5
###################################################
grl


###################################################
### code chunk number 51: GRangesList_range_based_ops6
###################################################
range(grl) 


###################################################
### code chunk number 52: GRangesList_range_based_ops7
###################################################
grl


###################################################
### code chunk number 53: GRangesList_range_based_ops8
###################################################
reduce(grl) 


###################################################
### code chunk number 54: GRangesList_range_based_ops9
###################################################
grl2 <- grl
grl2[[1]] <- grl2[[1]][3]; grl2[[2]] <- grl2[[2]][1]
grl3 <- unname(grl2)
grl3[[1]] <- narrow(unname(grl3[[1]]), start=5, end=-5)


###################################################
### code chunk number 55: GRangesList_range_based_ops10
###################################################
grl2
grl3


###################################################
### code chunk number 56: GRangesList_range_based_ops11
###################################################
setdiff(grl2, grl3)


