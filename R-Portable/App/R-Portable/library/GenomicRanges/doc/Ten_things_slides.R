### R code from vignette source 'Ten_things_slides.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=80)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(BSgenome)
library(hgu95av2probe)

example(GRangesList)

gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
              IRanges(1:10, width=10:1, names=head(letters, 10)),
              Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
              score=1:10, GC=seq(1, 0, length=10))

ir <- IRanges(c(11:13, 2, 7:6), width=3)
mcols(ir) <- DataFrame(id=letters[1:6], score=3:-2)

x <- GRanges(c("chr1:1-1000", "chr2:2000-3000"),
             score=c(0.45, 0.1), a1=c(5L, 7L), a2=c(6, 8))
mcols(x)$score[2] <- NA
y <- GRanges(c("chr2:150-151", "chr1:1-10", "chr2:2000-3000"),
             score=c(0.7, 0.82, 0.1), b1=c(0L, 5L, 1L), b2=c(1, -2, 1))


###################################################
### code chunk number 2: inner_outer_mcols
###################################################
mcols(grl)$id <- paste0("ID", seq_along(grl))
grl


###################################################
### code chunk number 3: inner_outer_mcols
###################################################
mcols(grl)  # outer mcols
mcols(unlist(grl, use.names=FALSE))  # inner mcols


###################################################
### code chunk number 4: Ten_things_slides.Rnw:84-85
###################################################
gr


###################################################
### code chunk number 5: Ten_things_slides.Rnw:95-96
###################################################
invertStrand(gr)


###################################################
### code chunk number 6: Ten_things_slides.Rnw:106-107
###################################################
grl


###################################################
### code chunk number 7: Ten_things_slides.Rnw:117-118
###################################################
invertStrand(grl)


###################################################
### code chunk number 8: Ten_things_slides.Rnw:131-135
###################################################
cvg <- Rle(c(0L, 2L, 5L, 1L, 0L), c(10, 6, 3, 4, 15))
cvg
i <- IRanges(c(16, 19, 9), width=5, names=letters[1:3])
i


###################################################
### code chunk number 9: Ten_things_slides.Rnw:143-144
###################################################
extractList(cvg, i)


###################################################
### code chunk number 10: Ten_things_slides.Rnw:154-157
###################################################
i <- IntegerList(c(25:20), NULL, seq(from=2, to=length(cvg), by=2))
i
extractList(cvg, i)


###################################################
### code chunk number 11: Ten_things_slides.Rnw:168-171
###################################################
ir
ir2 <- reduce(ir, with.revmap=TRUE)
ir2


###################################################
### code chunk number 12: Ten_things_slides.Rnw:180-186
###################################################
revmap <- mcols(ir2)$revmap
extractList(mcols(ir)$id, revmap)
extractList(mcols(ir)$score, revmap)
mcols(ir2) <- DataFrame(id=extractList(mcols(ir)$id, revmap),
                        score=extractList(mcols(ir)$score, revmap))
ir2


###################################################
### code chunk number 13: Ten_things_slides.Rnw:199-202
###################################################
sliding_query <- IRanges(1:6, width=0)
sliding_query
countOverlaps(sliding_query, IRanges(3, 4))


###################################################
### code chunk number 14: Ten_things_slides.Rnw:209-210
###################################################
countOverlaps(sliding_query, IRanges(3, 4), minoverlap=0)


###################################################
### code chunk number 15: Ten_things_slides.Rnw:223-227
###################################################
library(Biostrings)
library(hgu95av2probe)
probes <- DNAStringSet(hgu95av2probe)
probes


###################################################
### code chunk number 16: Ten_things_slides.Rnw:236-237
###################################################
replaceAt(probes, at=IRanges(3, 4), value="-++-")


###################################################
### code chunk number 17: Ten_things_slides.Rnw:246-247
###################################################
replaceAt(probes, at=IRanges(3, 4), value="")


###################################################
### code chunk number 18: Ten_things_slides.Rnw:256-257
###################################################
replaceAt(probes, at=IRanges(4, 3), value="-++-")


###################################################
### code chunk number 19: Ten_things_slides.Rnw:267-269
###################################################
midx <- vmatchPattern("VCGTT", probes, fixed=FALSE)
replaceAt(probes, at=midx, value="-++-")


###################################################
### code chunk number 20: GRanges_as_a_subscript_1
###################################################
cvg <- RleList(chr1=101:120, chr2=2:-8, chr3=31:40)
gr


###################################################
### code chunk number 21: GRanges_as_a_subscript_2
###################################################
cvg[gr]


###################################################
### code chunk number 22: Ten_things_slides.Rnw:304-310
###################################################
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
ex <- exons(txdb, columns=c("exon_id", "tx_name", "gene_id"))
v <- Views(genome, ex)


###################################################
### code chunk number 23: Ten_things_slides.Rnw:319-320
###################################################
v


###################################################
### code chunk number 24: Ten_things_slides.Rnw:329-331
###################################################
af <- alphabetFrequency(v, baseOnly=TRUE)
head(af)


###################################################
### code chunk number 25: Ten_things_slides.Rnw:341-351
###################################################
library(Rsamtools)
library(RNAseqData.HNRNPC.bam.chr14)
fl <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]
sbp <- ScanBamParam(which=GRanges("chr14", IRanges(1, 53674770)))
pp <- PileupParam(distinguish_nucleotides=FALSE,
                  distinguish_strands=FALSE,
                  min_mapq=13,
                  min_base_quality=10,
                  min_nucleotide_depth=4)
res <- pileup(fl, scanBamParam=sbp, pileupParam=pp)


###################################################
### code chunk number 26: Ten_things_slides.Rnw:359-361
###################################################
dim(res)
head(res)


###################################################
### code chunk number 27: Ten_things_slides.Rnw:372-374
###################################################
x
y


###################################################
### code chunk number 28: Ten_things_slides.Rnw:384-385
###################################################
merge(x, y)


###################################################
### code chunk number 29: Ten_things_slides.Rnw:395-396
###################################################
merge(x, y, all=TRUE)


