### R code from vignette source 'GenomicRangesHOWTOs.Rnw'

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
.precomputed_results_path <- "precomputed_results"


###################################################
### code chunk number 3: pasillaBamSubset
###################################################
library(pasillaBamSubset)
untreated1_chr4()
untreated3_chr4()


###################################################
### code chunk number 4: pasillaBamSubset_help (eval = FALSE)
###################################################
## ?pasillaBamSubset


###################################################
### code chunk number 5: readGAlignments_1
###################################################
library(pasillaBamSubset)
un1 <- untreated1_chr4()  # single-end reads


###################################################
### code chunk number 6: readGAlignments_2
###################################################
library(GenomicAlignments)
gal <- readGAlignments(un1)


###################################################
### code chunk number 7: readGAlignments_3
###################################################
what <- c("flag", "cigar") 
which <- GRanges("chr4", IRanges(1, 5000)) 
flag <- scanBamFlag(isMinusStrand = TRUE)
param <- ScanBamParam(which=which, what=what, flag=flag)
neg <- readGAlignments(un1, param=param)
neg


###################################################
### code chunk number 8: readGAlignmentPairs_1
###################################################
library(pasillaBamSubset)
un3 <- untreated3_chr4()  # paired-end reads


###################################################
### code chunk number 9: readGAlignmentPairs_2
###################################################
un3 <- untreated3_chr4()
gapairs <- readGAlignmentPairs(un3)


###################################################
### code chunk number 10: readGAlignmentPairs_3
###################################################
gapairs


###################################################
### code chunk number 11: readGAlignmentsList_1
###################################################
galist <- readGAlignmentsList(BamFile(un3, asMates=TRUE))


###################################################
### code chunk number 12: readGAlignmentsList_2
###################################################
galist


###################################################
### code chunk number 13: readGAlignmentsList_3
###################################################
non_mates <- galist[unlist(mcols(galist)$mate_status) == "unmated"]
table(elementNROWS(non_mates))


###################################################
### code chunk number 14: yieldSize
###################################################
library(pasillaBamSubset)
un1 <- untreated1_chr4()
bf <- BamFile(un1, yieldSize=100000)


###################################################
### code chunk number 15: readGAlignments_by_chunk
###################################################
open(bf)
cvg <- NULL
repeat {
    chunk <- readGAlignments(bf)
    if (length(chunk) == 0L)
        break
    chunk_cvg <- coverage(chunk)
    if (is.null(cvg)) {
        cvg <- chunk_cvg
    } else {
        cvg <- cvg + chunk_cvg
    }
}
close(bf)
cvg


###################################################
### code chunk number 16: coverage_1
###################################################
library(pasillaBamSubset)
un1 <- untreated1_chr4()  # single-end reads
library(GenomicAlignments)
reads1 <- readGAlignments(un1)
cvg1 <- coverage(reads1)
cvg1


###################################################
### code chunk number 17: coverage_2
###################################################
cvg1$chr4


###################################################
### code chunk number 18: coverage_3
###################################################
mean(cvg1$chr4)
max(cvg1$chr4)


###################################################
### code chunk number 19: peaks_1
###################################################
chr4_peaks <- slice(cvg1$chr4, lower=500)
chr4_peaks
length(chr4_peaks)  # nb of peaks


###################################################
### code chunk number 20: peaks_2
###################################################
sum(chr4_peaks)


###################################################
### code chunk number 21: makeTxDbFromUCSC_1 (eval = FALSE)
###################################################
## library(GenomicFeatures)
## ### Internet connection required! Can take several minutes...
## txdb <- makeTxDbFromUCSC(genome="sacCer2", tablename="ensGene")


###################################################
### code chunk number 22: TxDb.Hsapiens.UCSC.hg19.knownGene_1
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb


###################################################
### code chunk number 23: TxDb.Hsapiens.UCSC.hg19.knownGene_2
###################################################
transcripts(txdb)


###################################################
### code chunk number 24: makeTxDbFromBiomart_1 (eval = FALSE)
###################################################
## library(GenomicFeatures)
## ### Internet connection required! Can take several minutes...
## txdb <- makeTxDbFromBiomart(biomart="ensembl",
##                             dataset="hsapiens_gene_ensembl")


###################################################
### code chunk number 25: TxDb.Athaliana.BioMart.plantsmart22_1
###################################################
library(TxDb.Athaliana.BioMart.plantsmart22)
txdb <- TxDb.Athaliana.BioMart.plantsmart22
txdb


###################################################
### code chunk number 26: TxDb.Athaliana.BioMart.plantsmart22_2
###################################################
exons(txdb)


###################################################
### code chunk number 27: makeTxDbFromGFF_1
###################################################
library(GenomicFeatures)
gff_file <- system.file("extdata", "GFF3_files", "a.gff3",
                        package="GenomicFeatures")
txdb <- makeTxDbFromGFF(gff_file, format="gff3")
txdb


###################################################
### code chunk number 28: makeTxDbFromGFF_2
###################################################
exonsBy(txdb, by="gene")


###################################################
### code chunk number 29: hub_1
###################################################
library(AnnotationHub)
### Internet connection required!
hub <- AnnotationHub()
hub <- subset(hub, hub$species=='Drosophila melanogaster')


###################################################
### code chunk number 30: hub_2
###################################################
length(hub)
head(names(hub))
head(hub$title, n=10)
## then look at a specific slice of the hub object.
hub[7]


###################################################
### code chunk number 31: hub_3
###################################################
gr <- hub[[names(hub)[7]]]
summary(gr)


###################################################
### code chunk number 32: hub_4
###################################################
metadata(gr)


###################################################
### code chunk number 33: hub_5
###################################################
txbygn <- split(gr, gr$name)


###################################################
### code chunk number 34: count_1
###################################################
library(pasillaBamSubset)
reads <- c(untrt1=untreated1_chr4(),  # single-end reads
           untrt3=untreated3_chr4())  # paired-end reads


###################################################
### code chunk number 35: count_2
###################################################
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
exbygene <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")


###################################################
### code chunk number 36: count_3
###################################################
library(GenomicAlignments)
se <- summarizeOverlaps(exbygene, reads, mode="IntersectionNotEmpty")


###################################################
### code chunk number 37: count_4
###################################################
class(se)
head(table(assays(se)$counts))


###################################################
### code chunk number 38: count_5
###################################################
identical(length(exbygene), length(assays(se)$counts))


###################################################
### code chunk number 39: count_6
###################################################
rowRanges(se)


###################################################
### code chunk number 40: count_table
###################################################
colData(se)$trt <- factor(c("untrt", "untrt"), levels=c("trt", "untrt"))
colData(se)

library(DESeq2)
deseq <- DESeqDataSet(se, design= ~ 1)

library(edgeR)
edger <- DGEList(assays(se)$counts, group=rownames(colData(se)))


###################################################
### code chunk number 41: summarizeJunctions_1
###################################################
library(pasillaBamSubset)
un1 <- untreated1_chr4()  # single-end reads
library(GenomicAlignments)
reads1 <- readGAlignments(un1)
reads1


###################################################
### code chunk number 42: summarizeJunctions_2
###################################################
junc_summary <- summarizeJunctions(reads1)
junc_summary


###################################################
### code chunk number 43: trak_1
###################################################
trak2 <- "66008"


###################################################
### code chunk number 44: trak_2
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


###################################################
### code chunk number 45: trak_3
###################################################
library(GenomicFeatures)
trak2_txs <- transcriptsBy(txdb, by="gene")[[trak2]]
trak2_txs


###################################################
### code chunk number 46: trak_4
###################################################
trak2_tx_names <- mcols(trak2_txs)$tx_name
trak2_tx_names


###################################################
### code chunk number 47: trak_5
###################################################
trak2_exbytx <- exonsBy(txdb, "tx", use.names=TRUE)[trak2_tx_names]
elementNROWS(trak2_exbytx)


###################################################
### code chunk number 48: trak_7
###################################################
trak2_inbytx <- intronsByTranscript(txdb, use.names=TRUE)[trak2_tx_names]
elementNROWS(trak2_inbytx)


###################################################
### code chunk number 49: trak_8
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)


###################################################
### code chunk number 50: trak_9
###################################################
trak2_ex_seqs <- getSeq(Hsapiens, trak2_exbytx)
trak2_ex_seqs
trak2_ex_seqs[["uc002uyb.4"]]
trak2_ex_seqs[["uc002uyc.2"]]


###################################################
### code chunk number 51: trak_10
###################################################
trak2_in_seqs <- getSeq(Hsapiens, trak2_inbytx)
trak2_in_seqs
trak2_in_seqs[["uc002uyb.4"]]
trak2_in_seqs[["uc002uyc.2"]]


###################################################
### code chunk number 52: cancer_1
###################################################
library(KEGGREST)
li <- keggList("pathway", "hsa")
ptag <- names(grep("Colorectal cancer", li, value=TRUE))
ptag
tag <- gsub("path:hsa", "", ptag)


###################################################
### code chunk number 53: cancer_2
###################################################
library(KEGGgraph)
dest <- tempfile()
retrieveKGML(tag, "hsa", dest)


###################################################
### code chunk number 54: cancer_3
###################################################
crids <- as.character(parseKGML2DataFrame(dest)[,1])
crgenes <- unique(translateKEGGID2GeneID(crids))
head(crgenes)


###################################################
### code chunk number 55: cancer_4
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


###################################################
### code chunk number 56: cancer_5
###################################################
txbygene <- transcriptsBy(txdb, "gene")[crgenes] ## subset on colorectal genes
map <- relist(unlist(txbygene, use.names=FALSE)$tx_id, txbygene)
map


###################################################
### code chunk number 57: cancer_6
###################################################
cds <- cdsBy(txdb, "tx")
threeUTR <- threeUTRsByTranscript(txdb)
fiveUTR <- fiveUTRsByTranscript(txdb)


###################################################
### code chunk number 58: cancer_7
###################################################
txid <- unlist(map, use.names=FALSE)
cds <- cds[names(cds) %in% txid]
threeUTR <- threeUTR[names(threeUTR) %in% txid]
fiveUTR <- fiveUTR[names(fiveUTR) %in% txid]


###################################################
### code chunk number 59: cancer_8
###################################################
length(txid) ## all possible transcripts
length(cds)
length(threeUTR)
length(fiveUTR)


###################################################
### code chunk number 60: cancer_9
###################################################
cds


###################################################
### code chunk number 61: cancer_10
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19


###################################################
### code chunk number 62: cancer_11
###################################################
threeUTR_seqs <- extractTranscriptSeqs(genome, threeUTR) 
fiveUTR_seqs <- extractTranscriptSeqs(genome, fiveUTR) 
cds_seqs <- extractTranscriptSeqs(genome, cds) 


###################################################
### code chunk number 63: cancer_12
###################################################
cds_seqs


###################################################
### code chunk number 64: cancer_13
###################################################
lst3 <- relist(threeUTR_seqs, PartitioningByWidth(sum(map %in% names(threeUTR))))
lst5 <- relist(fiveUTR_seqs, PartitioningByWidth(sum(map %in% names(fiveUTR))))
lstc <- relist(cds_seqs, PartitioningByWidth(sum(map %in% names(cds))))


###################################################
### code chunk number 65: cancer_14
###################################################
length(map)
table(elementNROWS(map))


###################################################
### code chunk number 66: cancer_15
###################################################
table(elementNROWS(lstc))
table(elementNROWS(lst3))
names(lst3)[elementNROWS(lst3) == 0L] ## genes with no 3' UTR data
table(elementNROWS(lst5))
names(lst5)[elementNROWS(lst5) == 0L] ## genes with no 5' UTR data


###################################################
### code chunk number 67: cseq_1
###################################################
library(Rsamtools)
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
param <- ScanBamParam(what=c("seq", "qual"))
library(GenomicAlignments)
gal <- readGAlignments(bamfile, use.names=TRUE, param=param)


###################################################
### code chunk number 68: cseq_2
###################################################
qseq <- setNames(mcols(gal)$seq, names(gal))
qual <- setNames(mcols(gal)$qual, names(gal))
qseq_on_ref <- sequenceLayer(qseq, cigar(gal),
                             from="query", to="reference")
qual_on_ref <- sequenceLayer(qual, cigar(gal),
                             from="query", to="reference")


###################################################
### code chunk number 69: cseq_3
###################################################
qseq_on_ref_by_chrom <- splitAsList(qseq_on_ref, seqnames(gal))
qual_on_ref_by_chrom <- splitAsList(qual_on_ref, seqnames(gal))
pos_by_chrom <- splitAsList(start(gal), seqnames(gal))


###################################################
### code chunk number 70: cseq_4
###################################################
gr_by_chrom <- lapply(seqlevels(gal),
  function(seqname)
  {
    qseq_on_ref2 <- qseq_on_ref_by_chrom[[seqname]]
    qual_on_ref2 <- qual_on_ref_by_chrom[[seqname]]
    pos2 <- pos_by_chrom[[seqname]]
    qseq_on_ref_per_pos <- split(qseq_on_ref2, pos2)
    qual_on_ref_per_pos <- split(qual_on_ref2, pos2)
    nread <- elementNROWS(qseq_on_ref_per_pos)
    gr_mcols <- DataFrame(nread=unname(nread),
                          qseq_on_ref=unname(qseq_on_ref_per_pos),
                          qual_on_ref=unname(qual_on_ref_per_pos))
    gr <- GRanges(Rle(seqname, nrow(gr_mcols)),
                  IRanges(as.integer(names(nread)), width=1))
    mcols(gr) <- gr_mcols
    seqlevels(gr) <- seqlevels(gal)
    gr
  })


###################################################
### code chunk number 71: cseq_5
###################################################
gr <- do.call(c, gr_by_chrom)
seqinfo(gr) <- seqinfo(gal)


###################################################
### code chunk number 72: cseq_6
###################################################
gr[1:6]


###################################################
### code chunk number 73: cseq_7
###################################################
qseq_on_ref
qual_on_ref


###################################################
### code chunk number 74: cseq_8
###################################################
mcols(gr)$qseq_on_ref[[6]]


###################################################
### code chunk number 75: cseq_9
###################################################
mcols(gr)$qual_on_ref[[6]]


###################################################
### code chunk number 76: cseq_10
###################################################
qseq_on_ref <- mcols(gr)$qseq_on_ref
tmp <- unlist(qseq_on_ref, use.names=FALSE)
qseq_on_ref_id <- relist(match(tmp, tmp), qseq_on_ref)


###################################################
### code chunk number 77: cseq_11
###################################################
qseq_on_ref_id


###################################################
### code chunk number 78: cseq_12
###################################################
qseq_on_ref_id2 <- endoapply(qseq_on_ref_id,
    function(ids) ids[countMatches(ids, ids) >= 0.2 * length(ids)])


###################################################
### code chunk number 79: cseq_13
###################################################
tmp <- unlist(qseq_on_ref_id2, use.names=FALSE)
qseq_on_ref2 <- relist(unlist(qseq_on_ref, use.names=FALSE)[tmp],
                       qseq_on_ref_id2)


###################################################
### code chunk number 80: cseq_14
###################################################
split_factor <- rep.int(seqnames(gr), elementNROWS(qseq_on_ref2))
qseq_on_ref2 <- unlist(qseq_on_ref2, use.names=FALSE)
qseq_on_ref2_by_chrom <- splitAsList(qseq_on_ref2, split_factor)
qseq_pos_by_chrom <- splitAsList(start(gr), split_factor)

cm_by_chrom <- lapply(names(qseq_pos_by_chrom),
    function(seqname)
        consensusMatrix(qseq_on_ref2_by_chrom[[seqname]],
                        as.prob=TRUE,
                        shift=qseq_pos_by_chrom[[seqname]]-1,
                        width=seqlengths(gr)[[seqname]]))
names(cm_by_chrom) <- names(qseq_pos_by_chrom)


###################################################
### code chunk number 81: cseq_15
###################################################
lapply(cm_by_chrom, dim)


###################################################
### code chunk number 82: cseq_16
###################################################
cs_by_chrom <- lapply(cm_by_chrom,
    function(cm) {
        ## need to "fix" 'cm' because consensusString()
        ## doesn't like consensus matrices with columns
        ## that contain only zeroes (e.g., chromosome
        ## positions with no coverage)
        idx <- colSums(cm) == 0L
        cm["+", idx] <- 1
        DNAString(consensusString(cm, ambiguityMap="N"))
    })


###################################################
### code chunk number 83: cseq_17
###################################################
cs_by_chrom


###################################################
### code chunk number 84: bin_1
###################################################
library(BSgenome.Scerevisiae.UCSC.sacCer2)
set.seed(55)
my_var <- RleList(
    lapply(seqlengths(Scerevisiae),
        function(seqlen) {
            tmp <- sample(50L, seqlen, replace=TRUE) %/% 50L
            Rle(cumsum(tmp - rev(tmp)))
        }
    ),
    compress=FALSE)
my_var


###################################################
### code chunk number 85: bin_2
###################################################
bins <- tileGenome(seqinfo(Scerevisiae), tilewidth=100,
                   cut.last.tile.in.chrom=TRUE)


###################################################
### code chunk number 86: bin_3
###################################################
binnedAverage(bins, my_var, "binned_var")


###################################################
### code chunk number 87: SessionInfo
###################################################
sessionInfo()


