library(SCOPE)
library(BSgenome.Hsapiens.UCSC.hg38)
library(WGSmapp)
bamfolder <- "/home/fqin/SRA001/fastq/pair_bam"
bamFile <- as.character(read.table("/home/fqin/SRA001/fastq/pair.bam.files.txt")$V1)
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, 
                         hgref = "hg19", resolution = 100)
bamdir <- bambedObj$bamdir
sampname_raw <- bambedObj$sampname
ref_raw <- bambedObj$ref

mapp <- get_mapp(ref_raw, hgref = "hg19")
head(mapp)
gc <- get_gc(ref_raw, hgref = "hg19")
head(gc)
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
ref_raw

save(bambedObj, file="/home/fqin/SRA001/fastq/pair.bambedObj.RData")
save(ref_raw, file="/home/fqin/SRA001/fastq/pair.ref_raw.RData")

# Getting raw read depth
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, 
                                  seq = 'paired-end', hgref = "hg19")
save(coverageObj, file="/home/fqin/SRA001/fastq/pair.coverageObj.RData")

QCmetric_raw <- get_samp_QC(bambedObj)
save(QCmetric_raw, file="/home/fqin/SRA001/fastq/pair.QCmetric_raw.RData")