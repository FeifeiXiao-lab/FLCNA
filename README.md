# FLCNA
A statistical learning method for simultaneous copy number estimation and subclone clustering with single cell sequencing data.

## Author
Fei Qin, Guoshuai Cai, Feifei Xiao

## Description
Most CNA detection methods for scDNA-seq were designed to detect CNAs and identify subclones in separate ways, which may generate spurious information (e.g., false positive CNAs) during the procedure of CNA detection and thereafter diminish the accuracy of identifying subpopulations from large complex cell group. To overcome this limitation, we developed a fused lasso mode-based framework, FLCNA, for CNA detection and simultaneous subclone identification in scDNA-seq data. First, procedures including quality control (QC), normalization, logarithm transformation were used for pre-processing of the datasets. Subclone clustering was achieved based on a Gaussian Mixture Model (GMM), and breakpoints detection was conducted by adding a fused lasso penalty term to the typical GMM model. Finally, based on these shared breakpoints in each cluster, candidate CNA segments for each cell were clustered into different CNA states using a GMM-based clustering strategy. 

## Installation
```r
install.packages("devtools")
library(devtools)
install_github("FeifeiXiao-lab/FLCNA")
```

## Running FLCNA
### Examples

```r
# The example data have 3,000 markers and 93 cells.
library(FLCNA)
data(KTN126_data_3000)
data(KTN126_ref_3000)
RD <- KTN126_data_3000
dim(RD)
[1]  3000   93
```


```r
ref <- KTN126_ref_3000
head(ref)

GRanges object with 6 ranges and 2 metadata columns:
      seqnames          ranges strand |        gc      mapp
         <Rle>       <IRanges>  <Rle> | <numeric> <numeric>
  [1]     chr1   800001-900000      * |     56.65  0.920803
  [2]     chr1  900001-1000000      * |     62.13  0.973863
  [3]     chr1 1000001-1100000      * |     60.35  0.946335
  [4]     chr1 1100001-1200000      * |     61.37  0.976173
  [5]     chr1 1200001-1300000      * |     63.52  0.969312
  [6]     chr1 1300001-1400000      * |     56.19  0.946105
  -------
  seqinfo: 24 sequences from hg19 genome
```

```r
# Quality Control 
QCobject <- FLCNA_QC(Y_raw=RD, ref_raw=ref,
                     mapp_thresh = 0.9,
                     gc_thresh = c(20, 80))
```

```r
# Normalization
log2Rdata <- FLCNA_normalization(Y=QCobject$Y, gc=QCobject$ref$gc, map=QCobject$ref$mapp)
```
```r
# Simultaneous CNA detection and subclone clustering
output_FLCNA <- FLCNA(K=c(4,5,6), lambda=3, Y=t(log2Rdata), ref=QCobject$ref)
```

```r
# CNA clustering
CNA.output <- CNA.out(mean.matrix = res$mu.hat.best, Clusters=output$s.hat.best,
                      LRR=log2Rdata, QC_ref=QCobject$ref, cutoff=0.80, L=100)
```
