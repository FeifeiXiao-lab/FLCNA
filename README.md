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
install_github("FeifeiXiaoUSC/FLCNA")
```

## Running FLCNA
### Examples

```r
# The example data have 2,000 markers and 200 cells.
library(FLCNA)
data(Example_data_2000)
data(Example_ref_2000)
RD <- Example_data_2000
ref <- Example_ref_2000
dim(RD)
[1]  200 2000
```


```r
ref <- Example_ref_2000
head(ref)

GRanges object with 6 ranges and 2 metadata columns:
      seqnames          ranges strand |        gc      mapp
         <Rle>       <IRanges>  <Rle> | <numeric> <numeric>
  [1]     chr1 2000001-2100000      * |     56.96  0.984862
  [2]     chr1 2800001-2900000      * |     57.94  0.992544
  [3]     chr1 2900001-3000000      * |     55.43  0.984850
  [4]     chr1 3000001-3100000      * |     56.60  0.995182
  [5]     chr1 3100001-3200000      * |     58.16  0.989534
  [6]     chr1 3200001-3300000      * |     56.83  0.973831
  -------
  seqinfo: 24 sequences from hg38 genome
```

```r
# Quality Control 
QCobject <- FLCNA_QC(Y_raw=t(RD), 
                     ref_raw=ref,
                     mapp_thresh = 0.9,
                     gc_thresh = c(20, 80))
```

```r
# Normalization
log2Rdata <- FLCNA_normalization(Y=QCobject$Y, gc=QCobject$ref$gc, map=QCobject$ref$mapp)
```
```r
# Simultaneous CNA detection and subclone clustering
output_FLCNA <- FLCNA(K=c(4,5,6), lambda=3, Y=data.matrix(log2Rdata))
```
```r
# CNA clustering
CNA.output <- CNA.out(mean.matrix = output_FLCNA$mu.hat.best, LRR=log2Rdata, Clusters=output_FLCNA$s.hat.best, ref=ref, cutoff=0.35, L=100)
```
