# FLCNV
An accurate CNV detection and subclone identfication method for single cell DNA sequencing (scDNA-seq) Data

## Author
Fei Qin, Feifei Xiao

## Description
Most CNV detection methods for scDNA-seq were designed to detect CNVs and identify subclones in separate ways, which may generate spurious information (e.g., false positive CNVs) during the procedure of CNV detection and thereafter diminish the accuracy of identifying subpopulations from large cell group. To overcome this limitation, we developed a fused lasso mode-based framework, FLCNV, for CNV detection and simultaneous subclone identification in scDNA-seq data. 

## Installation
```r
install.packages("devtools")
library(devtools)
install_github("FeiQin92/FLCNV")
```

## Running FLCNV
### Examples

```r
library(FLCNV)
log2Rdata <- rbind(matrix(rnorm(10000, 0, 0.5), 10, 1000), matrix(rnorm(10000, 0.5, 0.5), 10, 1000)
output <- FLCNV(K=2, lambda=c(5,10), y=log2Rdata)
expre_data = exprs(EMTAB.eset)
pheno_data = pData(EMTAB.eset)
```
