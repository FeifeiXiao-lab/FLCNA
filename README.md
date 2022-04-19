# FLCNA
A statistical learning method for simultaneous copy number estimation and subclone clustering with single cell sequencing data.

## Author
Fei Qin, Guoshuai Cai, Feifei Xiao

## Description
Most CNA detection methods for scDNA-seq were designed to detect CNAs and identify subclones in separate ways, which may generate spurious information (e.g., false positive CNAs) during the procedure of CNA detection and thereafter diminish the accuracy of identifying subpopulations from large complex cell group. To overcome this limitation, we developed a fused lasso mode-based framework, FLCNA, for CNA detection and simultaneous subclone identification in scDNA-seq data. 

## Installation
```r
install.packages("devtools")
library(devtools)
install_github("FeifeiXiaoUSC/FLCNA")
```

## Running FLCNA
### Examples

```r
library(FLCNV)
```
