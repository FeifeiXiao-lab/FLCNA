library(HMMcopy)
setwd("/home/fqin/Fused.lasso/simulation_3clusters_5states/Simulation100Percent")
Y_sim <- t(get(load("RD_Ddup_medium_100.RData")))
qcObj <- get(load("qcObj_QC_clear.RData"))
ref <- qcObj$ref

copy.output <- NULL
for (i in 1:ncol(Y_sim)){
  normal_reads1 <- data.frame(chr=rep(ref@seqnames@values,ref@seqnames@lengths),
                              start=ref@ranges@start,
                              end=ref@ranges@start+ref@ranges@width,
                              reads=Y_sim[,i],
                              gc=ref$gc,
                              map=ref$mapp)
  normal_reads1$chr <- substr(normal_reads1$chr,4,nchar(as.character(normal_reads1$chr)))
  normal_reads1$gc <- normal_reads1$gc/100
  normal_copy <- correctReadcount(normal_reads1)

  default_param <- HMMsegment(normal_copy, getparam = TRUE)
  longseg_param <- default_param
  longseg_param$e <- 0.9999
  longseg_param$strength <- 1e30
  longseg_segments <- HMMsegment(normal_copy, longseg_param, verbose = FALSE)
  table(longseg_segments$state)
  copy.output <- cbind(copy.output, longseg_segments$state)
}
colnames(copy.output) <- colnames(Y_sim)

clusters <- hclust(dist(t(copy.output)))
clusterCut <- cutree(clusters, 3)

output <- list(copy.output=copy.output, clusters=clusters, clusterCut=clusterCut)

setwd("/home/fqin/Fused.lasso/simulation_3clusters_5states/Simulation100Percent/HMMcopy")
save(output, file="HMMcopy_Ddup_medium100_output.RData")
