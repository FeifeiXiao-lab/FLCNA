# Generate simulation data using spike-in strategy for 3 clusters
########################################################################

qcObj <- get(load("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Gen_simu_data\\qcObj_QC_clear.RData"))
Y <- qcObj$Y
rd2 <- Y[,21:220]
chr.index <- cumsum(qcObj$ref@seqnames@lengths)
set.seed(2021)
# factor <- matrix(1,nrow(Y_raw),ncol(Y_raw))
setwd("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states")
set.seed(2021)
random.Cell <- sample(1:3,ncol(rd2),replace = TRUE)
save(random.Cell, file="random.Cell.3clusters.RData")

nCluster=3
nCNV=50
Percent=100
Gen_simu <- function(rd2, nCluster, nCNV, Percent){
  ##Super short
  RD_Ddel_SuperShort=t(rd2)
  RD_Sdel_SuperShort=t(rd2)
  
  RD_Sdup_SuperShort=t(rd2)
  RD_Ddup_SuperShort=t(rd2)
  
  RD_mix_SuperShort=t(rd2)
  
  cpSuperShort.list <- NULL
  locations.SuperShort=NULL
  mix.CNVindex.list <- NULL
  for (j in 1:nCluster){
    
    ##Super short
    min=2
    max=5
    min1=200
    max1=480
    set.seed(j)
    cnv_len=as.integer(runif(n=nCNV,min = min,max = max))
    normal_len=as.integer(runif(n=nCNV+1,min = min1,max = max1))
    start=vector()
    end=vector()
    start[1]=normal_len[1]
    end[1]=normal_len[1]+cnv_len[1]
    for(i in 2:nCNV){
      start[i]=end[i-1]+normal_len[i]
      end[i]=start[i]+cnv_len[i]
      for (s in 1:22){
        s.chr <- chr.index[s]
        if (start[i] <= s.chr & end[i] > s.chr){
          start[i]=end[i]+normal_len[i]
          end[i]=start[i]+cnv_len[i]
        }
      }
    }
    cpSuperShort=cbind(start,end)
    cpSuperShort
    cpSuperShort.list <- cbind(cpSuperShort.list, cpSuperShort)
    
    
    
    random_Ddel_SuperShort=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      
      locate=data.frame(start=rep(start[i], length(random.index.CNVpercent)),
                        end=rep(end[i], length(random.index.CNVpercent)),
                        sample=random.index.CNVpercent,
                        cluster=j)
      locations.SuperShort=rbind(locations.SuperShort, locate)
      
      n1=dim(RD_Ddel_SuperShort[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddel_SuperShort[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddel=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2 
      random_Ddel_SuperShort[[i]]=random_Ddel
      RD_Ddel_SuperShort[random.index.CNVpercent, start[i]:end[i]]=RD_Ddel_SuperShort[random.index.CNVpercent, start[i]:end[i]]*random_Ddel
    }
    
    random_Sdel_SuperShort=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdel_SuperShort[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdel_SuperShort[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdel=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2 
      random_Sdel_SuperShort[[i]]=random_Sdel
      RD_Sdel_SuperShort[random.index.CNVpercent, start[i]:end[i]]=RD_Sdel_SuperShort[random.index.CNVpercent, start[i]:end[i]]*random_Sdel
    }
    
    
    random_Sdup_SuperShort=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdup_SuperShort[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdup_SuperShort[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdup=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2
      random_Sdup_SuperShort[[i]]=random_Sdup
      RD_Sdup_SuperShort[random.index.CNVpercent, start[i]:end[i]]=RD_Sdup_SuperShort[random.index.CNVpercent, start[i]:end[i]]*random_Sdup
    }
    
    random_Ddup_SuperShort=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Ddup_SuperShort[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddup_SuperShort[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddup=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2
      random_Ddup_SuperShort[[i]]=random_Ddup
      RD_Ddup_SuperShort[random.index.CNVpercent, start[i]:end[i]]=RD_Ddup_SuperShort[random.index.CNVpercent, start[i]:end[i]]*random_Ddup
    }
    
    random_mix_SuperShort=list()
    set.seed(j)
    CNV.index <- sample(c(0, 1, 3, 4), nCNV, prob=c(0.1, 0.5, 0.3, 0.1), replace=TRUE)
    mix.CNVindex.list <- rbind(mix.CNVindex.list, CNV.index)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_mix_SuperShort[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_mix_SuperShort[random.index.CNVpercent, start[i]:end[i]])[2]
      
      if (CNV.index[i]==0) {random_mix=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==1) {random_mix=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==3) {random_mix=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==4) {random_mix=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2}
      
      random_mix_SuperShort[[i]]=random_mix
      
      RD_mix_SuperShort[random.index.CNVpercent, start[i]:end[i]]=RD_mix_SuperShort[random.index.CNVpercent, start[i]:end[i]]*random_mix
    }
  }
  
  
  
  ## short
  RD_Ddel_short=t(rd2)
  RD_Sdel_short=t(rd2)
  
  RD_Sdup_short=t(rd2)
  RD_Ddup_short=t(rd2)
  
  RD_mix_short=t(rd2)
  
  cpshort.list <- NULL
  locations.short=NULL
  mix.CNVindex.list <- NULL
  for (j in 1:nCluster){
    ##short
    min=5
    max=10
    min1=200
    max1=500
    set.seed(j)
    cnv_len=as.integer(runif(n=nCNV,min = min,max = max))
    normal_len=as.integer(runif(n=nCNV+1,min = min1,max = max1))
    start=vector()
    end=vector()
    start[1]=normal_len[1]
    end[1]=normal_len[1]+cnv_len[1]
    for(i in 2:nCNV){
      start[i]=end[i-1]+normal_len[i]
      end[i]=start[i]+cnv_len[i]
      for (s in 1:22){
        s.chr <- chr.index[s]
        if (start[i] <= s.chr & end[i] > s.chr){
          start[i]=end[i]+normal_len[i]
          end[i]=start[i]+cnv_len[i]
        }
      }
    }
    cpshort=cbind(start,end)
    cpshort
    cpshort.list <- cbind(cpshort.list, cpshort)
    
    random_Ddel_short=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      
      locate=data.frame(start=rep(start[i], length(random.index.CNVpercent)),
                        end=rep(end[i], length(random.index.CNVpercent)),
                        sample=random.index.CNVpercent,
                        cluster=j)
      locations.short=rbind(locations.short, locate)
      
      n1=dim(RD_Ddel_short[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddel_short[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddel=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2 
      random_Ddel_short[[i]]=random_Ddel
      RD_Ddel_short[random.index.CNVpercent, start[i]:end[i]]=RD_Ddel_short[random.index.CNVpercent, start[i]:end[i]]*random_Ddel
    }
    
    random_Sdel_short=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdel_short[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdel_short[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdel=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2 
      random_Sdel_short[[i]]=random_Sdel
      RD_Sdel_short[random.index.CNVpercent, start[i]:end[i]]=RD_Sdel_short[random.index.CNVpercent, start[i]:end[i]]*random_Sdel
    }
    
    
    random_Sdup_short=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdup_short[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdup_short[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdup=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2
      random_Sdup_short[[i]]=random_Sdup
      RD_Sdup_short[random.index.CNVpercent, start[i]:end[i]]=RD_Sdup_short[random.index.CNVpercent, start[i]:end[i]]*random_Sdup
    }
    
    random_Ddup_short=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Ddup_short[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddup_short[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddup=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2
      random_Ddup_short[[i]]=random_Ddup
      RD_Ddup_short[random.index.CNVpercent, start[i]:end[i]]=RD_Ddup_short[random.index.CNVpercent, start[i]:end[i]]*random_Ddup
    }
    
    random_mix_short=list()
    set.seed(j)
    CNV.index <- sample(c(0, 1, 3, 4), nCNV, prob=c(0.1, 0.5, 0.3, 0.1), replace=TRUE)
    mix.CNVindex.list <- rbind(mix.CNVindex.list, CNV.index)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_mix_short[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_mix_short[random.index.CNVpercent, start[i]:end[i]])[2]
      
      if (CNV.index[i]==0) {random_mix=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==1) {random_mix=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==3) {random_mix=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==4) {random_mix=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2}
      
      random_mix_short[[i]]=random_mix
      
      RD_mix_short[random.index.CNVpercent, start[i]:end[i]]=RD_mix_short[random.index.CNVpercent, start[i]:end[i]]*random_mix
    }
  }
  
  
  
  ##medium
  RD_Ddel_medium=t(rd2)
  RD_Sdel_medium=t(rd2)
  
  RD_Sdup_medium=t(rd2)
  RD_Ddup_medium=t(rd2)
  RD_mix_medium=t(rd2)
  
  cpmedium.list <- NULL
  locations.medium=NULL
  mix.CNVindex.list <- NULL
  for (j in 1:nCluster){
    
    min=10
    max=20
    min1=200
    max1=460
    set.seed(j)
    cnv_len=as.integer(runif(n=nCNV,min = min,max = max))
    normal_len=as.integer(runif(n=nCNV+1,min = min1,max = max1))
    start=vector()
    end=vector()
    start[1]=normal_len[1]
    end[1]=normal_len[1]+cnv_len[1]
    for(i in 2:nCNV){
      start[i]=end[i-1]+normal_len[i]
      end[i]=start[i]+cnv_len[i]
      for (s in 1:22){
        s.chr <- chr.index[s]
        if (start[i] <= s.chr & end[i] > s.chr){
          start[i]=end[i]+normal_len[i]
          end[i]=start[i]+cnv_len[i]
        }
      }
    }
    cpmedium=cbind(start,end)
    cpmedium
    cpmedium.list <- cbind(cpmedium.list, cpmedium)
    
    random_Ddel_medium=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      
      locate=data.frame(start=rep(start[i], length(random.index.CNVpercent)),
                        end=rep(end[i], length(random.index.CNVpercent)),
                        sample=random.index.CNVpercent,
                        cluster=j)
      locations.medium=rbind(locations.medium, locate)
      
      n1=dim(RD_Ddel_medium[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddel_medium[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddel=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2 
      random_Ddel_medium[[i]]=random_Ddel
      RD_Ddel_medium[random.index.CNVpercent, start[i]:end[i]]=RD_Ddel_medium[random.index.CNVpercent, start[i]:end[i]]*random_Ddel
    }
    
    random_Sdel_medium=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdel_medium[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdel_medium[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdel=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2 
      random_Sdel_medium[[i]]=random_Sdel
      RD_Sdel_medium[random.index.CNVpercent, start[i]:end[i]]=RD_Sdel_medium[random.index.CNVpercent, start[i]:end[i]]*random_Sdel
    }
    
    
    random_Sdup_medium=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdup_medium[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdup_medium[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdup=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2
      random_Sdup_medium[[i]]=random_Sdup
      RD_Sdup_medium[random.index.CNVpercent, start[i]:end[i]]=RD_Sdup_medium[random.index.CNVpercent, start[i]:end[i]]*random_Sdup
    }
    
    random_Ddup_medium=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Ddup_medium[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddup_medium[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddup=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2
      random_Ddup_medium[[i]]=random_Ddup
      RD_Ddup_medium[random.index.CNVpercent, start[i]:end[i]]=RD_Ddup_medium[random.index.CNVpercent, start[i]:end[i]]*random_Ddup
    }
    
    random_mix_medium=list()
    set.seed(j)
    CNV.index <- sample(c(0, 1, 3, 4), nCNV, prob=c(0.1, 0.5, 0.3, 0.1), replace=TRUE)
    mix.CNVindex.list <- rbind(mix.CNVindex.list, CNV.index)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_mix_medium[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_mix_medium[random.index.CNVpercent, start[i]:end[i]])[2]
      
      if (CNV.index[i]==0) {random_mix=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==1) {random_mix=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==3) {random_mix=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==4) {random_mix=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2}
      
      random_mix_medium[[i]]=random_mix
      
      RD_mix_medium[random.index.CNVpercent, start[i]:end[i]]=RD_mix_medium[random.index.CNVpercent, start[i]:end[i]]*random_mix
    }
  }
  
  
  
  
  ##long
  RD_Ddel_long=t(rd2)
  RD_Sdel_long=t(rd2)
  
  RD_Sdup_long=t(rd2)
  RD_Ddup_long=t(rd2)
  RD_mix_long=t(rd2)
  
  cplong.list <- NULL
  locations.long=NULL
  mix.CNVindex.list <- NULL
  for (j in 1:nCluster){
    
    min=20
    max=35
    min1=200
    max1=450
    set.seed(j)
    cnv_len=as.integer(runif(n=nCNV,min = min,max = max))
    normal_len=as.integer(runif(n=nCNV+1,min = min1,max = max1))
    start=vector()
    end=vector()
    start[1]=normal_len[1]
    end[1]=normal_len[1]+cnv_len[1]
    for(i in 2:nCNV){
      start[i]=end[i-1]+normal_len[i]
      end[i]=start[i]+cnv_len[i]
      for (s in 1:22){
        s.chr <- chr.index[s]
        if (start[i] <= s.chr & end[i] > s.chr){
          start[i]=end[i]+normal_len[i]
          end[i]=start[i]+cnv_len[i]
        }
      }
    }
    end
    cplong=cbind(start,end)
    cplong
    cplong.list <- cbind(cplong.list, cplong)
    
    random_Ddel_long=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      
      locate=data.frame(start=rep(start[i], length(random.index.CNVpercent)),
                        end=rep(end[i], length(random.index.CNVpercent)),
                        sample=random.index.CNVpercent,
                        cluster=j)
      locations.long=rbind(locations.long, locate)
      
      n1=dim(RD_Ddel_long[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddel_long[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddel=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2 
      random_Ddel_long[[i]]=random_Ddel
      RD_Ddel_long[random.index.CNVpercent, start[i]:end[i]]=RD_Ddel_long[random.index.CNVpercent, start[i]:end[i]]*random_Ddel
    }
    
    random_Sdel_long=list()
    for(i in 1:nCNV){  
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdel_long[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdel_long[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdel=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2 
      random_Sdel_long[[i]]=random_Sdel
      RD_Sdel_long[random.index.CNVpercent, start[i]:end[i]]=RD_Sdel_long[random.index.CNVpercent, start[i]:end[i]]*random_Sdel
    }
    
    
    random_Sdup_long=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Sdup_long[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Sdup_long[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Sdup=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2
      random_Sdup_long[[i]]=random_Sdup
      RD_Sdup_long[random.index.CNVpercent, start[i]:end[i]]=RD_Sdup_long[random.index.CNVpercent, start[i]:end[i]]*random_Sdup
    }
    
    random_Ddup_long=list()
    set.seed(j)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_Ddup_long[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_Ddup_long[random.index.CNVpercent, start[i]:end[i]])[2]
      random_Ddup=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2
      random_Ddup_long[[i]]=random_Ddup
      RD_Ddup_long[random.index.CNVpercent, start[i]:end[i]]=RD_Ddup_long[random.index.CNVpercent, start[i]:end[i]]*random_Ddup
    }
    
    random_mix_long=list()
    set.seed(j)
    CNV.index <- sample(c(0, 1, 3, 4), nCNV, prob=c(0.1, 0.5, 0.3, 0.1), replace=TRUE)
    mix.CNVindex.list <- rbind(mix.CNVindex.list, CNV.index)
    for(i in 1:nCNV){
      set.seed((j-1)*nCNV+i)
      random.index.CNVpercent <- sample(which(random.Cell==j),round(sum(random.Cell==j)*Percent/100))
      random.index.CNVpercent <-  random.index.CNVpercent[order(random.index.CNVpercent)]
      n1=dim(RD_mix_long[random.index.CNVpercent, start[i]:end[i]])[1]
      n2=dim(RD_mix_long[random.index.CNVpercent, start[i]:end[i]])[2]
      
      if (CNV.index[i]==0) {random_mix=matrix(rnorm(n1*n2, mean=0.4, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==1) {random_mix=matrix(rnorm(n1*n2, mean=1.2, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==3) {random_mix=matrix(rnorm(n1*n2, mean=2.8, sd=0.1), n1, n2)/2}
      if (CNV.index[i]==4) {random_mix=matrix(rnorm(n1*n2, mean=4.2, sd=0.1), n1, n2)/2}
      
      random_mix_long[[i]]=random_mix
      
      RD_mix_long[random.index.CNVpercent, start[i]:end[i]]=RD_mix_long[random.index.CNVpercent, start[i]:end[i]]*random_mix
    }
  }
  
  
  res <- list(RD_Sdel_SuperShort=RD_Sdel_SuperShort, RD_Ddel_SuperShort=RD_Ddel_SuperShort,
              RD_Sdup_SuperShort=RD_Sdup_SuperShort, RD_Ddup_SuperShort=RD_Ddup_SuperShort,
              RD_mix_SuperShort=RD_mix_SuperShort, cpSuperShort.list=cpSuperShort.list,
              
              RD_Sdel_short=RD_Sdel_short, RD_Ddel_short=RD_Ddel_short,
              RD_Sdup_short=RD_Sdup_short, RD_Ddup_short=RD_Ddup_short,
              RD_mix_short=RD_mix_short, cpshort.list=cpshort.list,
              
              RD_Sdel_medium=RD_Sdel_medium, RD_Ddel_medium=RD_Ddel_medium,
              RD_Sdup_medium=RD_Sdup_medium, RD_Ddup_medium=RD_Ddup_medium,
              RD_mix_medium=RD_mix_medium, cpmedium.list=cpmedium.list,
              
              RD_Sdel_long=RD_Sdel_long, RD_Ddel_long=RD_Ddel_long,
              RD_Sdup_long=RD_Sdup_long, RD_Ddup_long=RD_Ddup_long,
              RD_mix_long=RD_mix_long, cplong.list=cplong.list,
              
              locations.SuperShort=locations.SuperShort, locations.short=locations.short, 
              locations.medium=locations.medium, locations.long=locations.long,
              
              mix.CNVindex.list = mix.CNVindex.list)
  return(res)
}






res <-Gen_simu(rd2=rd2, nCluster=3, nCNV=50, Percent=100)

setwd("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation100Percent")
mix.CNVindex.list <- res$mix.CNVindex.list
save(mix.CNVindex.list, file="mix.CNVindex.list.RData")

RD_Sdel_SuperShort <- res$RD_Sdel_SuperShort
RD_Ddel_SuperShort <- res$RD_Ddel_SuperShort
RD_Sdup_SuperShort <- res$RD_Sdup_SuperShort
RD_Ddup_SuperShort <- res$RD_Ddup_SuperShort
RD_mix_SuperShort <- res$RD_mix_SuperShort
cpSuperShort.list <- res$cpSuperShort.list
locations.SuperShort <- res$locations.SuperShort
save(RD_Sdel_SuperShort,file="RD_Sdel_SuperShort_100.RData")
save(RD_Ddel_SuperShort,file="RD_Ddel_SuperShort_100.RData")
save(RD_Sdup_SuperShort,file="RD_Sdup_SuperShort_100.RData")
save(RD_Ddup_SuperShort,file="RD_Ddup_SuperShort_100.RData")
save(RD_mix_SuperShort,file="RD_mix_SuperShort_100.RData")
save(cpSuperShort.list, file="cpSuperShort_100.RData")
save(locations.SuperShort, file="locations.SuperShort100.RData")

RD_Sdel_short <- res$RD_Sdel_short
RD_Ddel_short <- res$RD_Ddel_short
RD_Sdup_short <- res$RD_Sdup_short
RD_Ddup_short <- res$RD_Ddup_short
RD_mix_short <- res$RD_mix_short
cpshort.list <- res$cpshort.list
locations.short <- res$locations.short
save(RD_Sdel_short,file="RD_Sdel_short_100.RData")
save(RD_Ddel_short,file="RD_Ddel_short_100.RData")
save(RD_Sdup_short,file="RD_Sdup_short_100.RData")
save(RD_Ddup_short,file="RD_Ddup_short_100.RData")
save(RD_mix_short,file="RD_mix_short_100.RData")
save(cpshort.list, file="cpshort_100.RData")
save(locations.short, file="locations.short100.RData")

RD_Sdel_medium <- res$RD_Sdel_medium
RD_Ddel_medium <- res$RD_Ddel_medium
RD_Sdup_medium <- res$RD_Sdup_medium
RD_Ddup_medium <- res$RD_Ddup_medium
RD_mix_medium <- res$RD_mix_medium
cpmedium.list <- res$cpmedium.list
locations.medium <- res$locations.medium
save(RD_Sdel_medium,file="RD_Sdel_medium_100.RData")
save(RD_Ddel_medium,file="RD_Ddel_medium_100.RData")
save(RD_Sdup_medium,file="RD_Sdup_medium_100.RData")
save(RD_Ddup_medium,file="RD_Ddup_medium_100.RData")
save(RD_mix_medium,file="RD_mix_medium_100.RData")
save(cpmedium.list, file="cpmedium_100.RData")
save(locations.medium, file="locations.medium100.RData")

RD_Sdel_long <- res$RD_Sdel_long
RD_Ddel_long <- res$RD_Ddel_long
RD_Sdup_long <- res$RD_Sdup_long
RD_Ddup_long <- res$RD_Ddup_long
RD_mix_long <- res$RD_mix_long
cplong.list <- res$cplong.list
locations.long <- res$locations.long
save(RD_Sdel_long,file="RD_Sdel_long_100.RData")
save(RD_Ddel_long,file="RD_Ddel_long_100.RData")
save(RD_Sdup_long,file="RD_Sdup_long_100.RData")
save(RD_Ddup_long,file="RD_Ddup_long_100.RData")
save(RD_mix_long,file="RD_mix_long_100.RData")
save(cplong.list, file="cplong_100.RData")
save(locations.long, file="locations.long100.RData")





res <-Gen_simu(rd2=rd2, nCluster=3, nCNV=50, Percent=80)

setwd("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation80Percent")
mix.CNVindex.list <- res$mix.CNVindex.list
save(mix.CNVindex.list, file="mix.CNVindex.list.RData")

RD_Sdel_SuperShort <- res$RD_Sdel_SuperShort
RD_Ddel_SuperShort <- res$RD_Ddel_SuperShort
RD_Sdup_SuperShort <- res$RD_Sdup_SuperShort
RD_Ddup_SuperShort <- res$RD_Ddup_SuperShort
RD_mix_SuperShort <- res$RD_mix_SuperShort
cpSuperShort.list <- res$cpSuperShort.list
locations.SuperShort <- res$locations.SuperShort
save(RD_Sdel_SuperShort,file="RD_Sdel_SuperShort_80.RData")
save(RD_Ddel_SuperShort,file="RD_Ddel_SuperShort_80.RData")
save(RD_Sdup_SuperShort,file="RD_Sdup_SuperShort_80.RData")
save(RD_Ddup_SuperShort,file="RD_Ddup_SuperShort_80.RData")
save(RD_mix_SuperShort,file="RD_mix_SuperShort_80.RData")
save(cpSuperShort.list, file="cpSuperShort_80.RData")
save(locations.SuperShort, file="locations.SuperShort80.RData")

RD_Sdel_short <- res$RD_Sdel_short
RD_Ddel_short <- res$RD_Ddel_short
RD_Sdup_short <- res$RD_Sdup_short
RD_Ddup_short <- res$RD_Ddup_short
RD_mix_short <- res$RD_mix_short
cpshort.list <- res$cpshort.list
locations.short <- res$locations.short
save(RD_Sdel_short,file="RD_Sdel_short_80.RData")
save(RD_Ddel_short,file="RD_Ddel_short_80.RData")
save(RD_Sdup_short,file="RD_Sdup_short_80.RData")
save(RD_Ddup_short,file="RD_Ddup_short_80.RData")
save(RD_mix_short,file="RD_mix_short_80.RData")
save(cpshort.list, file="cpshort_80.RData")
save(locations.short, file="locations.short80.RData")

RD_Sdel_medium <- res$RD_Sdel_medium
RD_Ddel_medium <- res$RD_Ddel_medium
RD_Sdup_medium <- res$RD_Sdup_medium
RD_Ddup_medium <- res$RD_Ddup_medium
RD_mix_medium <- res$RD_mix_medium
cpmedium.list <- res$cpmedium.list
locations.medium <- res$locations.medium
save(RD_Sdel_medium,file="RD_Sdel_medium_80.RData")
save(RD_Ddel_medium,file="RD_Ddel_medium_80.RData")
save(RD_Sdup_medium,file="RD_Sdup_medium_80.RData")
save(RD_Ddup_medium,file="RD_Ddup_medium_80.RData")
save(RD_mix_medium,file="RD_mix_medium_80.RData")
save(cpmedium.list, file="cpmedium_80.RData")
save(locations.medium, file="locations.medium80.RData")

RD_Sdel_long <- res$RD_Sdel_long
RD_Ddel_long <- res$RD_Ddel_long
RD_Sdup_long <- res$RD_Sdup_long
RD_Ddup_long <- res$RD_Ddup_long
RD_mix_long <- res$RD_mix_long
cplong.list <- res$cplong.list
locations.long <- res$locations.long
save(RD_Sdel_long,file="RD_Sdel_long_80.RData")
save(RD_Ddel_long,file="RD_Ddel_long_80.RData")
save(RD_Sdup_long,file="RD_Sdup_long_80.RData")
save(RD_Ddup_long,file="RD_Ddup_long_80.RData")
save(RD_mix_long,file="RD_mix_long_80.RData")
save(cplong.list, file="cplong_80.RData")
save(locations.long, file="locations.long80.RData")




res <-Gen_simu(rd2=rd2, nCluster=3, nCNV=50, Percent=60)

setwd("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation60Percent")
mix.CNVindex.list <- res$mix.CNVindex.list
save(mix.CNVindex.list, file="mix.CNVindex.list.RData")

RD_Sdel_SuperShort <- res$RD_Sdel_SuperShort
RD_Ddel_SuperShort <- res$RD_Ddel_SuperShort
RD_Sdup_SuperShort <- res$RD_Sdup_SuperShort
RD_Ddup_SuperShort <- res$RD_Ddup_SuperShort
RD_mix_SuperShort <- res$RD_mix_SuperShort
cpSuperShort.list <- res$cpSuperShort.list
locations.SuperShort <- res$locations.SuperShort
save(RD_Sdel_SuperShort,file="RD_Sdel_SuperShort_60.RData")
save(RD_Ddel_SuperShort,file="RD_Ddel_SuperShort_60.RData")
save(RD_Sdup_SuperShort,file="RD_Sdup_SuperShort_60.RData")
save(RD_Ddup_SuperShort,file="RD_Ddup_SuperShort_60.RData")
save(RD_mix_SuperShort,file="RD_mix_SuperShort_60.RData")
save(cpSuperShort.list, file="cpSuperShort_60.RData")
save(locations.SuperShort, file="locations.SuperShort60.RData")

RD_Sdel_short <- res$RD_Sdel_short
RD_Ddel_short <- res$RD_Ddel_short
RD_Sdup_short <- res$RD_Sdup_short
RD_Ddup_short <- res$RD_Ddup_short
RD_mix_short <- res$RD_mix_short
cpshort.list <- res$cpshort.list
locations.short <- res$locations.short
save(RD_Sdel_short,file="RD_Sdel_short_60.RData")
save(RD_Ddel_short,file="RD_Ddel_short_60.RData")
save(RD_Sdup_short,file="RD_Sdup_short_60.RData")
save(RD_Ddup_short,file="RD_Ddup_short_60.RData")
save(RD_mix_short,file="RD_mix_short_60.RData")
save(cpshort.list, file="cpshort_60.RData")
save(locations.short, file="locations.short60.RData")

RD_Sdel_medium <- res$RD_Sdel_medium
RD_Ddel_medium <- res$RD_Ddel_medium
RD_Sdup_medium <- res$RD_Sdup_medium
RD_Ddup_medium <- res$RD_Ddup_medium
RD_mix_medium <- res$RD_mix_medium
cpmedium.list <- res$cpmedium.list
locations.medium <- res$locations.medium
save(RD_Sdel_medium,file="RD_Sdel_medium_60.RData")
save(RD_Ddel_medium,file="RD_Ddel_medium_60.RData")
save(RD_Sdup_medium,file="RD_Sdup_medium_60.RData")
save(RD_Ddup_medium,file="RD_Ddup_medium_60.RData")
save(RD_mix_medium,file="RD_mix_medium_60.RData")
save(cpmedium.list, file="cpmedium_60.RData")
save(locations.medium, file="locations.medium60.RData")

RD_Sdel_long <- res$RD_Sdel_long
RD_Ddel_long <- res$RD_Ddel_long
RD_Sdup_long <- res$RD_Sdup_long
RD_Ddup_long <- res$RD_Ddup_long
RD_mix_long <- res$RD_mix_long
cplong.list <- res$cplong.list
locations.long <- res$locations.long
save(RD_Sdel_long,file="RD_Sdel_long_60.RData")
save(RD_Ddel_long,file="RD_Ddel_long_60.RData")
save(RD_Sdup_long,file="RD_Sdup_long_60.RData")
save(RD_Ddup_long,file="RD_Ddup_long_60.RData")
save(RD_mix_long,file="RD_mix_long_60.RData")
save(cplong.list, file="cplong_60.RData")
save(locations.long, file="locations.long60.RData")




res <-Gen_simu(rd2=rd2, nCluster=3, nCNV=50, Percent=40)

setwd("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation40Percent")
mix.CNVindex.list <- res$mix.CNVindex.list
save(mix.CNVindex.list, file="mix.CNVindex.list.RData")

RD_Sdel_SuperShort <- res$RD_Sdel_SuperShort
RD_Ddel_SuperShort <- res$RD_Ddel_SuperShort
RD_Sdup_SuperShort <- res$RD_Sdup_SuperShort
RD_Ddup_SuperShort <- res$RD_Ddup_SuperShort
RD_mix_SuperShort <- res$RD_mix_SuperShort
cpSuperShort.list <- res$cpSuperShort.list
locations.SuperShort <- res$locations.SuperShort
save(RD_Sdel_SuperShort,file="RD_Sdel_SuperShort_40.RData")
save(RD_Ddel_SuperShort,file="RD_Ddel_SuperShort_40.RData")
save(RD_Sdup_SuperShort,file="RD_Sdup_SuperShort_40.RData")
save(RD_Ddup_SuperShort,file="RD_Ddup_SuperShort_40.RData")
save(RD_mix_SuperShort,file="RD_mix_SuperShort_40.RData")
save(cpSuperShort.list, file="cpSuperShort_40.RData")
save(locations.SuperShort, file="locations.SuperShort40.RData")

RD_Sdel_short <- res$RD_Sdel_short
RD_Ddel_short <- res$RD_Ddel_short
RD_Sdup_short <- res$RD_Sdup_short
RD_Ddup_short <- res$RD_Ddup_short
RD_mix_short <- res$RD_mix_short
cpshort.list <- res$cpshort.list
locations.short <- res$locations.short
save(RD_Sdel_short,file="RD_Sdel_short_40.RData")
save(RD_Ddel_short,file="RD_Ddel_short_40.RData")
save(RD_Sdup_short,file="RD_Sdup_short_40.RData")
save(RD_Ddup_short,file="RD_Ddup_short_40.RData")
save(RD_mix_short,file="RD_mix_short_40.RData")
save(cpshort.list, file="cpshort_40.RData")
save(locations.short, file="locations.short40.RData")

RD_Sdel_medium <- res$RD_Sdel_medium
RD_Ddel_medium <- res$RD_Ddel_medium
RD_Sdup_medium <- res$RD_Sdup_medium
RD_Ddup_medium <- res$RD_Ddup_medium
RD_mix_medium <- res$RD_mix_medium
cpmedium.list <- res$cpmedium.list
locations.medium <- res$locations.medium
save(RD_Sdel_medium,file="RD_Sdel_medium_40.RData")
save(RD_Ddel_medium,file="RD_Ddel_medium_40.RData")
save(RD_Sdup_medium,file="RD_Sdup_medium_40.RData")
save(RD_Ddup_medium,file="RD_Ddup_medium_40.RData")
save(RD_mix_medium,file="RD_mix_medium_40.RData")
save(cpmedium.list, file="cpmedium_40.RData")
save(locations.medium, file="locations.medium40.RData")

RD_Sdel_long <- res$RD_Sdel_long
RD_Ddel_long <- res$RD_Ddel_long
RD_Sdup_long <- res$RD_Sdup_long
RD_Ddup_long <- res$RD_Ddup_long
RD_mix_long <- res$RD_mix_long
cplong.list <- res$cplong.list
locations.long <- res$locations.long
save(RD_Sdel_long,file="RD_Sdel_long_40.RData")
save(RD_Ddel_long,file="RD_Ddel_long_40.RData")
save(RD_Sdup_long,file="RD_Sdup_long_40.RData")
save(RD_Ddup_long,file="RD_Ddup_long_40.RData")
save(RD_mix_long,file="RD_mix_long_40.RData")
save(cplong.list, file="cplong_40.RData")
save(locations.long, file="locations.long40.RData")





res <-Gen_simu(rd2=rd2, nCluster=3, nCNV=50, Percent=20)

setwd("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation20Percent")
mix.CNVindex.list <- res$mix.CNVindex.list
save(mix.CNVindex.list, file="mix.CNVindex.list.RData")

RD_Sdel_SuperShort <- res$RD_Sdel_SuperShort
RD_Ddel_SuperShort <- res$RD_Ddel_SuperShort
RD_Sdup_SuperShort <- res$RD_Sdup_SuperShort
RD_Ddup_SuperShort <- res$RD_Ddup_SuperShort
RD_mix_SuperShort <- res$RD_mix_SuperShort
cpSuperShort.list <- res$cpSuperShort.list
locations.SuperShort <- res$locations.SuperShort
save(RD_Sdel_SuperShort,file="RD_Sdel_SuperShort_20.RData")
save(RD_Ddel_SuperShort,file="RD_Ddel_SuperShort_20.RData")
save(RD_Sdup_SuperShort,file="RD_Sdup_SuperShort_20.RData")
save(RD_Ddup_SuperShort,file="RD_Ddup_SuperShort_20.RData")
save(RD_mix_SuperShort,file="RD_mix_SuperShort_20.RData")
save(cpSuperShort.list, file="cpSuperShort_20.RData")
save(locations.SuperShort, file="locations.SuperShort20.RData")

RD_Sdel_short <- res$RD_Sdel_short
RD_Ddel_short <- res$RD_Ddel_short
RD_Sdup_short <- res$RD_Sdup_short
RD_Ddup_short <- res$RD_Ddup_short
RD_mix_short <- res$RD_mix_short
cpshort.list <- res$cpshort.list
locations.short <- res$locations.short
save(RD_Sdel_short,file="RD_Sdel_short_20.RData")
save(RD_Ddel_short,file="RD_Ddel_short_20.RData")
save(RD_Sdup_short,file="RD_Sdup_short_20.RData")
save(RD_Ddup_short,file="RD_Ddup_short_20.RData")
save(RD_mix_short,file="RD_mix_short_20.RData")
save(cpshort.list, file="cpshort_20.RData")
save(locations.short, file="locations.short20.RData")

RD_Sdel_medium <- res$RD_Sdel_medium
RD_Ddel_medium <- res$RD_Ddel_medium
RD_Sdup_medium <- res$RD_Sdup_medium
RD_Ddup_medium <- res$RD_Ddup_medium
RD_mix_medium <- res$RD_mix_medium
cpmedium.list <- res$cpmedium.list
locations.medium <- res$locations.medium
save(RD_Sdel_medium,file="RD_Sdel_medium_20.RData")
save(RD_Ddel_medium,file="RD_Ddel_medium_20.RData")
save(RD_Sdup_medium,file="RD_Sdup_medium_20.RData")
save(RD_Ddup_medium,file="RD_Ddup_medium_20.RData")
save(RD_mix_medium,file="RD_mix_medium_20.RData")
save(cpmedium.list, file="cpmedium_20.RData")
save(locations.medium, file="locations.medium20.RData")

RD_Sdel_long <- res$RD_Sdel_long
RD_Ddel_long <- res$RD_Ddel_long
RD_Sdup_long <- res$RD_Sdup_long
RD_Ddup_long <- res$RD_Ddup_long
RD_mix_long <- res$RD_mix_long
cplong.list <- res$cplong.list
locations.long <- res$locations.long
save(RD_Sdel_long,file="RD_Sdel_long_20.RData")
save(RD_Ddel_long,file="RD_Ddel_long_20.RData")
save(RD_Sdup_long,file="RD_Sdup_long_20.RData")
save(RD_Ddup_long,file="RD_Ddup_long_20.RData")
save(RD_mix_long,file="RD_mix_long_20.RData")
save(cplong.list, file="cplong_20.RData")
save(locations.long, file="locations.long20.RData")





Percent=80
Final.output <- NULL
for (Percent in c(100, 80, 60, 40, 20)){
  data=data.frame(Percent=rep(Percent,4),
                  Length=c("SuperShort","Short","Medium","Long"))
  
  for (CNVtype in c("Ddel","Sdel","mix","Sdup","Ddup")){
    
    method <- "FLCNV"
    pre_cluster <- get(load("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\random.Cell.3clusters.RData"))
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_SuperShort",Percent,"_output.RData")))
    cluster <- output$s.hat.best
    # library(funtimes)
    # res <- purity(cluster,pre_cluster)
    # FLCNV_del_SuperShort <- res$pur
    
    library(MixGHD)
    FLCNV_SuperShort <- ARI(cluster,pre_cluster)
    
    
    
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_short",Percent,"_output.RData")))
    cluster <- output$s.hat.best
    # library(funtimes)
    # res <- purity(cluster,pre_cluster)
    # FLCNV_del_short <- res$pur
    FLCNV_short  <- ARI(cluster,pre_cluster)
    
    
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_medium",Percent,"_output.RData")))
    cluster <- output$s.hat.best
    # library(funtimes)
    # res <- purity(cluster,pre_cluster)
    # FLCNV_del_medium <- res$pur
    FLCNV_medium  <- ARI(cluster,pre_cluster)
    
    
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_long",Percent,"_output.RData")))
    cluster <- output$s.hat.best
    # library(funtimes)
    # res <- purity(cluster,pre_cluster)
    # FLCNV_del_long <- res$pur
    FLCNV_long  <- ARI(cluster,pre_cluster)
    
    
    FLCNV <- c(FLCNV_SuperShort, FLCNV_short, FLCNV_medium, FLCNV_long)
    
    
    
    
    #HMMcopy Evaluation of simulation results hierarchical clustering
    method <- "HMMcopy"
    pre_cluster <- get(load("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\random.Cell.3clusters.RData"))
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_SuperShort",Percent,"_output.RData")))
    copy.output <- output$copy.output
    table(output$copy.output[,6])
    clusters <- hclust(dist(t(copy.output)))
    clusterCut <- cutree(clusters, 3)
    pur <- ARI(clusterCut,pre_cluster)
    HMM_Hclus_SuperShort <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(copy.output),3) 
    # HMM_kmeans_del_SuperShort <- purity(kmeans_res$cluster,pre_cluster)$pur
    HMM_kmeans_SuperShort  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_short",Percent,"_output.RData")))
    copy.output <- output$copy.output
    table(output$copy.output[,6])
    clusters <- hclust(dist(t(copy.output)))
    clusterCut <- cutree(clusters, 3)
    pur <- ARI(clusterCut,pre_cluster)
    HMM_Hclus_short <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(copy.output),3) 
    # HMM_kmeans_del_short <- purity(kmeans_res$cluster,pre_cluster)$pur
    HMM_kmeans_short  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_medium",Percent,"_output.RData")))
    copy.output <- output$copy.output
    table(output$copy.output[,6])
    clusters <- hclust(dist(t(copy.output)))
    clusterCut <- cutree(clusters, 3)
    pur <- ARI(clusterCut,pre_cluster)
    HMM_Hclus_medium <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(copy.output),3) 
    # HMM_kmeans_del_medium <- purity(kmeans_res$cluster,pre_cluster)$pur
    HMM_kmeans_medium  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_long",Percent,"_output.RData")))
    copy.output <- output$copy.output
    table(output$copy.output[,6])
    clusters <- hclust(dist(t(copy.output)))
    clusterCut <- cutree(clusters, 3)
    pur <- ARI(clusterCut,pre_cluster)
    HMM_Hclus_long <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(copy.output),3) 
    # HMM_kmeans_del_long <- purity(kmeans_res$cluster,pre_cluster)$pur
    HMM_kmeans_long  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    HMM_Hclus <- c(HMM_Hclus_SuperShort, HMM_Hclus_short, HMM_Hclus_medium, HMM_Hclus_long)
    
    HMM_kmeans <- c(HMM_kmeans_SuperShort, HMM_kmeans_short, HMM_kmeans_medium, HMM_kmeans_long)
    
    HMM_Hclus
    HMM_kmeans
    
    
    
    #SCOPE Evaluation of simulation results hierarchical clustering
    method <- "SCOPE"
    pre_cluster <- get(load("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\random.Cell.3clusters.RData"))
    output <- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                              method,"\\",method,"_",CNVtype,"_SuperShort",Percent,"_alphaoutput.RData")))
    clusters <- hclust(dist(t(output[,1:200])))
    clusterCut <- cutree(clusters, 3)
    table(clusterCut)
    pur <- ARI(clusterCut,pre_cluster)
    SCOPE_Hclus_SuperShort <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(output[,1:200]),3) 
    # SCOPE_kmeans_del_SuperShort <- purity(kmeans_res$cluster,pre_cluster)$pur
    SCOPE_kmeans_SuperShort  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    output <- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                              method,"\\",method,"_",CNVtype,"_short",Percent,"_alphaoutput.RData")))
    clusters <- hclust(dist(t(output[,1:200])))
    clusterCut <- cutree(clusters, 3)
    table(clusterCut)
    pur <- ARI(clusterCut,pre_cluster)
    SCOPE_Hclus_short <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(output[,1:200]),3) 
    # SCOPE_kmeans_del_short <- purity(kmeans_res$cluster,pre_cluster)$pur
    SCOPE_kmeans_short  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    output<- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                             method,"\\",method,"_",CNVtype,"_medium",Percent,"_alphaoutput.RData")))
    clusters <- hclust(dist(t(output[,1:200])))
    clusterCut <- cutree(clusters, 3)
    table(clusterCut)
    pur <- ARI(clusterCut,pre_cluster)
    SCOPE_Hclus_medium <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(output[,1:200]),3) 
    # SCOPE_kmeans_del_medium <- purity(kmeans_res$cluster,pre_cluster)$pur
    SCOPE_kmeans_medium  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    
    output <- get(load(paste0("E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simulation",Percent,"Percent\\",
                              method,"\\",method,"_",CNVtype,"_long",Percent,"_alphaoutput.RData")))
    clusters <- hclust(dist(t(output[,1:200])))
    clusterCut <- cutree(clusters, 3)
    table(clusterCut)
    pur <- ARI(clusterCut,pre_cluster)
    SCOPE_Hclus_long <- pur
    set.seed(2022)
    kmeans_res <- kmeans(t(output[,1:200]),3) 
    # SCOPE_kmeans_del_long <- purity(kmeans_res$cluster,pre_cluster)$pur
    SCOPE_kmeans_long  <- ARI(kmeans_res$cluster,pre_cluster)
    
    
    SCOPE_Hclus <- c(SCOPE_Hclus_SuperShort, SCOPE_Hclus_short, SCOPE_Hclus_medium, SCOPE_Hclus_long)
    SCOPE_kmeans <- c(SCOPE_kmeans_SuperShort, SCOPE_kmeans_short, SCOPE_kmeans_medium, SCOPE_kmeans_long)
    
    SCOPE_Hclus
    SCOPE_kmeans
    
    
    
    
    data1 <- data.frame(
      FLCNV=c(FLCNV_SuperShort, FLCNV_short, FLCNV_medium, FLCNV_long),
      SCOPE_Hclus=c(SCOPE_Hclus_SuperShort, SCOPE_Hclus_short, SCOPE_Hclus_medium, SCOPE_Hclus_long),
      SCOPE_kmeans=c(SCOPE_kmeans_SuperShort, SCOPE_kmeans_short, SCOPE_kmeans_medium, SCOPE_kmeans_long),
      HMM_Hclus=c(HMM_Hclus_SuperShort, HMM_Hclus_short, HMM_Hclus_medium, HMM_Hclus_long),
      HMM_kmeans=c(HMM_kmeans_SuperShort, HMM_kmeans_short, HMM_kmeans_medium, HMM_kmeans_long))
    
    colnames(data1) <- paste0(colnames(data1), "_", CNVtype)
    
    data1[data1 < 0]<- 0
    data <- cbind(data, data1)                  
    
  }  
  Final.output <- rbind(Final.output, data)
  
}
Final.output


# write.csv(Final.output, file="E:\\10XGenomics\\breast_tissue_ABCDE_10k\\Simulation 3 clusters 5 states\\Simu_eval_output_3cluster_ARI.csv")






