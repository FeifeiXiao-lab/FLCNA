library(SCOPE)
library(DNAcopy)
library(doParallel)

setwd("/home/fqin/Fused.lasso/SCOPE-master/R_old")
files <- list.files()
for (i in files){
  source(i)
}

setwd("/home/fqin/Fused.lasso/simulation_3clusters_5states/Simulation100Percent")
Y_sim <- t(get(load("RD_Ddel_SuperShort_100.RData")))
qcObj <- get(load("qcObj_QC_clear.RData"))

Y_sim1 <- cbind(Y_sim, qcObj$Y[,1:20])+1
Normindex <- 201:220
ref_sim <- qcObj$ref
QCmetric_raw <- qcObj$QCmetric
setwd("/home/fqin/Fused.lasso/simulation_3clusters_5states/Simulation100Percent/SCOPE")

sampname_raw <- colnames(Y_sim1)
#qcObj <- perform_qc(Y_raw = Y_sim1, 
#                    sampname_raw = sampname_raw, ref_raw = ref_sim, 
#                    QCmetric_raw = QCmetric_raw, minCountQC = 0)

Y <- Y_sim1
gc_qc <- qcObj$ref$gc
ref <- ref_sim

# first-pass CODEX2 run with no latent factors
normObj.sim <- normalize_codex2_ns_noK(Y_qc = Y,
                                       gc_qc = gc_qc,
                                       norm_index = Normindex)

# Ploidy initialization
ploidy.sim <- initialize_ploidy(Y = Y, Yhat = normObj.sim$Yhat, ref = ref)

normObj.scope.sim <- normalize_scope_foreach(Y_qc = Y, gc_qc = gc_qc,
                                             K = 1, ploidyInt = ploidy.sim,
                                             norm_index = Normindex, T = 1:5,
                                             beta0 = normObj.sim$beta.hat)

save(normObj.scope.sim, file="SCOPE_Ddel_SuperShort100_output.RData")