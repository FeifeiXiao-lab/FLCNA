library("readxl")
setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\SC CNV Data\\SRP114962\\Data\\GSEA")
KTN126 <- read_excel("GSEA KTN126 output.xlsx", skip = 9)
KTN126 <- KTN126[1:20, 1:7]
KTN126 <- KTN126[,c(1:2)]
colnames(KTN126) <- c("SetName", "SetSize")
KTN126$Index <- 1
KTN126$Patients <- "KTN126"

library(openxlsx)
write.table(KTN126, file="Network\\KTN126.tsv", quote = F, row.names = F, sep="\t")


KTN129 <- read_excel("GSEA KTN129 output.xlsx", skip = 9)
KTN129 <- KTN129[1:19, 1:7]
KTN129 <- KTN129[,c(1:2)]
colnames(KTN129) <- c("SetName", "SetSize")
KTN129$Index <- 1
KTN129$Patients <- "KTN129"

KTN302 <- read_excel("GSEA KTN302 output.xlsx", skip = 9)
KTN302 <- KTN302[1:13, 1:7]
KTN302 <- KTN302[,c(1:2)]
colnames(KTN302) <- c("SetName", "SetSize")
KTN302$Index <- 1
KTN302$Patients <- "KTN302"

KTN <- merge(KTN126, KTN129, by.x="SetName", by.y="SetName", all.x=T, all.y=T)
KTN <- merge(KTN, KTN302, by.x="SetName", by.y="SetName", all.x=T, all.y=T)



KTN126_neg <-  read.table("KTN126.GseaPreranked.1678049944653//gsea_report_for_na_neg_1678049944653.tsv", sep="\t", header = T)[,-3]
KTN126_pos <-  read.table("KTN126.GseaPreranked.1678049944653//gsea_report_for_na_pos_1678049944653.tsv", sep="\t", header = T)[,-3]
KTN126 <- data.frame(NAME=unique(c(KTN126_neg$NAME, KTN126_pos$NAME)), INDEX=1)
write.csv(KTN126, file="KTN126.Geneindex.csv", row.names=F)

KTN129_neg <-  read.table("KTN129.GseaPreranked.1678050271482//gsea_report_for_na_neg_1678050271482.tsv", sep="\t", header = T)[,-3]
KTN129_neg$Index <- 1
colnames(KTN129_neg) <- c("NAME", "GS<br> follow link to MSigDB", "SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX", "LEADING EDGE", "X", "Index")
write.table(KTN129_neg, file="KTN129.GseaPreranked.1678050271482//gsea_report_for_na_neg_1678050271482_new.tsv", quote = F, sep="\t", row.names = F)

KTN129_pos <-  read.table("KTN129.GseaPreranked.1678050271482//gsea_report_for_na_pos_1678050271482.tsv", sep="\t", header = T)[,-3]
KTN129_pos$Index <- 1
colnames(KTN129_pos) <- c("NAME", "GS<br> follow link to MSigDB",  "SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX", "LEADING EDGE", "X", "Index")
write.table(KTN129_pos, file="KTN129.GseaPreranked.1678050271482//gsea_report_for_na_pos_1678050271482_new.tsv", quote = F, sep="\t", row.names = F)


KTN302_neg <-  read.table("KTN302.GseaPreranked.1678050226304//gsea_report_for_na_neg_1678050226304.tsv", sep="\t", header = T)[,-3]
KTN302_neg$Index <- 1
colnames(KTN302_neg) <- c("NAME", "GS<br> follow link to MSigDB", "SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX", "LEADING EDGE", "X", "Index")
write.table(KTN302_neg, file="KTN302.GseaPreranked.1678050226304//gsea_report_for_na_neg_1678050226304_new.tsv", quote = F, sep="\t", row.names = F)

KTN302_pos <-  read.table("KTN302.GseaPreranked.1678050226304//gsea_report_for_na_pos_1678050226304.tsv", sep="\t", header = T)[,-3]
KTN302_pos$Index <- 1
colnames(KTN302_pos) <- c("NAME", "GS<br> follow link to MSigDB", "SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX", "LEADING EDGE", "X", "Index")
write.table(KTN302_pos, file="KTN302.GseaPreranked.1678050226304//gsea_report_for_na_pos_1678050226304_new.tsv", quote = F, sep="\t", row.names = F)


