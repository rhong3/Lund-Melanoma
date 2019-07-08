library(dplyr)
Phospho <- read.delim("~/Documents/Lund_Melanoma/phospho/MM_Phosho_DIA_lund_Cohort_total_peptides_NORMALIZED.txt")
Phospho$label = paste(Phospho$EG.ModifiedSequence, Phospho$PG.ProteinAccessions, sep='~')
Phospho.cl = Phospho[,-c(1:6, 8)]
rownames(Phospho.cl) = Phospho.cl$label
Phospho.cl = Phospho.cl[, c(1:123)]
Phospho.cly = Phospho.cl[-which(rowMeans(is.na(Phospho.cl)) > 0.7), ]
Phospho.cly$MM790 = rowMeans(Phospho.cly[c('MM790.1.', 'MM790.2.')], na.rm=TRUE)
Phospho.cly$MM807 = rowMeans(Phospho.cly[c('MM807.1.', 'MM807.2.')], na.rm=TRUE)
Phospho.cly$MM808 = Phospho.cly$MM808_LG
Phospho.cly = select(Phospho.cly, -c(MM790.1., MM790.2., MM807.1., MM807.2., MM808_LG))
write.csv(Phospho.cly, "Clean_70_phospho.csv")

library(readxl)
clinical = read_excel('Copy of ClinicalData_144samplesLundMM_Query2018-10-01 (2).xlsx')
rownames(clinical) = toupper(clinical$sample)

clinical.phos = subset(clinical, rownames(clinical) %in% colnames(Phospho.cly))
rownames(clinical.phos)  = toupper(clinical.phos$sample)
Phospho.clz = Phospho.cly[, rownames(clinical.phos)]
clinical.phos = select(clinical.phos, -c(1))
write.table(Phospho.clz, file="Clean_phospho_WNA.tsv", quote=FALSE, sep='\t', col.names = NA)
clinical.phose = clinical.phos[,-which(colMeans(is.na(clinical.phos)) > 0.99)]
rownames(clinical.phose) = rownames(clinical.phos)
write.table(clinical.phose, file="Clean_clinical.tsv", quote=FALSE, sep='\t', col.names = NA)
Phospho.clNA = na.omit(Phospho.clz)
write.table(Phospho.clNA, file="Clean_phospho.tsv", quote=FALSE, sep='\t', col.names = NA)

# For imputed data only
library(dplyr)
Imputed_phos <- read.delim("~/Documents/Lund_Melanoma/phospho/Caret_imputed_data/Imputed_phos.tsv")
names(Imputed_phos)[1] <- "EG.ModifiedSequence"
Phospho <- read.delim("~/Documents/Lund_Melanoma/phospho/MM_Phosho_DIA_lund_Cohort_total_peptides_NORMALIZED.txt")
Phospho = Phospho[, c(1:8)]
Imphos<-merge(x=Imputed_phos,y=Phospho,by="EG.ModifiedSequence")
Imphos$label = paste(Imphos$EG.ModifiedSequence, Imphos$PG.ProteinAccessions, sep='~')
Imphos.cl = Imphos[,-c(1,124:130)]
Imphos.cl$MM790 = rowMeans(Imphos.cl[c('MM790.1.', 'MM790.2.')], na.rm=TRUE)
Imphos.cl$MM807 = rowMeans(Imphos.cl[c('MM807.1.', 'MM807.2.')], na.rm=TRUE)
Imphos.cl$MM808 = Imphos.cl$MM808_LG
Imphos.cl = select(Imphos.cl, -c(MM790.1., MM790.2., MM807.1., MM807.2., MM808_LG))
rownames(Imphos.cl) = Imphos.cl$label
Imphos.cl = Imphos.cl[,-118]

library(readxl)
clinical = read_excel('~/Documents/Lund_Melanoma/phospho/Copy of ClinicalData_144samplesLundMM_Query2018-10-01 (2).xlsx')
rownames(clinical) = toupper(clinical$sample)

clinical.phos = subset(clinical, rownames(clinical) %in% colnames(Imphos.cl))
rownames(clinical.phos)  = toupper(clinical.phos$sample)
Imphos.cl = Imphos.cl[, rownames(clinical.phos)]
clinical.phos = select(clinical.phos, -c(1))
clinical.phose = clinical.phos[,-which(colMeans(is.na(clinical.phos)) > 0.99)]
rownames(clinical.phose) = rownames(clinical.phos)
write.table(clinical.phose, file="Clean_clinical_ip.tsv", quote=FALSE, sep='\t', col.names = NA)
write.table(Imphos.cl, file="Clean_phospho_ip.tsv", quote=FALSE, sep='\t', col.names = NA)
