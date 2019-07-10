# 2 group comparison preparation
library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

# BRAF mut
ica <-read.delim("~/Documents/Lund_Melanoma/Transcriptome/ICA/Gene_level/J_Clean_transcriptome_IC_centroid.txt", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund_Melanoma/Transcriptome/T_wna_clinical.tsv", row.names=1)
data <- read.delim("~/Documents/Lund_Melanoma/Transcriptome/J_Clean_transcriptome.tsv", row.names=1)
wna_clinical = wna_clinical[order(wna_clinical['BRAF.status_V600A'],	wna_clinical['BRAF.status_V600E'],	wna_clinical['BRAF.status_V600K'],	wna_clinical['BRAF.status_WT'],	wna_clinical['BRAF.status_nan']),]
data = data[match(rownames(ica), rownames(data)), match(rownames(wna_clinical), colnames(data))]
categoryB = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['BRAF.status_nan']==1))),
                                                                    rep("WT", length(which(wna_clinical['BRAF.status_WT']==1))),
                                                                    rep("V600K", length(which(wna_clinical['BRAF.status_V600K']==1))),
                                                                    rep("V600E", length(which(wna_clinical['BRAF.status_V600E']==1))),
                                                                    rep("V600A", length(which(wna_clinical['BRAF.status_V600A']==1)))))
core_data <- data

wna_clinical = subset(wna_clinical, wna_clinical$BRAF.status_nan == 0)

B_data = transpose(core_data)
colnames(B_data) <- rownames(core_data)
rownames(B_data) <- colnames(core_data)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_data, keep.rownames = TRUE), rn)
B_data = B_data[wna_clinical, BRAF.status_WT := i.BRAF.status_WT][]
B_data = na.omit(B_data)
GP1 = subset(B_data, B_data$BRAF.status_WT == 1)
GP2 = subset(B_data, B_data$BRAF.status_WT == 0)
GP1.m = GP1$rn
GP2.m = GP2$rn
write.table(GP1.m, file="~/Documents/Lund_Melanoma/Transcriptome/OL/BRAF/G1.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\n')
write.table(GP2.m, file="~/Documents/Lund_Melanoma/Transcriptome/OL/BRAF/G2.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\n')


# Alive/dead binary
ica <-read.delim("~/Documents/Lund_Melanoma/Transcriptome/ICA/Gene_level/J_Clean_transcriptome_IC_centroid.txt", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund_Melanoma/Transcriptome/T_wna_clinical.tsv", row.names=1)
data <- read.delim("~/Documents/Lund_Melanoma/Transcriptome/J_Clean_transcriptome.tsv", row.names=1)
wna_clinical = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
data = data[match(rownames(ica), rownames(data)), match(rownames(wna_clinical), colnames(data))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['Alive.2016.12.05_nan']==1))), 
                                                                    rep("dead.unknown.reason", length(which(wna_clinical['Alive.2016.12.05_dead.unknown.reason']==1))),
                                                                    rep("dead.other.reason", length(which(wna_clinical['Alive.2016.12.05_dead.other.reason']==1))),
                                                                    rep("dead.(likely.melanoma)", length(which(wna_clinical['Alive.2016.12.05_dead..likely.melanoma.']==1))),
                                                                    rep("dead", length(which(wna_clinical['Alive.2016.12.05_dead']==1))),
                                                                    rep("alive", length(which(wna_clinical['Alive.2016.12.05_alive']==1)))))
core_data <- data

wna_clinical = subset(wna_clinical, wna_clinical$Alive.2016.12.05_nan == 0)

B_data = transpose(core_data)
colnames(B_data) <- rownames(core_data)
rownames(B_data) <- colnames(core_data)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_data, keep.rownames = TRUE), rn)
B_data = B_data[wna_clinical, Alive.2016.12.05_alive := i.Alive.2016.12.05_alive][]
B_data = na.omit(B_data)
GP1 = subset(B_data, B_data$Alive.2016.12.05_alive == 1)
GP2 = subset(B_data, B_data$Alive.2016.12.05_alive == 0)
GP1.m = GP1$rn
GP2.m = GP2$rn
write.table(GP1.m, file="~/Documents/Lund_Melanoma/Transcriptome/OL/Survival_binary/G1.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\n')
write.table(GP2.m, file="~/Documents/Lund_Melanoma/Transcriptome/OL/Survival_binary/G2.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\n')

# 5YR Alive/dead binary
ica <-read.delim("~/Documents/Lund_Melanoma/Transcriptome/ICA/Gene_level/J_Clean_transcriptome_IC_centroid.txt", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund_Melanoma/Transcriptome/T_wna_clinical.tsv", row.names=1)
data <- read.delim("~/Documents/Lund_Melanoma/Transcriptome/J_Clean_transcriptome.tsv", row.names=1)
wna_clinical = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
data = data[match(rownames(ica), rownames(data)), match(rownames(wna_clinical), colnames(data))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['Alive.2016.12.05_nan']==1))), 
                                                                    rep("dead.unknown.reason", length(which(wna_clinical['Alive.2016.12.05_dead.unknown.reason']==1))),
                                                                    rep("dead.other.reason", length(which(wna_clinical['Alive.2016.12.05_dead.other.reason']==1))),
                                                                    rep("dead.(likely.melanoma)", length(which(wna_clinical['Alive.2016.12.05_dead..likely.melanoma.']==1))),
                                                                    rep("dead", length(which(wna_clinical['Alive.2016.12.05_dead']==1))),
                                                                    rep("alive", length(which(wna_clinical['Alive.2016.12.05_alive']==1)))))

core_data <- data

wna_clinical = subset(wna_clinical, wna_clinical$Alive.2016.12.05_nan == 0)
B_data = transpose(core_data)
colnames(B_data) <- rownames(core_data)
rownames(B_data) <- colnames(core_data)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_data, keep.rownames = TRUE), rn)
B_data = B_data[wna_clinical, X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years := i.X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years][]
B_data = na.omit(B_data)
GP1 = subset(B_data, B_data$X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years > 0)
GP2 = subset(B_data, B_data$X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years < 0)
GP1.m = GP1$rn
GP2.m = GP2$rn
write.table(GP1.m, file="~/Documents/Lund_Melanoma/Transcriptome/OL/5YR_binary/G1.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\n')
write.table(GP2.m, file="~/Documents/Lund_Melanoma/Transcriptome/OL/5YR_binary/G2.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\n')


