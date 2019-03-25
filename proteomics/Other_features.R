# Other features BRAF

library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
ica <- read.delim("~/Documents/Lund/proteomics/ICA/Gene_level/J_Clean_proteomics_IC_centroid.txt", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[order(wna_clinical['BRAF.status_V600A'],	wna_clinical['BRAF.status_V600E'],	wna_clinical['BRAF.status_V600K'],	wna_clinical['BRAF.status_WT'],	wna_clinical['BRAF.status_nan']),]
proteomics = proteomics[match(rownames(ica), rownames(proteomics)), match(rownames(wna_clinical), colnames(proteomics))]
categoryB = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['BRAF.status_nan']==1))),
                                                                    rep("WT", length(which(wna_clinical['BRAF.status_WT']==1))),
                                                                    rep("V600K", length(which(wna_clinical['BRAF.status_V600K']==1))),
                                                                    rep("V600E", length(which(wna_clinical['BRAF.status_V600E']==1))),
                                                                    rep("V600A", length(which(wna_clinical['BRAF.status_V600A']==1)))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)
pdf("~/Documents/Lund/proteomics/BRAF/Extended_Core_gene_BRAF_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = F, cluster_rows = T, annotation_col = categoryB, fontsize_col = 6, main = 'Extended Core_gene vs BRAF')
dev.off()

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, BRAF.status_WT := i.BRAF.status_WT][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/BRAF/B_temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/BRAF/B_temp.csv", row.names=1)
WTBRAF = subset(B_prot, B_prot$BRAF.status_WT == 1)
WTBRAF = WTBRAF[,-c(26)]
MutBRAF = subset(B_prot, B_prot$BRAF.status_WT == 0)
MutBRAF = MutBRAF[,-c(26)]
t.test(WTBRAF, MutBRAF)
x = list('Overall'=t.test(WTBRAF, MutBRAF)$p.value)
for (i in 1:25){
  p = t.test(WTBRAF[,i], MutBRAF[,i])$p.value
  x[[colnames(WTBRAF)[i]]] = p
  if (p < 0.1){
    print(colnames(WTBRAF)[i], print(p))
  }
}
xdf=as.data.frame(x)
rownames(xdf) = c('T-Test p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/BRAF/T-Test.csv", row.names=TRUE)


# Alive/Dead binary

library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
ica <- read.delim("~/Documents/Lund/proteomics/ICA/Gene_level/J_Clean_proteomics_IC_centroid.txt", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
proteomics = proteomics[match(rownames(ica), rownames(proteomics)), match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['Alive.2016.12.05_nan']==1))), 
                                                                    rep("dead.unknown.reason", length(which(wna_clinical['Alive.2016.12.05_dead.unknown.reason']==1))),
                                                                    rep("dead.other.reason", length(which(wna_clinical['Alive.2016.12.05_dead.other.reason']==1))),
                                                                    rep("dead.(likely.melanoma)", length(which(wna_clinical['Alive.2016.12.05_dead..likely.melanoma.']==1))),
                                                                    rep("dead", length(which(wna_clinical['Alive.2016.12.05_dead']==1))),
                                                                    rep("alive", length(which(wna_clinical['Alive.2016.12.05_alive']==1)))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)

wna_clinical = subset(wna_clinical, wna_clinical$Alive.2016.12.05_nan == 0)

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, Alive.2016.12.05_alive := i.Alive.2016.12.05_alive][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/Survival/S_temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/Survival/S_temp.csv", row.names=1)
WTBRAF = subset(B_prot, B_prot$Alive.2016.12.05_alive == 1)
WTBRAF = WTBRAF[,-c(26)]
MutBRAF = subset(B_prot, B_prot$Alive.2016.12.05_alive == 0)
MutBRAF = MutBRAF[,-c(26)]
t.test(WTBRAF, MutBRAF)
x = list('Overall'=t.test(WTBRAF, MutBRAF)$p.value)
for (i in 1:25){
  p = t.test(WTBRAF[,i], MutBRAF[,i])$p.value
  x[[colnames(WTBRAF)[i]]] = p
  if (p < 0.1){
    print(colnames(WTBRAF)[i], print(p))
  }
}
xdf=as.data.frame(x)
rownames(xdf) = c('T-Test p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/Survival/Binary_T-Test.csv", row.names=TRUE)

# 5-yr binary
library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
ica <- read.delim("~/Documents/Lund/proteomics/ICA/Gene_level/J_Clean_proteomics_IC_centroid.txt", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
proteomics = proteomics[match(rownames(ica), rownames(proteomics)), match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['Alive.2016.12.05_nan']==1))), 
                                                                    rep("dead.unknown.reason", length(which(wna_clinical['Alive.2016.12.05_dead.unknown.reason']==1))),
                                                                    rep("dead.other.reason", length(which(wna_clinical['Alive.2016.12.05_dead.other.reason']==1))),
                                                                    rep("dead.(likely.melanoma)", length(which(wna_clinical['Alive.2016.12.05_dead..likely.melanoma.']==1))),
                                                                    rep("dead", length(which(wna_clinical['Alive.2016.12.05_dead']==1))),
                                                                    rep("alive", length(which(wna_clinical['Alive.2016.12.05_alive']==1)))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)

wna_clinical = subset(wna_clinical, wna_clinical$Alive.2016.12.05_nan == 0)

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years := i.X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/Survival/5S_temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/Survival/5S_temp.csv", row.names=1)
WTBRAF = subset(B_prot, B_prot$X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years > 0)
WTBRAF = WTBRAF[,-c(26)]
MutBRAF = subset(B_prot, X5.year.survival.from.primary.diagnosis...Date.death.date.prim.diagn..1825.days...Days.differing.from.5.years < 0)
MutBRAF = MutBRAF[,-c(26)]
t.test(WTBRAF, MutBRAF)
x = list('Overall'=t.test(WTBRAF, MutBRAF)$p.value)
for (i in 1:25){
  p = t.test(WTBRAF[,i], MutBRAF[,i])$p.value
  x[[colnames(WTBRAF)[i]]] = p
  if (p < 0.1){
    print(colnames(WTBRAF)[i], print(p))
  }
}
xdf=as.data.frame(x)
rownames(xdf) = c('T-Test p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/Survival/5YR_Binary_T-Test.csv", row.names=TRUE)


# dis.stage #
library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[order(wna_clinical['dis.stage']),]
proteomics = proteomics[, match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("1", length(which(wna_clinical['dis.stage']==1))), 
                                                                    rep("3", length(which(wna_clinical['dis.stage']==3))),
                                                                    rep("4", length(which(wna_clinical['dis.stage']==4))),
                                                                    rep("NA", 4)))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)
pdf("~/Documents/Lund/proteomics/dis_stage/Extended_Core_gene_dis_stage_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = F, cluster_rows = T, annotation_col = categoryS, fontsize_col = 6, main = 'Extended Core_gene vs dis_stage')
dev.off()

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, dis.stage := i.dis.stage][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/dis_stage/temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/dis_stage/temp.csv", row.names=1)

x = list()
for (i in 1:25){
  ff = aov(B_prot[, i] ~ dis.stage, data = B_prot)
  p = summary(ff)[[1]][["Pr(>F)"]][1]
  x[[colnames(B_prot)[i]]] = p
  if (p < 0.1){
    print(colnames(B_prot)[i], print(p))
  }
}

xdf=as.data.frame(x)
rownames(xdf) = c('ANOVA p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/dis_stage/ANOVA.csv", row.names=TRUE)

# prim breslow class
library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[!is.na(wna_clinical$prim.breslow.class),]
wna_clinical = wna_clinical[order(wna_clinical['prim.breslow.class']),]
proteomics = proteomics[, match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("1", length(which(wna_clinical['prim.breslow.class']==1))), 
                                                                    rep("2", length(which(wna_clinical['prim.breslow.class']==2))),
                                                                    rep("3", length(which(wna_clinical['prim.breslow.class']==3))),
                                                                    rep("4", length(which(wna_clinical['prim.breslow.class']==4)))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)
pdf("~/Documents/Lund/proteomics/prim_breslow_class/Extended_Core_gene_prim_breslow_class_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = F, cluster_rows = T, annotation_col = categoryS, fontsize_col = 6, main = 'Extended Core_gene vs prim_breslow_class')
dev.off()

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, prim.breslow.class := i.prim.breslow.class][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/prim_breslow_class/temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/prim_breslow_class/temp.csv", row.names=1)

x = list()
for (i in 1:25){
  ff = aov(B_prot[, i] ~ prim.breslow.class, data = B_prot)
  p = summary(ff)[[1]][["Pr(>F)"]][1]
  x[[colnames(B_prot)[i]]] = p
  if (p < 0.1){
    print(colnames(B_prot)[i], print(p))
  }
}

xdf=as.data.frame(x)
rownames(xdf) = c('ANOVA p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/prim_breslow_class/ANOVA.csv", row.names=TRUE)

# clark
library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[!is.na(wna_clinical$clark),]
wna_clinical = wna_clinical[order(wna_clinical['clark']),]
proteomics = proteomics[, match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("1", length(which(wna_clinical['clark']==1))), 
                                                                    rep("2", length(which(wna_clinical['clark']==2))),
                                                                    rep("3", length(which(wna_clinical['clark']==3))),
                                                                    rep("4", length(which(wna_clinical['clark']==4))),
                                                                    rep("5", length(which(wna_clinical['clark']==5)))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)
pdf("~/Documents/Lund/proteomics/clark/Extended_Core_gene_clark_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = F, cluster_rows = T, annotation_col = categoryS, fontsize_col = 6, main = 'Extended Core_gene vs clark')
dev.off()

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, clark := i.clark][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/clark/temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/clark/temp.csv", row.names=1)

x = list()
for (i in 1:25){
  ff = aov(B_prot[, i] ~ clark, data = B_prot)
  p = summary(ff)[[1]][["Pr(>F)"]][1]
  x[[colnames(B_prot)[i]]] = p
  if (p < 0.1){
    print(colnames(B_prot)[i], print(p))
  }
}

xdf=as.data.frame(x)
rownames(xdf) = c('ANOVA p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/clark/ANOVA.csv", row.names=TRUE)

#clin class
library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[!is.na(wna_clinical$clin.class),]
wna_clinical = wna_clinical[(wna_clinical$clin.class != ''),]
wna_clinical = wna_clinical[order(wna_clinical['clin.class']),]
proteomics = proteomics[, match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("0", length(which(wna_clinical['clin.class']==0))),
                                                                    rep("1", length(which(wna_clinical['clin.class']==1))), 
                                                                    rep("2", length(which(wna_clinical['clin.class']==2))),
                                                                    rep("3", length(which(wna_clinical['clin.class']==3))),
                                                                    rep("4", length(which(wna_clinical['clin.class']==4))),
                                                                    rep("5", length(which(wna_clinical['clin.class']==5))),
                                                                    rep("NM", length(which(wna_clinical['clin.class']=='NM')))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)
pdf("~/Documents/Lund/proteomics/clin_class/Extended_Core_gene_clin_class_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = F, cluster_rows = T, annotation_col = categoryS, fontsize_col = 6, main = 'Extended Core_gene vs clin_class')
dev.off()

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, clin.class := i.clin.class][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/clin_class/temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/clin_class/temp.csv", row.names=1)

x = list()
for (i in 1:25){
  ff = aov(B_prot[, i] ~ clin.class, data = B_prot)
  p = summary(ff)[[1]][["Pr(>F)"]][1]
  x[[colnames(B_prot)[i]]] = p
  if (p < 0.1){
    print(colnames(B_prot)[i], print(p))
  }
}

xdf=as.data.frame(x)
rownames(xdf) = c('ANOVA p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/clin_class/ANOVA.csv", row.names=TRUE)

# prim site
library("stats")
library('ggplot2')
library('pheatmap')
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[!is.na(wna_clinical$prim.site),]
wna_clinical = wna_clinical[(wna_clinical$prim.site != ''),]
wna_clinical = wna_clinical[order(wna_clinical['prim.site']),]
proteomics = proteomics[, match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("1", length(which(wna_clinical['prim.site']==1))), 
                                                                    rep("2", length(which(wna_clinical['prim.site']==2))),
                                                                    rep("3", length(which(wna_clinical['prim.site']==3))),
                                                                    rep("4", length(which(wna_clinical['prim.site']==4))),
                                                                    rep("5", length(which(wna_clinical['prim.site']==5))),
                                                                    rep("f", length(which(wna_clinical['prim.site']=='f')))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)
pdf("~/Documents/Lund/proteomics/prim_site/Extended_Core_gene_prim_site_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = F, cluster_rows = T, annotation_col = categoryS, fontsize_col = 6, main = 'Extended Core_gene vs prim_site')
dev.off()

B_proteomics = transpose(core_proteomics)
colnames(B_proteomics) <- rownames(core_proteomics)
rownames(B_proteomics) <- colnames(core_proteomics)
setDT(wna_clinical, keep.rownames = TRUE)
setkey(setDT(B_proteomics, keep.rownames = TRUE), rn)
B_proteomics = B_proteomics[wna_clinical, prim.site := i.prim.site][]
B_proteomics = na.omit(B_proteomics)
write.csv(B_proteomics, file = "~/Documents/Lund/proteomics/prim_site/temp.csv", row.names=FALSE)
B_prot <- read.csv("~/Documents/Lund/proteomics/prim_site/temp.csv", row.names=1)

x = list()
for (i in 1:25){
  ff = aov(B_prot[, i] ~ prim.site, data = B_prot)
  p = summary(ff)[[1]][["Pr(>F)"]][1]
  x[[colnames(B_prot)[i]]] = p
  if (p < 0.1){
    print(colnames(B_prot)[i], print(p))
  }
}

xdf=as.data.frame(x)
rownames(xdf) = c('ANOVA p-value')
oxdf = transpose(xdf)
colnames(oxdf) <- rownames(xdf)
rownames(oxdf) <- colnames(xdf)
write.csv(oxdf, file = "~/Documents/Lund/proteomics/prim_site/ANOVA.csv", row.names=TRUE)

