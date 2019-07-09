Survivalmd = 'DSS.DAYS.v2...days.from.sample.collection..surgery..to.death.or.censoring'
md = 'Surgery'
# OS.DAYS.....days.from.primary.diagnosis.to.death.or.censoring
# DSS.DAYS.v1.days.from.first.metastasis.to.death.or.censoring
# DSS.DAYS.v2...days.from.sample.collection..surgery..to.death.or.censoring

# Survival based on clusters
library("pheatmap")
library("ggplot2")

core_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1')
ica <-read.csv("~/Documents/Lund_Melanoma/phospho/ICA/Gene_phospho_ip_IC_Centroid.csv", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund_Melanoma/phospho/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund_Melanoma/phospho/J_Clean_phospho_ip.tsv", row.names=1)
wna_clinical = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
proteomics = proteomics[match(rownames(ica), rownames(proteomics)), match(rownames(wna_clinical), colnames(proteomics))]
categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['Alive.2016.12.05_nan']==1))), 
                                                                    rep("dead.unknown.reason", length(which(wna_clinical['Alive.2016.12.05_dead.unknown.reason']==1))),
                                                                    rep("dead.other.reason", length(which(wna_clinical['Alive.2016.12.05_dead.other.reason']==1))),
                                                                    rep("dead.(likely.melanoma)", length(which(wna_clinical['Alive.2016.12.05_dead..likely.melanoma.']==1))),
                                                                    rep("dead", length(which(wna_clinical['Alive.2016.12.05_dead']==1))),
                                                                    rep("alive", length(which(wna_clinical['Alive.2016.12.05_alive']==1)))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% core_genes)
pdf("~/Documents/Lund_Melanoma/phospho/Cluster_Core_gene_survival_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = T, cluster_rows = T, annotation_col = categoryS, fontsize_col = 6, main = 'Cluster_Core_gene vs 5yr Survival')
dev.off()

# BRAF WT vs Mutation
# OS.DAYS.....days.from.primary.diagnosis.to.death.or.censoring
# DSS.DAYS.v1.days.from.first.metastasis.to.death.or.censoring
# DSS.DAYS.v2...days.from.sample.collection..surgery..to.death.or.censoring
library("survival")
library("survminer")
BRAF_clinical = subset(wna_clinical, wna_clinical$BRAF.status_nan == 0)
BRAF_clinical = subset(BRAF_clinical, BRAF_clinical$Alive.2016.12.05_nan == 0)
BRAF_clinical$censor = abs(BRAF_clinical$Alive.2016.12.05_alive - 2)
BRAF_clinical = BRAF_clinical[!is.na(BRAF_clinical[[Survivalmd]]),]
fit <- survfit(Surv(as.numeric(BRAF_clinical[[Survivalmd]]), BRAF_clinical$censor) ~ BRAF_clinical$BRAF.status_WT, data = BRAF_clinical)
print(fit)
summary(fit)$table
ggsurvplot(fit, title = paste('BRAF survival', md), pval = TRUE, conf.int = TRUE, risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

# Core genes survival
library("stats")
kms = kmeans(t(core_proteomics), 2, iter.max = 100, nstart = 1,
       algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                     "MacQueen"), trace=FALSE)

Core_clinical = wna_clinical
Core_clinical$cluster = kms$cluster
Core_clinical = subset(Core_clinical, Core_clinical$Alive.2016.12.05_nan == 0)
Core_clinical$censor = abs(Core_clinical$Alive.2016.12.05_alive - 2)
Core_clinical = Core_clinical[!is.na(Core_clinical[[Survivalmd]]),]
fit <- survfit(Surv(as.numeric(Core_clinical[[Survivalmd]]), Core_clinical$censor) ~ Core_clinical$cluster, data = Core_clinical)
print(fit)
summary(fit)$table
ggsurvplot(fit, title = paste('Core-gene survival', md), pval = TRUE, conf.int = TRUE, risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


# IC ranked-gene survival
# OS.DAYS.....days.from.primary.diagnosis.to.death.or.censoring
# DSS.DAYS.v1.days.from.first.metastasis.to.death.or.censoring
# DSS.DAYS.v2...days.from.sample.collection..surgery..to.death.or.censoring

ica$survival = ica$X84
ica = ica[order(ica$survival),]
wna_clinical <- read.delim("~/Documents/Lund_Melanoma/phospho/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund_Melanoma/phospho/J_Clean_phospho_ip.tsv", row.names=1)
wna_clinicalS = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
proteomics = proteomics[match(rownames(ica), rownames(proteomics)), match(rownames(wna_clinical), colnames(proteomics))]

proteomics = data.matrix(proteomics[-c(26:1568),])
kms = kmeans(t(proteomics), 2, iter.max = 100, nstart = 1,
             algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                           "MacQueen"), trace=FALSE)

IC_clinical = wna_clinical
IC_clinical$cluster = kms$cluster
IC_clinical = subset(IC_clinical, IC_clinical$Alive.2016.12.05_nan == 0)
IC_clinical$censor = abs(IC_clinical$Alive.2016.12.05_alive - 2)
IC_clinical = IC_clinical[!is.na(IC_clinical[[Survivalmd]]),]
fit <- survfit(Surv(as.numeric(IC_clinical[[Survivalmd]]), IC_clinical$censor) ~ IC_clinical$cluster, data = IC_clinical)
print(fit)
summary(fit)$table
  ggsurvplot(fit, title = paste('IC-84 survival', md), pval = TRUE, conf.int = TRUE, risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


#Extended Core genes
# OS.DAYS.....days.from.primary.diagnosis.to.death.or.censoring
# DSS.DAYS.v1.days.from.first.metastasis.to.death.or.censoring
# DSS.DAYS.v2...days.from.sample.collection..surgery..to.death.or.censoring

library("survival")
library("survminer")
library("stats")
library('ggplot2')
library('pheatmap')

Ecore_genes = c('BRAF', 'NRAS', 'TP53', 'NF1', 'CDKN2A', 'ARID2', 'PTEN', 'PPP6C', 'RAC1', 'IDH1', 'DDX3X', 'MAP2K1', 'RB1', 'CTNNB1', 'CASP8', 'PCDHGA1', 'SERPINB1', 'IRF7', 'HRAS', 'PTPN11',
                'ITGA4', 'FAM113B', 'MSR1', 'RPS27', 'SIRPB1', 'MRPS31', 'NOTCH2NL', 'KNSTRN', 'ZFX', 'RAPGEFS', 'RCAN2', 'PPIAL4G', 'ACD', 'WDR12', 'COL9A2', 'STK19', 'CCDC28A', 'LRRC37A3', 'OXA1L',
                'NDUFB9', 'EMG1', 'TMEM216', 'RQCD1', 'TBC1D3B', 'GNAI2', 'B2M', 'FAM58A', 'C3orf71')
wna_clinical = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
proteomics = proteomics[match(rownames(ica), rownames(proteomics)), match(rownames(wna_clinical), colnames(proteomics))]
category = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['Alive.2016.12.05_nan']==1))), 
                                                                    rep("dead.unknown.reason", length(which(wna_clinical['Alive.2016.12.05_dead.unknown.reason']==1))),
                                                                    rep("dead.other.reason", length(which(wna_clinical['Alive.2016.12.05_dead.other.reason']==1))),
                                                                    rep("dead.(likely.melanoma)", length(which(wna_clinical['Alive.2016.12.05_dead..likely.melanoma.']==1))),
                                                                    rep("dead", length(which(wna_clinical['Alive.2016.12.05_dead']==1))),
                                                                    rep("alive", length(which(wna_clinical['Alive.2016.12.05_alive']==1)))))
core_proteomics <- subset(proteomics, rownames(proteomics) %in% Ecore_genes)
pdf("~/Documents/Lund_Melanoma/phospho/Cluster_Extended_Core_gene_survival_HM.pdf", width = 15, paper = 'a4r')
pheatmap(core_proteomics, cluster_cols = T, cluster_rows = T, annotation_col = category, fontsize_col = 6, main = 'Extended Core_gene vs 5yr Survival')
dev.off()

kms = kmeans(t(core_proteomics), 2, iter.max = 100, nstart = 1,
             algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                           "MacQueen"), trace=FALSE)

Core_clinical = wna_clinical
Core_clinical$cluster = kms$cluster
Core_clinical = subset(Core_clinical, Core_clinical$Alive.2016.12.05_nan == 0)
Core_clinical$censor = abs(Core_clinical$Alive.2016.12.05_alive - 2)
Core_clinical = Core_clinical[!is.na(Core_clinical[[Survivalmd]]),]
fit <- survfit(Surv(as.numeric(Core_clinical[[Survivalmd]]), Core_clinical$censor) ~ Core_clinical$cluster, data = Core_clinical)
print(fit)
summary(fit)$table
ggsurvplot(fit, title = paste('Extended survival', md), pval = TRUE, conf.int = TRUE, risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))




