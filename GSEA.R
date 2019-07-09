ICX = 'X4'
ICM = 'IC4'
ICN = 4
Var = 'Local'
#########################################
library('org.Hs.eg.db')
# source("https://bioconductor.org/biocLite.R")
# biocLite("reactome.db")
library("reactome.db", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("ReactomePA", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("fgsea", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
J_Clean_proteomics_IC_centroid <- read.csv("~/Documents/Lund_Melanoma/phospho/ICA/Gene_phospho_ip_IC_Centroid.csv", row.names=1)
Gene_order = J_Clean_proteomics_IC_centroid[order(J_Clean_proteomics_IC_centroid[[ICN]]),]
ENTREZID = mapIds(org.Hs.eg.db, row.names(Gene_order), 'ENTREZID', 'SYMBOL')
Gene_order.ICN <- setNames(as.numeric(Gene_order[[ICX]]), unname(ENTREZID))
my_pathways <- reactomePathways(names(Gene_order.ICN))
summary(sapply(my_pathways, length))
fgsea_reactome <- fgsea(pathways = my_pathways, 
                        stats = Gene_order.ICN,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)
fgsea_reactome <- na.omit(fgsea_reactome[order(pval), ])
fgsea_reactome$leadingEdge = as.character(fgsea_reactome$leadingEdge)
fgsea_reactome.sig = fgsea_reactome[fgsea_reactome$padj < 0.01,]
write.csv(fgsea_reactome, file = paste("~/Documents/Lund_Melanoma/phospho/ICA/GSEA_ip/", ICM, ".csv"))
# pdf("~/Documents/Lund_Melanoma/proteomics/ICA/0403ICA/GSEA/Most_IC8_sig.pdf", paper = 'letter')
# plotEnrichment(my_pathways[[head(fgsea_reactome, 1)$pathway]], Gene_order.8) + labs(title=head(fgsea_reactome, 1)$pathway)
# dev.off()

topPathwaysUp <- fgsea_reactome.sig[ES > 0][head(order(padj), n=20), pathway]
topPathwaysDown <- fgsea_reactome.sig[ES < 0][head(order(padj), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf(paste("~/Documents/Lund_Melanoma/phospho/ICA/GSEA_ip/", ICM, "_sig.pdf"), width = 15, paper = 'a4r')
plotGseaTable(my_pathways[topPathways], Gene_order.ICN, fgsea_reactome, 
              gseaParam = 0.5)
dev.off()

###################################################
library("pheatmap", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
# Heatmap
ica <-read.csv("~/Documents/Lund_Melanoma/phospho/ICA/Gene_phospho_ip_IC_Centroid.csv", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund_Melanoma/phospho/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund_Melanoma/phospho/J_Clean_phospho_ip.tsv", row.names=1)
ica$survival = ica[[ICX]]
ica = ica[order(ica$survival),]
wna_clinicalC = wna_clinical[order(wna_clinical['clin.class.det_ALM'],	wna_clinical['clin.class.det_LMM'],	wna_clinical['clin.class.det_Mucosal'],	wna_clinical['clin.class.det_NM'],	wna_clinical['clin.class.det_Other'],	wna_clinical['clin.class.det_SSM'],	wna_clinical['clin.class.det_Unknownprimary']),]
wna_clinicalS = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]
wna_clinicalST = wna_clinical[order(wna_clinical['stage_General'],	wna_clinical['stage_Local'],	wna_clinical['stage_Regional'],	wna_clinical['stage_In.transit']),]
wna_clinicalL = wna_clinical[order(wna_clinical['local_Cutaneous'],	wna_clinical['local_Lymph.node'],	wna_clinical['local_Subcutaneous'],	wna_clinical['local_Visceral']),]
wna_clinicalB = wna_clinical[order(wna_clinical['BRAF.status_V600A'],	wna_clinical['BRAF.status_V600E'],	wna_clinical['BRAF.status_V600K'],	wna_clinical['BRAF.status_WT']),]

wna_clinical = wna_clinicalL

GSEA_proteomics = proteomics[match(rownames(ica), rownames(proteomics)), match(rownames(wna_clinical), colnames(proteomics))]
categoryC = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['clin.class.det_nan']==1))), 
                                                                    rep("Unknownprimary", length(which(wna_clinical['clin.class.det_Unknownprimary']==1))),
                                                                    rep("SSM", length(which(wna_clinical['clin.class.det_SSM']==1))),
                                                                    rep("Other", length(which(wna_clinical['clin.class.det_Other']==1))),
                                                                    rep("NM", length(which(wna_clinical['clin.class.det_NM']==1))),
                                                                    rep("Mucosal", length(which(wna_clinical['clin.class.det_Mucosal']==1))),
                                                                   rep("LMM", length(which(wna_clinical['clin.class.det_LMM']==1))),
                                                                   rep("ALM", length(which(wna_clinical['clin.class.det_ALM']==1)))))

categoryS = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['Alive.2016.12.05_nan']==1))), 
                                                                   rep("dead.unknown.reason", length(which(wna_clinical['Alive.2016.12.05_dead.unknown.reason']==1))),
                                                                   rep("dead.other.reason", length(which(wna_clinical['Alive.2016.12.05_dead.other.reason']==1))),
                                                                   rep("dead.(likely.melanoma)", length(which(wna_clinical['Alive.2016.12.05_dead..likely.melanoma.']==1))),
                                                                   rep("dead", length(which(wna_clinical['Alive.2016.12.05_dead']==1))),
                                                                   rep("alive", length(which(wna_clinical['Alive.2016.12.05_alive']==1)))))

categoryST = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['stage_nan']==1))), 
                                                                    rep("stage.in.transit", length(which(wna_clinical['stage_In.transit']==1))),
                                                                    rep("stage.regional", length(which(wna_clinical['stage_Regional']==1))),
                                                                    rep("stage.local", length(which(wna_clinical['stage_Local']==1))),
                                                                    rep("stage.general", length(which(wna_clinical['stage_General']==1)))))

categoryL = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['local_nan']==1))), 
                                                                     rep("local.Visceral", length(which(wna_clinical['local_Visceral']==1))),
                                                                     rep("local.Subcutaneous", length(which(wna_clinical['local_Subcutaneous']==1))),
                                                                     rep("local.Lymph.node", length(which(wna_clinical['local_Lymph.node']==1))),
                                                                     rep("local.Curaneous", length(which(wna_clinical['local_Cutaneous']==1)))))

categoryB = data.frame(row.names=rownames(wna_clinical), category=c(rep("NA", length(which(wna_clinical['BRAF.status_nan']==1))), 
                                                                    rep("BRAF.status.WT", length(which(wna_clinical['BRAF.status_WT']==1))),
                                                                    rep("BRAF.status.V600K", length(which(wna_clinical['BRAF.status_V600K']==1))),
                                                                    rep("BRAF.status.V600E", length(which(wna_clinical['BRAF.status_V600E']==1))),
                                                                    rep("BRAF.status.V600A", length(which(wna_clinical['BRAF.status_V600A']==1)))))


GSEA_proteomics = data.matrix(GSEA_proteomics[-c(26:1568),])
GSEA_proteomics = subset(GSEA_proteomics, nchar(as.character(rownames(GSEA_proteomics))) <= 10)
pdf(paste("~/Documents/Lund_Melanoma/phospho/ICA/GSEA_ip/Cluster_", ICM, Var, "_HM.pdf"), width = 15, paper = 'a4r')
pheatmap(GSEA_proteomics, cluster_cols = T, cluster_rows = F, annotation_col = categoryL, fontsize_col = 6, main = paste(ICM, ' vs ', Var))
dev.off()

