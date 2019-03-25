# Heatmap
ica <- read.delim("~/Documents/Lund/proteomics/ICA/Gene_level/J_Clean_proteomics_IC_centroid.txt", row.names=1)
wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv", row.names=1)
proteomics <- read.delim("~/Documents/Lund/proteomics/J_Clean_proteomics.tsv", row.names=1)
ica$survival = ica$X92
ica = ica[order(ica$survival),]
wna_clinicalC = wna_clinical[order(wna_clinical['clin.class.det_ALM'],	wna_clinical['clin.class.det_LMM'],	wna_clinical['clin.class.det_Mucosal'],	wna_clinical['clin.class.det_NM'],	wna_clinical['clin.class.det_Other'],	wna_clinical['clin.class.det_SSM'],	wna_clinical['clin.class.det_Unknownprimary']),]
wna_clinicalS = wna_clinical[order(wna_clinical['Alive.2016.12.05_alive'],	wna_clinical['Alive.2016.12.05_dead'],	wna_clinical['Alive.2016.12.05_dead..likely.melanoma.'],	wna_clinical['Alive.2016.12.05_dead.other.reason'],	wna_clinical['Alive.2016.12.05_dead.unknown.reason']),]

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

GSEA_proteomics = data.matrix(GSEA_proteomics[-c(26:7984),])
GSEA_proteomics = subset(GSEA_proteomics, nchar(as.character(rownames(GSEA_proteomics))) <= 10)
pdf("~/Documents/Lund/proteomics/ICA/Gene_level/GSEA/Cluster_IC92survival_HM.pdf", width = 15, paper = 'a4r')
pheatmap(GSEA_proteomics, cluster_cols = T, cluster_rows = F, annotation_col = categoryS, fontsize_col = 6, main = 'IC92 vs 5yr Survival')
dev.off()
#########################################

J_Clean_proteomics_IC_centroid <- read.delim("~/Documents/Lund/proteomics/ICA/Gene_level/J_Clean_proteomics_IC_centroid.txt", row.names=1)[17:11051,]
Gene_order = J_Clean_proteomics_IC_centroid[order(J_Clean_proteomics_IC_centroid$X96),]
library('org.Hs.eg.db')
ENTREZID = mapIds(org.Hs.eg.db, row.names(Gene_order), 'ENTREZID', 'SYMBOL')
Gene_order.96 <- setNames(as.numeric(Gene_order$X96), unname(ENTREZID))
# source("https://bioconductor.org/biocLite.R")
# biocLite("reactome.db")
my_pathways <- reactomePathways(names(Gene_order.96))
summary(sapply(my_pathways, length))
fgsea_reactome <- fgsea(pathways = my_pathways, 
                        stats = Gene_order.96,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)
fgsea_reactome <- na.omit(fgsea_reactome[order(pval), ])
fgsea_reactome$leadingEdge = as.character(fgsea_reactome$leadingEdge)
fgsea_reactome.sig = fgsea_reactome[fgsea_reactome$padj < 0.01,]
write.csv(fgsea_reactome, file = "~/Documents/Lund/proteomics/ICA/Gene_level/GSEA/IC96.csv")
# pdf("~/Documents/Lund/proteomics/ICA/Gene_level/GSEA/IC96_most_enriched.pdf", paper = 'letter')
# plotEnrichment(my_pathways[[head(fgsea_reactome, 1)$pathway]], Gene_order.96) + labs(title=head(fgsea_reactome, 1)$pathway)
# dev.off()

topPathwaysUp <- fgsea_reactome.sig[ES > 0][head(order(padj), n=20), pathway]
topPathwaysDown <- fgsea_reactome.sig[ES < 0][head(order(padj), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("~/Documents/Lund/proteomics/ICA/Gene_level/GSEA/IC96_sig.pdf", width = 15, paper = 'a4r')
plotGseaTable(my_pathways[topPathways], Gene_order.96, fgsea_reactome, 
              gseaParam = 0.5)
dev.off()

###################################################
