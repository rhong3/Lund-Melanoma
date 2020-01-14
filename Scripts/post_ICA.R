## Post-ICA formating
#' 
#' Created on 1/14/2020
#' 
#' @author: RH
#' 

# proteomics
protICA <- read.delim("~/documents/Lund_Melanoma/Results/proteomics/ICA/ICA_proteomics_IC_centroid.txt")
protmix = read.delim("~/documents/Lund_Melanoma/Results/proteomics/ICA/ICA_proteomics_IC_mean_mixing_score.txt")
prot = read.csv("~/Documents/Lund_Melanoma/Data/proteomics/proteomics.csv")
prot = prot[,111:113]
rownames(protmix) = protmix$X
colnames(protICA) = c("Accession", rownames(protmix))
mgprotICA = merge(prot, protICA, by="Accession", all=FALSE, suffixes = c("",""))
write.csv(mgprotICA, "~/documents/Lund_Melanoma/Results/proteomics/ICA/MG_ICA_proteomics_IC_centroid.csv", row.names = FALSE)

# transcriptomics
transICA <- read.delim("~/documents/Lund_Melanoma/Results/transcriptomics/ICA/ICA_transcriptomics_IC_centroid.txt")
transmix = read.delim("~/documents/Lund_Melanoma/Results/transcriptomics/ICA/ICA_transcriptomics_IC_mean_mixing_score.txt")
trans = read.csv("~/Documents/Lund_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans = trans[,1:12]
rownames(transmix) = transmix$X
colnames(transICA) = c("Gene.name", rownames(transmix))
mgtransICA = merge(trans, transICA, by="Gene.name", all=FALSE, suffixes = c("",""))
write.csv(mgtransICA, "~/documents/Lund_Melanoma/Results/transcriptomics/ICA/MG_ICA_transcriptomics_IC_centroid.csv", row.names = FALSE)

# phospho
phosphoICA <- read.delim("~/documents/Lund_Melanoma/Results/phospho/ICA/ICA_phospho_IC_centroid.txt")
phosphomix = read.delim("~/documents/Lund_Melanoma/Results/phospho/ICA/ICA_phospho_IC_mean_mixing_score.txt")
phospho = read.csv("~/Documents/Lund_Melanoma/Data/phospho/phospho.csv")
phospho = phospho[,1:4]
rownames(phosphomix) = phosphomix$X
colnames(phosphoICA) = c("Modified_sequence", rownames(phosphomix))
mgphosphoICA = merge(phospho, phosphoICA, by="Modified_sequence", all=FALSE, suffixes = c("",""))
write.csv(mgphosphoICA, "~/documents/Lund_Melanoma/Results/phospho/ICA/MG_ICA_phospho_IC_centroid.csv", row.names = FALSE)

