# Match back to gene name
phospho <- read.csv("~/Documents/Lund_Melanoma/phospho/Clean_70_phospho.csv")
Clean_phospho_IC_centroid <- read.delim("~/Documents/Lund_Melanoma/phospho/ICA/Clean_phospho_IC_centroid.txt")

rownames(phospho) = phospho$X
phospho = phospho[,c(1,2)]
rownames(Clean_phospho_IC_centroid) = Clean_phospho_IC_centroid$X
Clean_phospho_IC_centroid = Clean_phospho_IC_centroid[,-c(1)]
new = merge(Clean_phospho_IC_centroid, phospho)
new.w = aggregate(new, by=list(new$Gene), FUN = mean)
rownames(new.w) = new.w$Group.1
new.w = new.w[,c(3:120)]
write.csv(new.w, "~/Documents/Lund_Melanoma/phospho/ICA/Gene_phospho_ICA_Centroid.csv")
