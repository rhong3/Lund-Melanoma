# Match back to gene name
Phospho <- read.csv("~/Documents/Lund_Melanoma/phospho/Clean_70_phospho.csv")
Phospho = Phospho[, c(1:8)]
Clean_phospho_IC_centroid <- read.delim("~/Documents/Lund_Melanoma/phospho/Clean_phospho_ip.tsv")

rownames(Phospho) = Phospho$X
Phospho = Phospho[,c(1,2)]
rownames(Clean_phospho_IC_centroid) = Clean_phospho_IC_centroid$X
Clean_phospho_IC_centroid = Clean_phospho_IC_centroid[,-c(1)]
new = merge(Clean_phospho_IC_centroid, Phospho, by="row.names")
new.w = aggregate(new, by=list(new$Gene), FUN = mean)
rownames(new.w) = new.w$Group.1
new.w = new.w[,c(3:120)]
write.table(new.w, file="~/Documents/Lund_Melanoma/phospho/J_Clean_phospho_ip.tsv", quote=FALSE, sep='\t')
