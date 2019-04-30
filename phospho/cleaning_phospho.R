Phospho <- read.delim("~/Documents/Lund_Melanoma/phospho/PhosphoDIA_peptides_MM_Lund_Cohort_NORMALIZED.txt")
Phospho$label = paste(Phospho$EG.ModifiedSequence, Phospho$PG.ProteinAccessions, sep='~')
Phospho.cl = Phospho[,-c(1:7)]
rownames(Phospho.cl) = Phospho.cl$label
Phospho.cl = Phospho.cl[,c(1:122)]
Phospho.cly = Phospho.cl[-which(rowMeans(is.na(Phospho.cl)) > 0.5), ]
