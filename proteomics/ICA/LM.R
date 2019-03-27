wna_clinical <- read.delim("~/Documents/Lund/proteomics/wna_clinical.tsv")
Clean_proteomics_IC_mean_mixing_score <- read.delim("~/Documents/Lund/proteomics/ICA/Initial/Clean_proteomics_IC_mean_mixing_score.txt")
survival1 = wna_clinical[,10]
survival2 = wna_clinical[,11]
IC49 = as.matrix(Clean_proteomics_IC_mean_mixing_score[49,2:137])
IC60 = as.matrix(Clean_proteomics_IC_mean_mixing_score[60,2:137])
IC62 = as.matrix(Clean_proteomics_IC_mean_mixing_score[62,2:137])
plot(IC62[1,-c(117)],survival1[-c(117)],xlim=c(-0.1,0.2), xlab="IC62", ylab="5 year survival from primary diagnosis")
title(main = 'IC62 vs 5 year survival from primary diagnosis ')
plot(IC62[1,-c(117)],survival2[-c(117)],xlim=c(-0.1,0.2), xlab="IC62", ylab="5 year survival from first metastasis")
title(main = 'IC62 vs 5 year survival from first metastasis ')
x1=lm(formula = survival1[-c(117)]~IC49[1,-c(117)]+IC60[1,-c(117)]+IC62[1,-c(117)], na.action=na.omit)
x2=lm(formula = survival2[-c(117)]~IC49[1,-c(117)]+IC60[1,-c(117)]+IC62[1,-c(117)], na.action=na.omit)
summary(x1)

summary(x2)



