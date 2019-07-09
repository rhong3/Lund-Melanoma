ICs = list(list('8', 'class'), list('40', 'class'), list('45', 'survival'), list('85', 'survival'), list('102', 'survival'), list('105', 'stage'), list('105', 'local'), list('107','survival'))
TICs = list(list('38', 'class'), list('96', 'class'), list('96', 'BRAF'))
# Aggregate significant GSEA pathways
finalfile <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c('pathway', 'pval',	'padj',	'ES',	'NES',	'nMoreExtreme',	'size',	'leadingEdge', 'cliniccal'))
for(i in TICs){
  nname = paste('IC', i[1], sep='')
  nname = paste(nname, '.csv', sep='')
  nname = paste("~/Documents/Lund_Melanoma/Transcriptome/ICA/0405ICA/GSEA/", nname, sep='')
  IC.file = read.csv(nname, row.names=1)
  IC.file = IC.file[IC.file$padj < 0.01, ]
  modess = i[2]
  IC.file$clinical = modess
  finalfile <- rbind(finalfile, IC.file)
}
finalfile$clinical = as.character(finalfile$clinical)
finalfile = finalfile[order(finalfile$padj), ]
write.csv(finalfile, file = "~/Documents/Lund_Melanoma/Transcriptome/ICA/0405ICA/GSEA/Aggregated.csv", row.names=TRUE)

# Aggregate significant genes
genedict = list(list('5YR_Binary', 'T-Test'), list('binary', 'T-Test'), list('BRAF', 'T-Test'), list('clark', 'ANOVA'), list('clin_class', 'ANOVA'), list('dis_stage', 'ANOVA'), list('prim_breslow_class', 'ANOVA'), list('prim_site', 'ANOVA'))
finalfile <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c('Gene', 'Pvalue', 'Test', 'Clinical'))
for(i in genedict){
  nname = paste('~/Documents/Lund_Melanoma/phospho/Analysis0709/', i[1], '/', i[2], '.csv', sep='')
  filess = read.csv(nname)
  colnames(filess) = c('Gene', 'Pvalue')
  filess = filess[filess$Pvalue < 0.01, ]
  filess$Test = i[2]
  filess$Clinical = i[1]
  finalfile = rbind(finalfile, filess)
}
finalfile$Test = as.character(finalfile$Test)
finalfile$Clinical = as.character(finalfile$Clinical)
finalfile = finalfile[order(finalfile$Pvalue), ]
write.csv(finalfile, file = "~/Documents/Lund_Melanoma/phospho/Analysis0709/Gene_Aggregated.csv", row.names=TRUE)


