#' Proteomics Clean-up
#' 
#' Created on 1/31/2019
#' 
#' @author: RH


library(readxl)
Original <- read_excel("Segundo TMT data to Krzysztof.xlsx")
library(tibble)
library(data.table)
colnames(Original) = Original[3,]
Original.cl1 = Original[-c(1:9),-1]
Original.cl1 = Original.cl1[,1:151]
colnames(Original.cl1)[151] = 'Accession'
Original.cl1 = data.frame(column_to_rownames(Original.cl1, var = "Accession"))
Original.clt = transpose(Original.cl1)
colnames(Original.clt) = row.names(Original.cl1)
rownames(Original.clt) = colnames(Original.cl1)
remove = c('Fake_sample', 'Fake_sample.1', 'Fake_sample.2', 'Fake_sample.3', 'REFERENCE2.3', 'REFERENCE2.2', 'REFERENCE2.1', 'REFERENCE2', 'MM1069')
Original.cltf = subset(Original.clt, !(rownames(Original.clt) %in% remove))
write.csv(Original.cltf, "rm_Clean1.csv")

library(BBmisc)
library(dplyr)
samp = rownames(Original.cltf)
Original.cltf <- mutate_all(Original.cltf, function(x) as.numeric(as.character(x)))
sapply(Original.cltf, class)
rownames(Original.cltf) = samp
write.csv(Original.cltf, "rm_Clean1_5.csv")

Original.norm = normalize(Original.cltf, method = "standardize", margin = 1)
write.csv(Original.norm, "rm_Clean2.csv")

library(reshape2)
library(ggplot2)

Original.normt = transpose(Original.norm)
colnames(Original.normt) = row.names(Original.norm)
rownames(Original.normt) = colnames(Original.norm)
Original.normt$accession = rownames(Original.normt)
long = melt(Original.normt, id.vars = "accession")
p <- ggplot(aes(x=value, colour=variable), data=long)
p + geom_density()
summary(Original.normt)
write.csv(Original.normt, "rm_Clean3.csv")
ggsave('rm_Density.png', scale = 1, width = 30, height = 12, units = "in", limitsize = TRUE)


