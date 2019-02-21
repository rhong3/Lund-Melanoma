#' transcriptomics Clean-up
#' 
#' Created on 2/20/2019
#' 
#' @author: RH
setwd("~documents/proteomics/Transcriptome")
Clean_transcriptome <- read.delim("~/Documents/proteomics/Transcriptome/Clean_transcriptome.tsv", row.names=1)
library(tibble)
library(data.table)
Original.cltf = Clean_transcriptome

library(BBmisc)
library(dplyr)
samp = rownames(Original.cltf)
Original.cltf <- mutate_all(Original.cltf, function(x) as.numeric(as.character(x)))
sapply(Original.cltf, class)
rownames(Original.cltf) = samp

Original.norm = normalize(Original.cltf, method = "standardize", margin = 1)

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
ggsave('rm_Density.png', scale = 1, width = 30, height = 12, units = "in", limitsize = TRUE)


