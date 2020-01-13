## Data clearning
#' 
#' Created on 1/13/2020
#' 
#' @author: RH
library(readxl)
library(readr)
library(plyr)
library(data.table)
library(dplyr)
# Clinical and histology
histology = read_excel("~/Documents/Lund_Melanoma/Data/Sumary Histology evaluation, Segundo study.xlsx")
histology[1,9] = (histology[1,9] + histology[113,9])/2
histology[109,1] = "MM808"
histology = histology[histology$`Case ID` != "MM814",]
histology[1,1] = "MM814"
rownames(histology) = histology$`Case ID`

clinical = read_excel("~/Documents/Lund_Melanoma/Data/20190130ClinicalData_144samplesLundMM_Query2018-12HE-10.xlsx")
clinical = clinical[-1,-c(12,21,22,24,25,32:42)]
clinical$sample = toupper(clinical$sample)
rownames(clinical) = clinical$sample
clinical$stage <- mapvalues(clinical$stage, from=c("Local", "In transit", "Regional", "General"), to=c(1, 2, 3, 4))
clinical["dist.met.location.Sc/Lgl"] = as.numeric(clinical$dist.met.location %like% "0")
clinical["dist.met.location.Liver"] = as.numeric(clinical$dist.met.location %like% "1")
clinical["dist.met.location.Lung"] = as.numeric(clinical$dist.met.location %like% "2")
clinical["dist.met.location.Bone"] = as.numeric(clinical$dist.met.location %like% "3")
clinical["dist.met.location.CNS"] = as.numeric(clinical$dist.met.location %like% "4")
clinical["dist.met.location.OthersViceral"] = as.numeric(clinical$dist.met.location %like% "5")

clinical["Local.Cutaneous"] = as.numeric(clinical$local == "Cutaneous")
clinical["Local.Lymph.node"] = as.numeric(clinical$local == "Lymph node")
clinical["Local.Subcutaneous"] = as.numeric(clinical$local == "Subcutaneous")
clinical["Local.Visceral"] = as.numeric(clinical$local == "Visceral")

clinical$`Type of first metastasis` = mapvalues(clinical$`Type of first metastasis`, from=c("Local", "In-transit", "Regional", 
                                                                                            "General", "In-transit, regional, distant", 
                                                                                            "In transit, Regional", "Regiona, distant", 
                                                                                            "Regional, Distant", "Satellite", "In-transit, Regional", 
                                                                                            "Regional, distant", "In-transit, regional"), to=c(1, 2, 3, 4, 4, 3, 4, 4, 3, 3, 4, 3))

clinical$`No. Involved lymph nodes (first met)` = mapvalues(clinical$`No. Involved lymph nodes (first met)`, from=c("unknown", "1", "2-3", ">3"), to=c("NA",1,2,3))
clinical$gender = mapvalues(clinical$gender, from=c('Male', "Female"), to=c(0,1))
clinical = clinical[,-which(colnames(clinical) %in% c("dist.met.location", "Type of first metastasis", "No. Involved lymph nodes (first met)", "local"))]
rownames(clinical) = clinical$sample

his_clin = merge(clinical, histology, by=0, all=TRUE, suffixes = c("",""))
rownames(his_clin) = his_clin$sample
his_clin = his_clin[,-which(colnames(his_clin) %in% c("Row.names", "Case ID", "sample"))]
write.csv(his_clin, "~/Documents/Lund_Melanoma/Data/clinial&histology.csv")


# Proteomics
prot <- read_delim("~/Documents/Lund_Melanoma/Data/proteomics/20191211_Segundo_MM_TMT11_isotope_correction_normalized_ratios_Final.txt", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 3)
prot['MM692'] = (prot['MM692']+prot['MM692_1'])/2
prot['MM778'] = (prot['MM778']+prot['MM778_1'])/2
prot['MM807'] = (prot['MM807']+prot['MM807_1'])/2
prot['MM790'] = (prot['MM790']*77.76+prot['MM790_1']*78.1923076923077)/(77.76+78.1923076923077)*2
prot['MM814'] = (prot['LG7438']*83.4+prot['MM814']*96.2)/(96.2+83.4)*2
tumor <- read_excel("~/Documents/Lund_Melanoma/Data/Sample tumor content.xlsx")
tumor.rm = tumor[tumor$`% Tumor content`<30,]$TumorID
remove = c('Fake_sample', 'Fake_sample_1', 'Fake_sample_2', 'Fake_sample_3', 'REFERENCE2_3', 'REFERENCE2_2', 'REFERENCE2_1', 'REFERENCE2', 'MM1069', 'MM692_1', 'MM778_1', 'MM807_1', 'MM790_1', 'LG7438')
prot.sample = subset(prot, select = !names(prot) %in% remove)
prot.sample = subset(prot.sample, select = !names(prot.sample) %in% tumor.rm)
prot_his_clin <- read.csv("~/Documents/Lund_Melanoma/Data/clinial&histology.csv", row.names=1)
prot.sample = subset(prot.sample, select = names(prot.sample) %in% c(intersect(colnames(prot.sample), rownames(prot_his_clin)), "Accession", "Gene name", "Description"))
prot.sample = na.omit(prot.sample)
write.csv(prot.sample, "~/Documents/Lund_Melanoma/Data/proteomics/proteomics.csv", row.names = FALSE)
prot.sample.ICA = select(prot.sample, -c("Gene name", "Description", "Accession"))
rownames(prot.sample.ICA) = prot.sample$Accession
write.csv(prot.sample.ICA, "~/Documents/Lund_Melanoma/Data/proteomics/ICA_proteomics.csv", row.names = TRUE)
prot_his_clin = prot_his_clin[intersect(colnames(prot.sample), rownames(prot_his_clin)), ]
write.csv(prot_his_clin, "~/Documents/Lund_Melanoma/Data/proteomics/proteomics_clinial.csv")


# Transcriptomics
Transcriptomics = read_excel("~/Documents/Lund_Melanoma/Data/transcriptomics/Transcriptomics_final.214.allMm.not.meancentered_rounded.xlsx")
colnames(Transcriptomics) = toupper(colnames(Transcriptomics))
Transcriptomics = Transcriptomics[,-c(1:4,6,8,9)]
colnames(Transcriptomics) = c("Gene name", "Description", c(colnames(Transcriptomics)[3:216]))
Trans_his_clin <- read.csv("~/Documents/Lund_Melanoma/Data/clinial&histology.csv", row.names=1)
Trans.sample = subset(Transcriptomics , select = names(Transcriptomics ) %in% c(intersect(colnames(Transcriptomics ), rownames(Trans_his_clin)), "Gene name", "Description"))
Trans.sample = na.omit(Trans.sample)
write.csv(Trans.sample, "~/Documents/Lund_Melanoma/Data/transcriptomics/transcriptomics.csv", row.names = FALSE)
Trans.sample.ICA = select(Trans.sample, -c("Gene name", "Description"))
rownames(Trans.sample.ICA) = Trans.sample$`Gene name`
write.csv(Trans.sample.ICA, "~/Documents/Lund_Melanoma/Data/transcriptomics/ICA_transcriptomics.csv", row.names = TRUE)
Trans_his_clin = Trans_his_clin[intersect(colnames(Trans.sample), rownames(Trans_his_clin)), ]
write.csv(Trans_his_clin, "~/Documents/Lund_Melanoma/Data/transcriptomics/transcriptomics_clinial.csv")


# Phospho
# filepath is path to 'PhosphoDIA_peptides_MM_Lund_Cohort_NORMALIZED.txt'

readAndImputePhosphoData <- function(filepath) {
  library(RANN)
  inputdata <- read.delim(filepath)
  phosphoExpression <- inputdata[, 8:ncol(inputdata)]
  phosphoExpression <- apply(phosphoExpression, 2, as.numeric)
  colnames(phosphoExpression) <- colnames(inputdata[8:ncol(inputdata)])
  rownames(phosphoExpression) <- inputdata$EG.ModifiedSequence
  rownames(inputdata) = inputdata$EG.ModifiedSequence
  # Filter out features found in too few samples, using 95% presence here
  rMiss <- function(x) { sum(is.na(x)) / length(x) }
  phosMiss <- apply(phosphoExpression, 1, rMiss)
  
  # KNN impute with caret (alternatively, use the knnImpute function from the PAMR package?)
  x <- t(phosphoExpression[phosMiss < 0.05, ])
  set.seed(1011)
  library(caret)
  x.imputeModel = preProcess(x, "knnImpute", k = 2)
  
  phosExpr = predict(x.imputeModel, x)
  phosExpr = data.frame(t(phosExpr))
  
  phosExpr = merge(inputdata[, 1:7], phosExpr, by=0, all=FALSE, suffixes = c("",""))
  
  return(phosExpr)
  
}

Phospho = readAndImputePhosphoData("~/Documents/Lund_Melanoma/Data/phospho/MM_Phosho_DIA_lund_Cohort_total_peptides_NORMALIZED.txt")
Phospho = Phospho[,-c(1,3,6,7,9)]
Phospho$MM790 = rowMeans(Phospho[c('MM790.1.', 'MM790.2.')], na.rm=TRUE)
Phospho$MM807 = rowMeans(Phospho[c('MM807.1.', 'MM807.2.')], na.rm=TRUE)
Phospho$MM808 = Phospho$MM808_LG
Phospho = select(Phospho, -c(MM790.1., MM790.2., MM807.1., MM807.2., MM808_LG))
colnames(Phospho) = c("Modified_sequence", "Accession", "Description", "Gene name", colnames(Phospho)[5:124])
Phospho_his_clin <- read.csv("~/Documents/Lund_Melanoma/Data/clinial&histology.csv", row.names=1)
Phospho.sample = subset(Phospho , select = names(Phospho) %in% c(intersect(colnames(Phospho), rownames(Phospho_his_clin)), "Modified_sequence", "Accession", "Description", "Gene name"))
Phospho.sample = na.omit(Phospho.sample)
write.csv(Phospho.sample, "~/Documents/Lund_Melanoma/Data/phospho/phospho.csv", row.names = FALSE)
Phospho.sample.ICA = select(Phospho.sample, -c("Modified_sequence", "Accession", "Description", "Gene name"))
rownames(Phospho.sample.ICA) = Phospho.sample$`Modified_sequence`
write.csv(Phospho.sample.ICA, "~/Documents/Lund_Melanoma/Data/phospho/ICA_phospho.csv", row.names = TRUE)
Phospho_his_clin = Phospho_his_clin[intersect(colnames(Phospho.sample), rownames(Phospho_his_clin)), ]
write.csv(Phospho_his_clin, "~/Documents/Lund_Melanoma/Data/phospho/phospho_clinial.csv")

  
  
  
  
  
  