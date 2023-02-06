library(DESeq2)
library(dplyr)
library(limma)
library(readxl)
library(ggplot2)
library(forcats)
library(genefilter)
#Load data
cts <- read.table("data/AllSamples_featureCounts_output.Rmatrix.txt", header=T, sep='\t')
annotation <- read_excel("data/elzar_morphine_RNAseq_sample_sheet.xlsx")
cData <- data.frame(annotation)
rownames(cts)<-cts$Geneid
cts$Geneid<-NULL

#give name to rows based on their facility_id and patient_id
new_id_mouse<-sub("N","C",cData$patient_id)