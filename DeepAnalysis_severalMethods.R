library(DESeq2)
library(dplyr)
library(limma)
library(readxl)
library(ggplot2)
library(forcats)
library(genefilter)
library(tidyverse)
library(edgeR)
library(limma)
library(ggpubr)
library(factoextra)
library(biomaRt)
library(EDASeq)


#Loading data
cts <- read.table("output/GeneCountMatrix_GENCODE.txt", header=T, sep='\t')
annotation <- read_excel("data/elzar_morphine_RNAseq_sample_sheet.xlsx")
cData <- data.frame(annotation)
cts <- cts %>% rownames_to_column(var="gene") %>% separate(gene, into=c("gene", "version"), sep = "\\.")


#Getting annotation
ensembl <-  useDataset("mmusculus_gene_ensembl",mart=useMart("ensembl"))
annotations <- getBM(attributes=c("ensembl_gene_id","mgi_symbol", "description", "gene_biotype","chromosome_name","start_position","end_position"), 
                     filters = "ensembl_gene_id", 
                     values=cts$gene, 
                     mart=ensembl)

mybiotypes <- data.frame(annotations$ensembl_gene_id,annotations$gene_biotype) 
mychroms <-  data.frame(annotations$ensembl_gene_id,annotations$chromosome_name,annotations$start_position,annotations$end_position)
my_gene_set<-rownames(geneLegth)
mylength <- data.frame(my_gene_set,geneLegth$length)
mygc <- data.frame(my_gene_set,geneLegth$gc)

myData<-readData(data=cts, length=mylength, gc=mygc, biotype=mybiotypes, chromosome=mychroms, factors=cData)
#geneLegth<-getGeneLengthAndGCContent(cts$gene,"mmusculus_gene_ensembl")

