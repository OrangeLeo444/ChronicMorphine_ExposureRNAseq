library(DESeq2)
library(dplyr)
library(limma)
library(readxl)
library(ggplot2)
library(forcats)
library(genefilter)
#Load data
cts <- read.table("data/AllSamplesCountsPerGene.txt", header=T, sep='\t')
annotation <- read_excel("data/elzar_morphine_RNAseq_sample_sheet.xlsx")
cData <- data.frame(annotation)
rownames(cts)<-cts$Geneid #only necessary if rows doesn't have a name
cts$Geneid<-NULL #only necessary if rows doesn't have a name

#give name to rows based on their sample_id and file_name_aligned
order_new <- factor(colnames(cts))
#cData<-cData[match(order_new,cData$file_name_aligned),] #for feature counts
cData<-cData[match(order_new,cData$sample_name),]
colnames(cts)<-cData$sample_id
rownames(cData)<-cData$sample_id

#filter data in case is necessary for a given condition use function "filter"

cData_filtered <- cData
cData_filtered<-filter(cData, brain_region == 'CeA' | brain_region == 'NAc')

#cData_filtered<-cData

#missing samples from count matrix 
noMissing<-rownames(cData_filtered)[rownames(cData_filtered) %in% colnames(cts)]

#excluding those samples 
cts_clean <- cts[,noMissing]
cData_filtered<-cData_filtered[colnames(cts_clean),]

#check the order
all(rownames(cData_filtered) %in% colnames(cts_clean))
all(rownames(cData_filtered)==colnames(cts_clean))

# strip transcripts with expression lower than avg 10 counts per sample
expressed_genes <- cts_clean[rowSums(cts_clean)>40,]

cData_filtered$brain_region<-factor(cData_filtered$brain_region)
cData_filtered$group<-factor(cData_filtered$group)


# create a DESeq data structure from raw counts for statistical testing
dds = DESeqDataSetFromMatrix(countData=expressed_genes, colData=cData_filtered, design= ~group)
# choose and assign reference samples
#dds$brain_region = relevel(dds$brain_region, ref="dStr")
dds$group = relevel(dds$group, ref="control")
# run DESeq to test data
dds = DESeq(dds)



# report results to res variable
#res = results(dds,alpha=0.05, contrast=c("brain_region","NAc","dStr"))
res = results(dds,alpha=0.05, contrast=c("group","morphine","control"))

summary(res)
sum(res$padj < 0.05, na.rm = TRUE)

#LFC shrinkage 
resLFC <- lfcShrink(dds, coef="brain_region_dStr_vs_CeA", type="apeglm")
DESeq2::plotMA(res,0.05,main='padj < 0.05 brain region',ylim=c(-7,7))
DESeq2::plotMA(resLFC,0.01,main='padj < 0.01 log fold change shrinkage brain region', ylim=c(-4,4))

# another normalization method with a variance stabilizing transform, a type of log transform
# good for comparing between datasets
vst_norm = vst(dds)


vst_corr_for_table <- assay(vst_norm)


plotPCA(vst_norm,"brain_region")
plotPCA(vst_norm,"group")



### Format res report
#order data based on logfold change 
res_no_format = res
res_no_format_log_fold_ordered = res_no_format[order(-res_no_format$log2FoldChange),]
#filter for significant and abundant genes
#get a boolean vector
res_of_interest_bool = !is.na(res_no_format_log_fold_ordered$padj) & res_no_format_log_fold_ordered$padj < 0.05 & (res_no_format_log_fold_ordered$log2FoldChange >= 1.5 | res_no_format_log_fold_ordered$log2FoldChange <= -1.5)
#filter on boolean vector
res_of_interest = res_no_format_log_fold_ordered[res_of_interest_bool,]
#filter on boolean vector
DESeq2::plotMA(res_of_interest,0.05,main='padj < 0.05 tissue+batch',ylim=c(-4,4))

###### Realize the batch correction 
#actual batch correction
assay(vst_norm) = limma::removeBatchEffect(assay(vst_norm), vst_norm$seq_batch)

# check to look at the corrected structure of the data

plotPCA(vst_norm,"group")




pcaData <- plotPCA(vst_norm, intgroup=c("passage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=factor(passage,levels=c('P1','P2','P3','P4','P7','P8','P9','P12')))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(color = "Passage") +scale_color_brewer(palette="Dark2")

library(RColorBrewer)
mycolors = c(brewer.pal(name="Set1", n = 9), brewer.pal(name="Pastel1", n = 5))
pcaData <- plotPCA(vst_norm, intgroup=c("patient_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=factor(patient_id,levels=c('R1N','R1T','R5N','R5T','R6N','R6T','R7T','R8T','R11N', 'R11T','R12N','R12T','R13T','R14T')))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(color = "Patient") +scale_color_manual(values=mycolors)


