library(DESeq2)
cts <- read.table("output/Filtered_GeneCountMatrix_GENCODE.txt", header=T, sep='\t')
annotation <- read_excel("data/elzar_morphine_RNAseq_sample_sheet.xlsx") 
#annotations <- read.table("output/GeneCountMatrix_GENCODE_annotations.txt", header=T, sep='\t')
annotations <- old_annotations

cData <- data.frame(annotation)
cts <- cts %>% rownames_to_column(var="gene") %>% separate(gene, into=c("gene", "version"), sep = "\\.")
cts$version<-NULL
rownames(cts)<-cts$gene
cts$gene<-NULL


### filter brain region 
order_new <- factor(colnames(cts))
cData<-cData[match(order_new,cData$sample_condition ),]
colnames(cts)<-cData$sample_condition
rownames(cData)<-cData$sample_condition

#cData_filtered <- filter(cData,brain_region == 'dStr')
cData_filtered <- cData

#cData_filtered<-cData

#missing samples from count matrix 
noMissing<-rownames(cData_filtered)[rownames(cData_filtered) %in% colnames(cts)]

#excluding those samples 
cts_clean <- cts[,noMissing]
cData_filtered<-cData_filtered[colnames(cts_clean),]

#check the order
all(rownames(cData_filtered) %in% colnames(cts_clean))
all(rownames(cData_filtered)==colnames(cts_clean))
cts_clean<-cts_clean[rowSums(cts_clean)>20,]



#Reorder annotations
annotations<-annotations[annotations$ensembl_gene_id %in% rownames(cts_clean),]
annotations<-annotations[match(rownames(cts_clean),annotations$ensembl_gene_id),]
id_symbol<-tibble(ensembl_gene_id=annotations$ensembl_gene_id,mgi_symbol=annotations$mgi_symbol)

anno <- AnnotationDbi::select(org.Mm.eg.db, 
                             keys=rownames(cts_clean), 
                             columns=c("SYMBOL", "GENENAME"),
                             keytype="ENSEMBL")
anno<-anno %>% distinct(ENSEMBL, .keep_all = TRUE)
colnames(anno)<-c('ensembl_gene_id','SYMBOL','GENENAME')


cts_clean$ensembl_gene_id<-rownames(cts_clean)
cts_clean_merge<-left_join(cts_clean,id_symbol, by='ensembl_gene_id')
cts_clean_merge<-left_join(cts_clean_merge,anno, by='ensembl_gene_id')
cts_clean<-cts_clean_merge[!is.na(cts_clean_merge$SYMBOL),]
rownames(cts_clean)<-cts_clean$mgi_symbol


#cts_clean[rownames(cts_clean)=='ENSMUSG00000079737',] 
cts_clean<-cts_clean[!rownames(cts_clean)=='ENSMUSG00000079737',] 
rownames(cts_clean)<-cts_clean$mgi_symbol
#write.table(cts_clean, "output/cts_clean_annotated.txt",sep="\t",quote=F)
cts_clean$mgi_symbol<-NULL
cts_clean$ensembl_gene_id<-NULL
cts_clean$GENENAME<-NULL
cts_clean$SYMBOL<-NULL


expressed_genes<-cts_clean
# create a DESeq data structure from raw counts for statistical testing
cData_filtered$brain_region<-factor(cData_filtered$brain_region)
cData_filtered$group<-factor(cData_filtered$group)
dds = DESeqDataSetFromMatrix(countData=expressed_genes, colData=cData_filtered, design= ~group)
# choose and assign reference samples
#dds$brain_region = relevel(dds$brain_region, ref="dStr")
dds$group = relevel(dds$group, ref="control")
# run DESeq to test data
dds = DESeq(dds)

test<-plotCounts(dds,gene = "Drd1",intgroup = "brain_region", returnData = TRUE)
normalized_cts<-counts(dds, normalized=TRUE)
#plot_genes<-normalized_cts[rownames(normalized_cts) %in% c('Drd1','Drd2','Snc4b','Prkcd','Peg10','Tac2','Camk1g','Foxp1','Calb1','Oprm1','Nr4a1'),]
plot_genes<-normalized_cts[rownames(normalized_cts) %in% c('Drd1','Dlk1','Pdyn','Drd2','Snc4b','Prkcd','Nnat','Peg10','Tac2','Camk1g','Calb1','Oprm1','Oprk1','Nr4a1','Gabrg1','Wfs1','Acvrl1','Cartpt'),]


library(reshape2) 
#install.packages("viridis")  # Install
library("viridis")    
melted_plot_genes<-melt(plot_genes)

melted_plot_genes<-melted_plot_genes %>%  separate(Var2, into=c("brain_region", "condition"), sep = "\\_")

df.summary2 <- melted_plot_genes %>%
  group_by(brain_region, Var1) %>%
  summarise(
    sd = sd(value),
    counts.mean.normalized = mean(value)
  )

ggbarplot(df.summary2, x = "Var1", y = "counts.mean.normalized", 
          color = "brain_region", palette = c("#E11F26", "#387EB8",'#4DAE49'),
          position = position_dodge(0.8)) +   xlab("Genes") +ylab("Normalized Counts") + grids()
  
ggplot(df.summary2) +
  geom_bar(aes(x = Var1, y = counts.mean.normalized, color = brain_region )) +
  scale_y_log10(limits = c(1, 1000), breaks = c(1, 10, 100, 1000)) +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Expression") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
