library(ggrepel)

#Calculate DEG for one sample per group base on filtered genes
cts <- read.table("output/Filtered_GeneCountMatrix_GENCODE.txt", header=T, sep='\t')
annotation <- read_excel("data/elzar_morphine_RNAseq_sample_sheet.xlsx") 
annotations <- read.table("output/GeneCountMatrix_GENCODE_annotations.txt", header=T, sep='\t')
geneLegth <- read.table("output/GeneCountMatrix_GENCODE_genLength_GC.txt", header=T, sep='\t')

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

#filter data in case is necessary for a given condition use function "filter"

#cData_filtered <- filter(cData,brain_region == 'CeA')
cData_filtered <- filter(cData,(brain_region == 'NAc' | brain_region == 'dStr') & group == 'control' )


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

#old_annotations<-annotations
#annotations <- annotations %>% distinct(ensembl_gene_id, .keep_all=TRUE)

#Reorder annotations
annotations<-annotations[annotations$ensembl_gene_id %in% rownames(cts_clean),]
annotations<-annotations[match(rownames(cts_clean),annotations$ensembl_gene_id),]

mybiotypes <- data.frame(annotations$ensembl_gene_id,annotations$gene_biotype) 
mychroms <-  data.frame(annotations$ensembl_gene_id,annotations$chromosome_name,annotations$start_position,annotations$end_position)
rownames(mychroms)<-mychroms$annotations.ensembl_gene_id
mychroms$annotations.ensembl_gene_id<-NULL

##reorder geneLegth
my_gene_set<-geneLegth[rownames(geneLegth) %in% annotations$ensembl_gene_id,]
my_gene_set<-my_gene_set[match(annotations$ensembl_gene_id,rownames(my_gene_set)),]

mylength <- data.frame(rownames(my_gene_set),my_gene_set$length)

mygc <- data.frame(rownames(my_gene_set),my_gene_set$gc)

colnames(cts_clean)<-cData_filtered$sample_condition
rownames(cData_filtered)<-cData_filtered$sample_condition

myData<-readData(data=cts_clean, length=mylength, gc=mygc, biotype=mybiotypes, chromosome=mychroms, factors=cData_filtered)

##NORMALIZE AND RUN NOISEQ
#myfilt <- filtered.data(cts_clean, factor = myData$group, norm=FALSE, depth=NULL, method=1, cv.cutoff=100, cpm=1, p.adj="fdr")




#mynoiseqSim <- noiseq(myData_normalizedTMM, k = 0.5, norm = "n", factor = "group", pnr = 0.2, nss = 5, v = 0.02,lc=0,replicates = "no")
mynoiseqSim <- noiseq(myData, k = NULL, norm = "tmm", factor = "brain_region", pnr = 0.2, nss = 6, v = 0.02,lc=1,replicates = "no")
mynoiseqSim.deg <- degenes(mynoiseqSim, q=0.9, M=NULL)
#Plot expression with DEGs highlighted in red
DE.plot(mynoiseqSim, q = 0.9, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqSim, q = 0.9, graphic = "MD")



significant_genes_indx<-annotations[annotations$ensembl_gene_id %in% rownames(mynoiseqSim.deg),]
#reorder_indx<-significant_genes_indx[match(rownames(mynoiseqSim.deg),significant_genes_id),]
reorder_indx<-significant_genes_indx[match(rownames(mynoiseqSim.deg),significant_genes_indx$ensembl_gene_id),]
all(rownames(mynoiseqSim.deg)==reorder_indx$ensembl_gene_id)
mynoiseqSim.deg$mgi_symbol<-reorder_indx$mgi_symbol
mynoiseqSim.deg.filtered_biotype<-mynoiseqSim.deg %>% filter(Biotype=='protein_coding')
#CeA_mynoiseqsim.deg<-mynoiseqSim.deg.filtered_biotype
#write.table(mynoiseqSim.deg, "output/CeA_mynoiseqsim_DGEA.txt",sep="\t",quote=F)
#CeAvsNAc_mynoiseqsim.deg<-mynoiseqSim.deg.filtered_biotype
#write.table(mynoiseqSim.deg, "output/CeAvsNAc_mynoiseqsim_DGEA.txt",sep="\t",quote=F)
#CeAvsdStr_mynoiseqsim.deg<-mynoiseqSim.deg.filtered_biotype
#write.table(mynoiseqSim.deg, "output/CeAvsdStr_mynoiseqsim_DGEA.txt",sep="\t",quote=F)


NAcvsdStr_mynoiseqsim.deg<-mynoiseqSim.deg.filtered_biotype
write.table(mynoiseqSim.deg, "output/NAcvsdStr_mynoiseqsim_DGEA.txt",sep="\t",quote=F)
saveRDS(mynoiseqSim,"output/NAcvsdStr_mynoiseqsim_DGEA.rds")

data_q_09<-mynoiseqSim.deg.filtered_biotype
##### new volcano plot######
logFold_threshold<-1
p_value_threshold<-0.05
biotype<-'protein_coding'
data_Exploratory<-data.frame(pvalue=1-mynoiseqSim@results[[1]]$prob,log2FoldChange=mynoiseqSim@results[[1]]$M)
data_Exploratory$D<-mynoiseqSim@results[[1]]$D
rownames(data_Exploratory)<-rownames(mynoiseqSim@results[[1]])
id_gene<-rownames(data_Exploratory)
filtered_annotations<-annotations %>% filter (ensembl_gene_id %in% id_gene)
filtered_annotations<-filtered_annotations[match(id_gene,filtered_annotations$ensembl_gene_id),]
data_Exploratory$gene_symbol<-filtered_annotations$mgi_symbol
data_Exploratory$biotype<-filtered_annotations$gene_biotype
data_Exploratory$chroms<-filtered_annotations$chromosome_name

p <- ggplot(data=data_Exploratory, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
# add a column of NAs
data_Exploratory$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
data_Exploratory$diffexpressed[data_Exploratory$log2FoldChange > logFold_threshold & data_Exploratory$pvalue < p_value_threshold] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data_Exploratory$diffexpressed[data_Exploratory$log2FoldChange < -(logFold_threshold) & data_Exploratory$pvalue < p_value_threshold] <- "DOWN"
p <- ggplot(data=data_Exploratory, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-logFold_threshold, logFold_threshold), col="red") +
  geom_hline(yintercept=-log10(p_value_threshold), col="red")
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
data_Exploratory$delabel <- NA
#data_Exploratory$delabel[data_Exploratory$gene_symbol %in% data_q_09$mgi_symbol] <- data_Exploratory$gene_symbol[data_Exploratory$gene_symbol %in% data_q_09$mgi_symbol]
data_Exploratory$delabel[data_Exploratory$diffexpressed != "NO" & data_Exploratory$biotype=="protein_coding"] <- data_Exploratory$gene_symbol[data_Exploratory$diffexpressed != "NO" & data_Exploratory$biotype=="protein_coding" & data_Exploratory$chroms != "MT"]
data_Exploratory$delabel[data_Exploratory$diffexpressed != "NO" & data_Exploratory$biotype=="protein_coding"] <- data_Exploratory$gene_symbol[data_Exploratory$diffexpressed != "NO" & data_Exploratory$biotype=="protein_coding"]

ggplot(data=data_Exploratory, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# volcano with names in nice way
ggplot(data=data_Exploratory, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-logFold_threshold, logFold_threshold), col="red") +
  geom_hline(yintercept=-log10(p_value_threshold), col="red")

ggplot(data=data_Exploratory, aes(x=log2FoldChange, y=D, col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) + 
  scale_y_log10(limits = c(1, 10000), breaks = c(1, 10, 100, 1000,10000))

ggplot(data=data_Exploratory, aes(x=log2FoldChange, y=D, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) + 
  scale_y_log10(limits = c(1, 10000), breaks = c(1, 10, 100, 1000,10000))