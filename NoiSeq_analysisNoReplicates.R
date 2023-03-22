#Calculate DEG for one sample per group base on filtered genes
cts <- read.table("output/filtered_GeneCountMatrix_GENCODE.txt", header=T, sep='\t')
annotation <- read_excel("data/elzar_morphine_RNAseq_sample_sheet.xlsx")
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

cData_filtered <- filter(cData,brain_region == 'CeA')


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

old_annotations<-annotations
annotations <- annotations %>% distinct(ensembl_gene_id, .keep_all=TRUE)

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



myTMM <- tmm(assayData(myData)$exprs,long=1000,lc=0) #featureData(myData)$Length
myData_normalizedTMM<-readData(data=myTMM, length=mylength, gc=mygc, biotype=mybiotypes, chromosome=mychroms, factors=cData_filtered)

#mynoiseqSim <- noiseq(myData_normalizedTMM, k = 0.5, norm = "n", factor = "group", pnr = 0.2, nss = 5, v = 0.02,lc=0,replicates = "no")
mynoiseqSim <- noiseq(myData, k = NULL, norm = "tmm", factor = "group", pnr = 0.2, nss = 6, v = 0.02,lc=1,replicates = "no")
mynoiseqSim.deg <- degenes(mynoiseqSim, q=0.9, M=NULL)
#Plot expression with DEGs highlighted in red
DE.plot(mynoiseqSim, q = 0.9, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqSim, q = 0.9, graphic = "MD")



significant_genes_indx<-annotations[annotations$ensembl_gene_id %in% rownames(mynoiseqSim.deg),]
reorder_indx<-significant_genes_indx[match(rownames(mynoiseqSim.deg),significant_genes_id),]
all(rownames(mynoiseqSim.deg)==reorder_indx$ensembl_gene_id)
mynoiseqSim.deg$mgi_symbol<-reorder_indx$mgi_symbol
mynoiseqSim.deg.filtered_biotype<-mynoiseqSim.deg %>% filter(Biotype=='protein_coding')
CeA_mynoiseqsim.deg<-mynoiseqSim.deg.filtered_biotype
write.table(mynoiseqSim.deg, "output/CeA_mynoiseqsim_DGEA.txt",sep="\t",quote=F)

