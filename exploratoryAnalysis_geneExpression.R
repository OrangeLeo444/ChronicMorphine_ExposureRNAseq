library(tidyverse)
library(edgeR)
library(limma)
library(ggpubr)
library(factoextra)
library(biomaRt)

files <- dir(path = "data_021723/", pattern = "*.tab", full.names = T) 

counttablefull <- files %>%
  map(read_tsv,  skip = 4, col_names = FALSE ) %>%
  reduce(cbind) 

#Rename data sets
datasets <-
  files %>%
  stringr::str_replace("data_021723/", "") %>% # replace data/ at the beginning of filename
  stringr::str_replace("ReadsPerGene.out.tab", "") # replace ReadsPerGene.out.tab end of name
datasets

#Rename Columns
columnnames <- c()
for (i in datasets) {
  columnnames <- c(columnnames,
                   paste0("gene_", i),
                   paste0(i, "_unstranded"),
                   paste0(i, "_forwardstrand"),
                   paste0(i, "_reversestrand")
  )
}

names(counttablefull) <- columnnames
rm(columnnames, datasets)

#Column for gene identifier, and remove duplicate gene information
counttablefull <- counttablefull %>%
  mutate(ensembl_gene =gene_CeA_Ctrl) %>%
  select(-starts_with("gene"))

counttablefull %>% head()

#check if the library is stranded
stranded_library<-counttablefull %>% # remove gene name column
  select(-ensembl_gene) %>% # find column sum
  summarise_each(funs(sum)) %>%  # gather wide data (columns) into a long table
  gather(library, counts) %>% # split library column into dataset and protocol 
  separate(library, into = c("region", "treatment", "stranding"), sep="_") %>% # split library into 3 columns
  spread(stranding,counts) %>% # make columns for each of the strandings
  mutate(propF = forwardstrand/unstranded, propR = reversestrand/unstranded) # assess rev/unst and for/unst

counttablefull %>% # remove gene name column
  select(ends_with("forwardstrand")) %>%
  summarise_each(funs(sum)) %>%  # gather wide data (columns) into a long table
  gather(library, counts) %>%  # split library column into dataset and protocol
  separate(library, into = c("region", "treatment", "stranding"), sep="_") %>%
  select(-stranding) %>%
  ggplot(aes(y = counts, x = region, fill = treatment)) +geom_bar(stat = "identity", position = "dodge") 

#Reduce library to forward stranded. 
counttablefull <- counttablefull %>% 
  select(-ends_with("unstranded"), -ends_with("reversestrand")) 
names(counttablefull) <- str_replace(names(counttablefull), "_forwardstrand", "")

dim(counttablefull)

#Check the counts per gene in sample fromm the library
counttablefull %>%
  select(-ensembl_gene) %>%
  gather(library, value) %>%
  filter(value > 0) %>% 
  ggplot(aes(x = value)) + geom_density() + facet_wrap(~library) + scale_x_log10()

#Convert count table into matrix
# best to convert to a matrix, as edgeR expects a matrix of counts or a DGEList object as input
counttablematrix <- counttablefull %>%
  select(-ensembl_gene) %>%
  as.matrix()

row.names(counttablematrix) <- counttablefull$ensembl_gene
#Save matrix
write.table(counttablematrix, "output/GeneCountMatrix_GENCODE.txt",sep="\t",quote=F)

#Transform into counts per million 
counttable_cpm <- cpm(counttablematrix)

#check possible thresholds
counttablefull %>%
  select(-ensembl_gene) %>%
  purrr::modify(~sum(.)) %>% 
  distinct() %>%
  gather(library, counts) %>%  # split library column into dataset and protocol
  mutate(millionreads = counts/(10^6)) %>%
  mutate(cpmThreshold = 10/millionreads) %>%
  filter(str_detect(library, 'Ctrl')) %>%
  mutate(mycutoff = 3.5 * millionreads)

subsetting_matrix <- counttable_cpm > 3.5
head(subsetting_matrix)

subsetting_vector <- rowSums(subsetting_matrix) >= 2
head(subsetting_vector)

counttablematrix_filt <- counttablematrix[subsetting_vector,]
dim(counttablematrix_filt)
write.table(counttablematrix, "output/Filtered_GeneCountMatrix_GENCODE.txt",sep="\t",quote=F)

mydgelist <- DGEList(counttablematrix_filt)
mydgelist_cpm <- cpm(mydgelist, log = TRUE)
boxplot(mydgelist_cpm)
abline(h=median(mydgelist_cpm),col="blue")
title("Boxplots of logCPMs (unnormalised)")

mydgelist_cpm %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(library, log2cpm, -rowname) %>%
  separate(library, into = c("region", "treatment"), sep = "_") %>%
  ggplot(aes(x = region, y = log2cpm, fill = treatment)) + geom_boxplot() +
  geom_hline(yintercept = median(mydgelist_cpm), lty = 2, color = 'blue') + theme_minimal()

#PCA
res_pca <- prcomp(t(mydgelist_cpm))
colors <- as.factor(str_split(rownames(mydgelist$samples), "_", simplify = TRUE)[,1])
fviz_pca_ind(res_pca,
             col.ind = colors,
             palette = "Set1",
             repel = TRUE,    # Avoid text overlapping
             title = "PCA plot of MOR nuerons data"
)

#DGEA
mydgelist <- calcNormFactors(mydgelist)
mydgelist$samples$norm.factors
group <- colors
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
mydgelist_voomed <- voom(mydgelist,design,plot = TRUE)

# Fit the linear model
mydgelist_fit <- lmFit(mydgelist_voomed)
# set up the contrast matrix
contrast_matrix <- makeContrasts(CeA_vs_dStr = CeA - dStr,
                                 NAc_vs_dStr = NAc - dStr,
                                 CeA_vs_NAc =  CeA - NAc,
                                 levels=design)

mydgelist_fit_contrasts <- contrasts.fit(mydgelist_fit, contrast_matrix)
mydgelist_fit_contrasts<- eBayes(mydgelist_fit_contrasts)
mydgelist_de <- decideTests(mydgelist_fit_contrasts)
summary(mydgelist_de)
vennDiagram(mydgelist_de,
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3", "orange"))



mydgelist_treat <- treat(mydgelist_fit_contrasts,lfc=1)
mydgelist_de_treat <- decideTests(mydgelist_treat)
summary(mydgelist_de_treat)


par(mfrow=c(1,2))
plotMD(mydgelist_fit_contrasts,column = 1,status=mydgelist_de[,"CeA_vs_dStr"], values = c(-1, 1), main = "CeA_vs_dStr")
volcanoplot(mydgelist_fit_contrasts,coef = 1,highlight=800, main = "CeA_vs_dStr ")

par(mfrow=c(1,2))
plotMD(mydgelist_fit_contrasts,column = 2,status=mydgelist_de[,"NAc_vs_dStr"], values = c(-1, 1), main = "NAc_vs_dStr")
volcanoplot(mydgelist_fit_contrasts,coef = 2,highlight=800, main = "NAc_vs_dStr")

par(mfrow=c(1,2))
plotMD(mydgelist_fit_contrasts,column = 3,status=mydgelist_de[,"CeA_vs_NAc"], values = c(-1, 1), main = "CeA_vs_NAc")
volcanoplot(mydgelist_fit_contrasts,coef = 3,highlight=800, main = "CeA_vs_NAc")

##get genes
final_DGEA_df <- topTreat(mydgelist_treat, coef = 1, n="Inf") %>% rownames_to_column()  %>% dplyr::rename(gene = rowname)
final_DGEA_df$significant <- final_DGEA_df$adj.P.Val <= 0.05
# use the table command to verify that the above matches the output in the venn diagram :) 
table(final_DGEA_df$significant)
# need to remove decimal point after gene id
final_DGEA_df <- final_DGEA_df %>% separate(gene, into=c("gene", "version"), sep = "\\.")

#GENE ANNOTATION
ensembl <-  useDataset("mmusculus_gene_ensembl",mart=useMart("ensembl"))
annotations <- getBM(attributes=c("ensembl_gene_id","mgi_symbol", "description", "gene_biotype"), 
                     filters = "ensembl_gene_id", 
                     values=final_DGEA_df$gene, 
                     mart=ensembl)

final_DGEA_df_anno <- inner_join(final_DGEA_df, annotations, by = c("gene" = "ensembl_gene_id"))
