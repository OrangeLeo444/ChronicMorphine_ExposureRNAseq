#Script for getting gene length and GC proportions from a given data set of genes with esenmbl ID format
ensembl <-  useDataset("mmusculus_gene_ensembl",mart=useMart("ensembl"))
annotations <- getBM(attributes=c("ensembl_gene_id","mgi_symbol", "description", "gene_biotype","chromosome_name","start_position","end_position"), 
                     filters = "ensembl_gene_id", 
                     values=rownames(cts), 
                     mart=ensembl)

annotations <- annotations %>% distinct(ensembl_gene_id, .keep_all=TRUE)
write.table(annotations, "output/GeneCountMatrix_GENCODE_annotations.txt",sep="\t",quote=F)


total_features<-length(cts$gene)
upper_n <- 1000
upper_cycle<-1:floor(total_features/upper_n)
geneLegth<-data.frame()
end_cycle<-1
cycle<-1
geneLegth<-data.frame()
repeat{
  end_inx<-upper_cycle[cycle]*upper_n
  temp_length <- getGeneLengthAndGCContent(cts$gene[end_cycle:end_inx],"mmusculus_gene_ensembl")
  #temp_length<-cts$gene[end_cycle:end_inx]
  geneLegth<-rbind(geneLegth,temp_length)
  end_cycle<-end_inx+1
  cycle<-cycle+1
  if (cycle > length(upper_cycle)) break
}
temp_length <- getGeneLengthAndGCContent(cts$gene[end_cycle:total_features],"mmusculus_gene_ensembl")
geneLegth<-rbind(geneLegth,data.frame(temp_length))

write.table(geneLegth, "output/GeneCountMatrix_GENCODE_genLength_GC.txt",sep="\t",quote=F)

