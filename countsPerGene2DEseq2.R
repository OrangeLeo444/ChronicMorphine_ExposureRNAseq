ff <- list.files( path = "F:/RNA-Seq/ChronicMorphine_Exposure_raw/12-22_Batch/gene_counts", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
ff <- gsub( "ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "F:/RNA-Seq/ChronicMorphine_Exposure_raw/12-22_Batch/gene_counts/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
counts
write.table(res, "output_human/DESeq2_results_batch_correct.txt",sep="\t",quote=F)