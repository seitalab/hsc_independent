library(DESeq2)
library(rtracklayer)
ref_gene = "../data/mm10/genes/genes.gtf"
first_name <- "right_end"
second_name <- "left_end"

fc <- featureCounts (files = c(
                        #right end
                        "./data/bam/E11A/E11A_13_S151Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_08_S210Aligned.sortedByCoord.out.bam",
                    
                        #left end
                        "./data/bam/E11A/E11A_29_S167Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_13_S215Aligned.sortedByCoord.out.bam"
              
                     ),
             annot.ext=ref_gene,
             isGTFAnnotationFile=TRUE,
             nthreads=8,
             allowMultiOverlap = TRUE)

print(head(fc$counts))
  
firstgroup <- paste(first_name,"_counts",sep="")
secondgroup <- paste(second_name,"_counts",sep="")

count_file_name <- paste(first_name,"_vs_",second_name,"_counts.txt",sep="")
deseq_file_name <- paste(first_name,"_vs_",second_name,"_DESeq2_result.txt", sep="")
deseq_add_genename_file_name <- paste(first_name,"_vs_",second_name,"_DESeq2_result_add_gene_name.txt", sep="")

write.table(file = count_file_name, fc$counts, quote=FALSE, sep="\t")
count_t <- t(fc$counts)
print(head(count_t))

dt <- as.matrix(fc$counts)
group <- data.frame(con = factor(c(firstgroup,firstgroup,
                                   secondgroup,secondgroup)))
dds <- DESeqDataSetFromMatrix(countData = dt, colData = group, design = ~ con)
View(counts(dds))


dds <- DESeq(dds)
res <- results(dds)
print(head(res))
res$gene_id <- row.names(res)
write.table(res,file=deseq_file_name,sep="\t",row.names=F,col.names=T,quote=F)
gtf <- readGFF(ref_gene)
gtf_gene <- subset(gtf, gtf$type == "gene")
colnames(gtf_gene)
gtf_gene_1 <- gtf_gene[,c("gene_id","gene_name","gene_biotype")]
DEseq2_result <- data.frame(res)
DESeq2_result_ref <- merge(DEseq2_result,gtf_gene_1)
write.table(DESeq2_result_ref,file=deseq_add_genename_file_name,sep="\t",row.names=F,col.names=T,quote=F)
print(head(DESeq2_result_ref))


