library(DESeq2)
library(rtracklayer)
ref_gene = "../data/mm10/genes/genes.gtf"
first_name <- "right"
second_name <- "left"

fc <- featureCounts (files = c(
                        #right
                        "./data/bam/E11A/E11A_01_S139Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_02_S140Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_03_S141Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_04_S142Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_06_S144Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_09_S147Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_10_S148Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_13_S151Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_14_S152Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_15_S153Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_18_S156Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_19_S157Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_20_S158Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_24_S162Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_26_S164Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_33_S171Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_34_S172Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_36_S174Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_37_S175Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_38_S176Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_40_S178Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_44_S182Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_48_S186Aligned.sortedByCoord.out.bam", #E11A limit
                        "./data/bam/E12F/E12F_01_S203Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_02_S204Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_03_S205Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_04_S206Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_05_S207Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_06_S208Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_07_S209Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_08_S210Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_09_S211Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_10_S212Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_12_S214Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_15_S217Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_16_S218Aligned.sortedByCoord.out.bam", #E12F limit
                        "./data/bam/E14F/E14F_05_S223Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_06_S224Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_09_S227Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_14_S232Aligned.sortedByCoord.out.bam",
    
                        #top
                        "./data/bam/E14F/E14F_03_S221Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_08_S226Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_13_S231Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_16_S234Aligned.sortedByCoord.out.bam", #right limit
                        
                        #left
                        "./data/bam/E11A/E11A_05_S143Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_07_S145Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_08_S146Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_11_S149Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_16_S154Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_17_S155Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_21_S159Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_22_S160Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_23_S161Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_25_S163Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_27_S165Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_28_S166Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_29_S167Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_30_S168Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_31_S169Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_32_S170Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_02_S220Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_35_S173Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_39_S177Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_41_S179Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_42_S180Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_43_S181Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_45_S183Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_46_S184Aligned.sortedByCoord.out.bam",
                        "./data/bam/E11A/E11A_47_S185Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_11_S213Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_13_S215Aligned.sortedByCoord.out.bam",
                        "./data/bam/E12F/E12F_14_S216Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_01_S219Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_04_S222Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_07_S225Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_10_S228Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_11_S229Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_12_S230Aligned.sortedByCoord.out.bam",
                        "./data/bam/E14F/E14F_15_S233Aligned.sortedByCoord.out.bam"
              
                        
                        #center
#                         "./data/bam/E11A/E11A_12_S150Aligned.sortedByCoord.out.bam",
                        

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
group <- data.frame(con = factor(c(firstgroup,firstgroup,firstgroup,firstgroup,firstgroup,
                                   firstgroup,firstgroup,firstgroup,firstgroup,firstgroup,
                                   firstgroup,firstgroup,firstgroup,firstgroup,firstgroup,
                                   firstgroup,firstgroup,firstgroup,firstgroup,firstgroup,
                                   firstgroup,firstgroup,firstgroup,firstgroup,firstgroup,
                                   firstgroup,firstgroup,firstgroup,firstgroup,firstgroup,
                                   firstgroup,firstgroup,firstgroup,firstgroup,firstgroup,
                                   firstgroup,firstgroup,firstgroup, #38
                                   firstgroup,firstgroup,firstgroup,firstgroup, #top4
                                   
                                   secondgroup,secondgroup,secondgroup,secondgroup,secondgroup,
                                   secondgroup,secondgroup,secondgroup,secondgroup,secondgroup,
                                   secondgroup,secondgroup,secondgroup,secondgroup,secondgroup,
                                   secondgroup,secondgroup,secondgroup,secondgroup,secondgroup,
                                   secondgroup,secondgroup,secondgroup,secondgroup,secondgroup,
                                   secondgroup,secondgroup,secondgroup,secondgroup,secondgroup,
                                   secondgroup,secondgroup,secondgroup,secondgroup,secondgroup,
                                   secondgroup,secondgroup))) #37
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


