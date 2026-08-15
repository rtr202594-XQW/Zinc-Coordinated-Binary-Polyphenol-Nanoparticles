library(DESeq2)
library(edgeR) 
setwd("D:/")
count_matrix <- read.table("My_Analysis_Results.isoform.counts.matrix.txt", header = TRUE, row.names = 1, sep = "\t")
count_matrix <- round(count_matrix)
sample_info <- read.table("metadata.txt", header = TRUE, row.names = 1, sep = "\t")

if(!all(rownames(sample_info) == colnames(count_matrix))){
  stop("Error: Sample names in the metadata do not match the column names of the count matrix, or the sample order is inconsistent!")
}
sample_info$Group <- factor(sample_info$Group, levels = c("Con", "DSS", "ZCA"))
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Group)
mycounts <- counts(dds)
cpm_matrix <- cpm(mycounts)
min_samples_with_cpm_gt_1 <- 4
keep <- rowSums(cpm_matrix > 1) >= min_samples_with_cpm_gt_1
dds <- dds[keep, ]
cat(paste("After filtering,", nrow(dds), "genes remain for downstream analysis.\n"))
dds <- DESeq(dds)
comparisons <- list(
  "DSS_vs_Con"  = c("Group", "DSS", "Con"),
  "ZCA_vs_DSS"  = c("Group", "ZCA", "DSS"),
  "ZCA_vs_Con"  = c("Group", "ZCA", "Con")
)
if (!dir.exists("DEG_Results_pvalue")) {
  dir.create("DEG_Results_pvalue")
}

for (comp_name in names(comparisons)) {
  
  
  res <- results(dds, contrast = comparisons[[comp_name]])
  
  
  res_ordered <- res[order(res$pvalue), ]
  
  
  res_df <- as.data.frame(res_ordered)
  
  
  file_name <- paste0("DEG_Results_pvalue/DEGs_", comp_name, ".csv")
  write.csv(res_df, file = file_name)
 
  
  significant_degs <- subset(res_df, pvalue < 0.05 & abs(log2FoldChange) > 1)
  sig_file_name <- paste0("DEG_Results_pvalue/Significant_DEGs_", comp_name, ".csv")
  write.csv(significant_degs, file = sig_file_name)
  

  up_regulated_degs <- subset(significant_degs, log2FoldChange > 1)
  up_file_name <- paste0("DEG_Results_pvalue/Up_regulated_DEGs_", comp_name, ".csv")
  write.csv(up_regulated_degs, file = up_file_name)
  

  down_regulated_degs <- subset(significant_degs, log2FoldChange < -1)
  down_file_name <- paste0("DEG_Results_pvalue/Down_regulated_DEGs_", comp_name, ".csv")
  write.csv(down_regulated_degs, file = down_file_name)
  

  cat(paste("Completed comparison:", comp_name, "\n"))
  cat(paste("- Total significant DEGs:", nrow(significant_degs), "\n"))
  cat(paste("- Upregulated genes:", nrow(up_regulated_degs), "\n"))
  cat(paste("- Downregulated genes:", nrow(down_regulated_degs), "\n\n"))
}
