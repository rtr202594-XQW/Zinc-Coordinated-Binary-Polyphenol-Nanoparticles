setwd("D:/")

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tibble)

tpm_df <- read.table("all_DEG_TPM_data.txt", header = TRUE, sep = "\t", check.names = FALSE)

up_list_df   <- read.table("UP_ZCA_VS_DSS.txt", header = TRUE)
down_list_df <- read.table("DOWN_ZCA_VS_DSS.txt", header = TRUE)
process_gene_list <- function(id_list_vector, tpm_matrix, top_n) {

  sub_data <- tpm_matrix[match(id_list_vector, tpm_matrix$geneID), ]
  sub_data <- na.omit(sub_data) 
  sub_data <- sub_data[, 1:10] 

  gene_map <- bitr(sub_data$geneID,
                   fromType = "ENSEMBL",
                   toType   = "SYMBOL",
                   OrgDb    = "org.Mm.eg.db")

  merged_data <- merge(sub_data, gene_map, by.x = "geneID", by.y = "ENSEMBL")

  merged_data <- merged_data[match(id_list_vector, merged_data$geneID), ]
  merged_data <- na.omit(merged_data) 

  merged_data <- merged_data[!duplicated(merged_data$SYMBOL), ]

  final_data <- head(merged_data, top_n)
  return(final_data)
}

final_up <- process_gene_list(up_list_df$geneID, tpm_df, 25)
final_up$Direction <- "Up" 

final_down <- process_gene_list(down_list_df$geneID, tpm_df, 25)
final_down$Direction <- "Down" 

plot_data_final <- rbind(final_up, final_down)

final_matrix <- as.data.frame(plot_data_final)
rownames(final_matrix) <- final_matrix$SYMBOL

final_matrix <- final_matrix[, c("Control_1","Control_2","Control_3",
                                 "DSS_1","DSS_2","DSS_3",
                                 "ZnCA_1","ZnCA_2","ZnCA_3")]

plot_matrix <- log2(final_matrix + 1)

annotation_col <- data.frame(
  Group = factor(rep(c("Control", "DSS", "ZnCA"), each = 3),
                 levels = c("Control", "DSS", "ZnCA"))
)
rownames(annotation_col) <- colnames(plot_matrix)

annotation_row <- data.frame(
  Cluster = factor(plot_data_final$Direction, levels = c("Up", "Down"))
)
rownames(annotation_row) <- rownames(plot_matrix)

ann_colors <- list(
  Group = c(Control = "#74add1", DSS = "#d73027", ZnCA = "#4575b4"),
  Cluster = c(Up = "#d73027", Down = "#4575b4") 
)


pheatmap(plot_matrix,
         scale = "row",
         cluster_rows = FALSE,    
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_row = annotation_row, 
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         fontsize_row = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 Genes",
         filename = "Heatmap_Top50_ZCA_VS_DSS.pdf",
         width = 6, height = 10)
