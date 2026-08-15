library(ggplot2)

setwd("D:/")

if (!dir.exists("Volcano_Plots")) {
  dir.create("Volcano_Plots")
}

log2FC_cutoff <- 1
padj_cutoff   <- 0.05

file_list <- c("DEGs_DSS_vs_Con.csv", 
               "DEGs_ZCA_vs_Con.csv", 
               "DEGs_ZCA_vs_DSS.csv")

for (file_path in file_list) {

  if (!file.exists(file_path)) {

    next
  }
  

  res_df <- read.csv(file_path, row.names = 1, stringsAsFactors = FALSE)
  
  comp_name <- gsub(".*DEGs_|.csv", "", file_path)

  res_df$Significant <- "No Sig."
  res_df$Significant[res_df$padj < padj_cutoff & res_df$log2FoldChange > log2FC_cutoff] <- "Up"
  res_df$Significant[res_df$padj < padj_cutoff & res_df$log2FoldChange < -log2FC_cutoff] <- "Down"

  res_df$Significant <- factor(res_df$Significant, levels = c("Up", "Down", "No Sig."))

  up_count   <- sum(res_df$Significant == "Up", na.rm = TRUE)
  down_count <- sum(res_df$Significant == "Down", na.rm = TRUE)
  legend_labels <- c(paste0("Up:", up_count), paste0("Down:", down_count))

  plot_data <- res_df[!is.na(res_df$padj), ]

  plot_data$padj[plot_data$padj == 0] <- 1e-100

  max_x <- max(abs(plot_data$log2FoldChange), na.rm = TRUE)
  x_limit <- ceiling(max_x) + 1 

  p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    
    geom_point(alpha = 0.6, size = 1.2) +
    
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), 
               linetype = "dashed", color = "black", linewidth = 0.4) +

    geom_hline(yintercept = -log10(padj_cutoff), 
               linetype = "dashed", color = "black", linewidth = 0.4) +

    scale_color_manual(
      name = "Significant",
      values = c("Up" = "#E19F41",      
                 "Down" = "#1F77B4",    
                 "No Sig." = "grey85"), 
      breaks = c("Up", "Down"),         
      labels = legend_labels            
    ) +

    labs(
      title = comp_name,
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~italic(padj))  
    ) +

    theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),

      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = c(0.85, 0.85),                     
      legend.background = element_rect(fill = "transparent", color = NA), 
      legend.key = element_blank(),                       
      panel.grid.minor = element_blank()                   
    ) +

    scale_x_continuous(limits = c(-x_limit, x_limit)) +     
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 

  ggsave(filename = paste0("Volcano_Plots/Volcano_", comp_name, ".pdf"), plot = p, width = 7, height = 6)
  ggsave(filename = paste0("Volcano_Plots/Volcano_", comp_name, ".png"), plot = p, width = 7, height = 6, dpi = 300)

}

