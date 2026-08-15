library(qiime2R)
library(dplyr)
library(tidyverse)
library(vegan)
library(ggExtra)
library(ggplot2)
#Fig9c
perform_pcoa_analysis <- function(metadata_file, distance_matrix_file, output_file) {
  group <- read.delim(metadata_file, sep = '', stringsAsFactors = FALSE, fileEncoding = "UTF-16LE")
  bray_curtis_distance_matrix <- read_qza(distance_matrix_file)
  bray_curtis_distance_matrix_data <- bray_curtis_distance_matrix$data
  pcoa <- cmdscale(bray_curtis_distance_matrix_data, k = 3, eig = TRUE)
  pcoa_data <- data.frame(pcoa$points)
  pcoa_data$SampleID <- rownames(pcoa_data)
  names(pcoa_data)[1:3] <- paste0("PCoA", 1:3)
  eig <- pcoa$eig
  eig_percent <- round(pcoa$eig / sum(eig) * 100, 1)
  pcoa_result <- merge(pcoa_data, group, by = "SampleID")
  dune.div <- adonis2(bray_curtis_distance_matrix_data ~ Group, data = group, permutations = 999)
  dune_adonis <- paste0("adonis R2: ", round(dune.div$R2, 2), "; p_value: ", dune.div$`Pr(>F)`)
  pcoa_result$Group <- ifelse(pcoa_result$Group == "Con", "Control", pcoa_result$Group)
  pcoa_result$Group <- ifelse(pcoa_result$Group == "ZnCA", "Zn+CA", pcoa_result$Group)
  pcoa_result$Group <- ifelse(pcoa_result$Group == "ZCAm", "ZnCA (M)", pcoa_result$Group)
  pcoa_result$Group <- ifelse(pcoa_result$Group == "ZCAh", "ZnCA (H)", pcoa_result$Group)
  pcoa_result$Group <- factor(pcoa_result$Group, levels = c("Control", "DSS", "Zn+CA", "ZnCA (M)", "ZnCA (H)"))
  p <- ggplot(pcoa_result, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(aes(color = Group)) +
    theme_bw() +
    labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
         y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
         caption = dune_adonis) +
    scale_colour_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb", "#a6d854"),
                        guide = guide_legend(override.aes = list(size = 3))) +
    theme(legend.position = c(0.2, 0.85),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.25, "cm"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black", size = 10)) +
    geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed") +
    geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")
  p1 <- p + stat_ellipse(data = pcoa_result,
                         geom = "polygon",
                         level = 0.9,
                         linetype = 0,
                         linewidth = 0.5,
                         aes(fill = Group),
                         alpha = 0.1,
                         show.legend = FALSE) +
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb", "#a6d854"))
  ggsave(output_file, p1, width = 4, height = 4,dpi = 300)
}

#jaccard
perform_pcoa_analysis(metadata_file = 'metadata.txt',
                      distance_matrix_file = 'jaccard_distance_matrix.qza',
                      output_file = 'Fig9c.png')