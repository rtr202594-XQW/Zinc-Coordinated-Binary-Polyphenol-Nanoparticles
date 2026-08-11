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
#Fig9d
library(ggpubr)
data <- read.table("β_diversity_index_raw_data.tsv",header = T,sep = "\t",stringsAsFactors = F,check.names = FALSE)
data$Group2 <- ifelse(data$Group2 == "Con", "Control", 
                      ifelse(data$Group2 == "ZCAh", "ZnCA(H)",
                             ifelse(data$Group2 == "ZCAm", "ZnCA(M)", 
                                    ifelse(data$Group2 == "ZnCA", "Zn+CA", data$Group2))))
data_filter <- data[c(1:110),]
data_filter$Group2 <- factor(data_filter$Group2,levels=c("Control","DSS","Zn+CA","ZnCA(M)","ZnCA(H)"))
Fig9d <- ggplot(data_filter, aes(x = Group2, y = Distance,fill = Group2)) +
  geom_boxplot() +  
  geom_jitter(width = 0.005) +
  labs(title = "Distance to Control",x = "", y = "Bray-Curtis Distance") +
  theme_bw() + scale_fill_manual(values = c("#66c2a5","#fc8d62","#e78ac3","#8da0cb","#a6d854")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),legend.position = "none") +
  ylim(0, 1.2) +
  stat_compare_means(
  comparisons = list(c(1,2), c(2,3), c(2,4), c(2,5),
                     c(3,4), c(3,5)),  
  method = "wilcox.test",  
  label = "p.signif",  
  hide.ns = FALSE,  
  tip.length = 0.01,  
  vjust = 0.5
)
ggsave("Fig9d.png",Fig9d,width = 3.5,height = 3.5)

#Fig9e
library(phyloseq)
create_stacked_barplot <- function(features_qza, taxonomy_qza, tax_level, metadata_file, top_n, output_file) {
  metadata <- read.delim(metadata_file, sep = '', stringsAsFactors = FALSE, fileEncoding = "UTF-16LE")
  physeq <- qza_to_phyloseq(
    features = features_qza,
    taxonomy = taxonomy_qza
  )
  physeq_tax_level <- tax_glom(physeq, tax_level, NArm = TRUE)
  physeq_tax_level_rel <- transform_sample_counts(physeq_tax_level, function(x) x / sum(x) * 100)
  tax_level_data <- psmelt(physeq_tax_level_rel) %>%
    group_by(Sample, !!sym(tax_level)) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup()
  top_taxa <- tax_level_data %>%
    group_by(!!sym(tax_level)) %>%
    summarize(TotalAbundance = mean(Abundance)) %>%
    arrange(desc(TotalAbundance)) %>%
    head(top_n) %>%
    pull(!!sym(tax_level))
  
  tax_level_data <- tax_level_data %>%
    mutate(!!sym(tax_level) := if_else(!!sym(tax_level) %in% top_taxa, !!sym(tax_level), "Others")) %>%
    group_by(Sample, !!sym(tax_level)) %>%
    summarize(Abundance = sum(Abundance))
  tax_level_data <- tax_level_data %>%
    left_join(metadata, by = c("Sample" = "SampleID"))
  tax_level_order <- tax_level_data %>%
    group_by(!!sym(tax_level)) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    arrange(desc(TotalAbundance))
  
  if ("Others" %in% tax_level_order[[tax_level]]) {
    other_row <- tax_level_order %>% filter(!!sym(tax_level) == "Others")
    tax_level_order <- tax_level_order %>% filter(!!sym(tax_level) != "Others")
    tax_level_order <- bind_rows(tax_level_order, other_row)
  }
  
  tax_level_data[[tax_level]] <- factor(tax_level_data[[tax_level]], levels = tax_level_order[[tax_level]])
  tax_level_data$Group <- ifelse(tax_level_data$Group == "Con", "Control", tax_level_data$Group)
  tax_level_data$Group <- ifelse(tax_level_data$Group == "ZnCA", "Zn+CA", tax_level_data$Group)
  tax_level_data$Group <- ifelse(tax_level_data$Group == "ZCAm", "ZnCA (M)", tax_level_data$Group)
  tax_level_data$Group <- ifelse(tax_level_data$Group == "ZCAh", "ZnCA (H)", tax_level_data$Group)
  tax_level_data$Group <- factor(tax_level_data$Group, levels = c("Control", "DSS", "Zn+CA", "ZnCA (M)", "ZnCA (H)"))
  tax_level_levels <- unique(tax_level_data[[tax_level]])
  color_map <- c(
    "Muribaculaceae"="#e05759",
    "Bacteroidaceae"="#8dd3c7",
    "Lachnospiraceae"="#b3de69",
    "Oscillospiraceae_88309"="#8cd17d",
    "CAG-508"="#4d79a6",
    "Ruminococcaceae"="#ff9d99",
    "Rikenellaceae"="#ccebc5",
    "Desulfonisporaceae"="#ffffb3",
    "Lactobacillaceae"="#fccde5",
    "Clostridiaceae_222000"="#80b1d3",
    "Nanosyncoccaceae"="#fdb462",
    "Acutalibacteraceae"="#fb8072",
    "Tannerellaceae"="#ffed6f",
    "Coprobacillaceae"="#bebada",
    "Eggerthellaceae"="#6eb800",
    "Others"="#d9d9d9"
  )
  if ("Others" %in% tax_level_levels) {
    color_map["Others"] <- "#d9d9d9"
  }
  p1 <- ggplot(tax_level_data, aes(x = Sample, y = Abundance, fill = !!sym(tax_level))) +
    geom_col(position = 'stack', width = 0.6) +
    scale_y_continuous(expand = c(0.002, 0)) +
    scale_fill_manual(values = color_map) +
    labs(x = 'Samples', y = 'Relative Abundance(%)') +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          strip.text = element_text(size = 12)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), 
          axis.title = element_text(size = 13), legend.title = element_text(size = 13, face = "bold"), 
          legend.text = element_text(size = 11,face = "italic"), legend.position = "right",
          legend.key = element_rect(fill = NA, color = NA))
  
  p2 <- p1 + facet_wrap(~Group, scales = 'free_x', ncol = 5) +
    theme(strip.text = element_text(color = "black", size = 12),
          strip.background = element_rect(color = "black", fill = "grey90"))
  ggsave(output_file, plot = p2, width = 10, height = 6, dpi = 300)
  
  return(p2)
}
create_stacked_barplot(features_qza = "table.qza", taxonomy_qza = "taxonomy-new.qza", tax_level = "Family", metadata_file = "metadata.txt", top_n = 15, output_file = "Fig9e.png")

#Fig9f
library(tidyverse)
library(ggplot2)
metadata <- read.delim('metadata.txt', sep = '', stringsAsFactors = FALSE,fileEncoding = "UTF-16LE")
rownames(metadata) <- metadata$SampleID
plot_otu <- read.csv("plot_otu.csv",row.names = 1)
c <- data.frame(t(plot_otu), check.names = FALSE)
plot_otu_intersect_table <- merge(metadata,c,by.x = "row.names",by.y = "row.names",all.x = TRUE)
plot_otu_intersect_table_gather <- plot_otu_intersect_table %>%
  gather(key=type, value=Relative_abundance, -SampleID, -Group,-Row.names)
table(plot_otu_intersect_table_gather$type)
plot_otu_intersect_table_gather$Group <- ifelse(plot_otu_intersect_table_gather$Group == "Con", "Control", 
                                                ifelse(plot_otu_intersect_table_gather$Group == "ZCAh", "ZnCA(H)",
                                                       ifelse(plot_otu_intersect_table_gather$Group == "ZCAm", "ZnCA(M)", 
                                                              ifelse(plot_otu_intersect_table_gather$Group == "ZnCA", "Zn+CA",
                                                                     ifelse(plot_otu_intersect_table_gather$Group == "DSS", "DSS",plot_otu_intersect_table_gather$Group)))))
plot_otu_intersect_table_gather$Group <- factor(plot_otu_intersect_table_gather$Group,levels=c("Control","DSS","Zn+CA","ZnCA(M)","ZnCA(H)"))
plot_otu_intersect_table_gather$type <- sub("s__", "\ns__", plot_otu_intersect_table_gather$type)
plot_otu_intersect_table_gather$type <- factor(plot_otu_intersect_table_gather$type,levels=c("f__Muribaculaceae__1","f__Muribaculaceae; g__Sodaliphilus; \ns__Sodaliphilus pleomorphus__2",
                                                                                             "f__Muribaculaceae__3","f__Muribaculaceae__4","f__Muribaculaceae__5","f__Nanosyncoccaceae","f__Oscillospiraceae_88309__1",
                                                                                             "f__Oscillospiraceae_88309__2","f__CAG-314; g__CAG-314; \ns__CAG-314 sp900551395"))
table(plot_otu_intersect_table_gather$type)
Fig9f <- ggplot(plot_otu_intersect_table_gather, aes(x = Group, y = Relative_abundance, fill = Group)) +
  geom_boxplot() +  
  geom_jitter(width = 0.005) +
  stat_summary(
    fun = median,  
    geom = "line", 
    aes(group = 1), 
    color = "gray50", 
    size = 1,
    linetype = "solid"
  ) +
  facet_wrap(~ type, scales = "free_y", ncol = 3) +
  labs(title = "OTU",x = "", y = "Relative_abundance(%)") +
  theme_bw() + scale_fill_manual(values = c("#66c2a5","#fc8d62","#e78ac3","#8da0cb","#a6d854"),guide = "none") +
  scale_y_continuous(
    labels = function(x) x * 100  
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),strip.text = element_text(size = 8))
ggsave("Fig9f.png",Fig9f,width = 9,height = 9)

#Fig9g
library(psych)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(readxl)
library(dplyr)
library(readxl)
library(tibble)
df <- read_excel("Cytokine.xlsx") %>%
  column_to_rownames(var = names(.)[1])
df <- t(df)
plot_otu <- read.csv("Cytokine.csv",row.names = 1)
plot_otu <- t(plot_otu)
plot_otu <- plot_otu[rownames(df),]
cor <-corr.test(plot_otu, df, method = "pearson",adjust="none")
cmt <-cor$r
pmt <- cor$p
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}
mycol<-colorRampPalette(c("#1E4A75","white","#E41A1C"))(800)
Fig9g <- pheatmap(cmt,scale = "none",cluster_row = T, cluster_col = T, border=NA,
              display_numbers = pmt,fontsize_number = 12, number_color = "white",
              cellwidth = 20, cellheight =20,color=mycol)
ggsave("Fig9g.png",Fig9g,width = 8,height = 5)