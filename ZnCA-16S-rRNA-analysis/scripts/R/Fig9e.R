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
