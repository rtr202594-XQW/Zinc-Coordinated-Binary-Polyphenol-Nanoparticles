library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(tibble)
library(ggplot2)

setwd("D:/")

filename <- "DEGs_ZCA_vs_DSS.csv" 

df <- read.csv(filename, header = TRUE, row.names = 1) 
df <- df %>% rownames_to_column(var = "ENSEMBL_ID")


gene_map <- bitr(df$ENSEMBL_ID, 
                 fromType = "ENSEMBL", 
                 toType   = c("ENTREZID", "SYMBOL"), 
                 OrgDb    = org.Mm.eg.db)

df_merge <- merge(df, gene_map, by.x = "ENSEMBL_ID", by.y = "ENSEMBL")

df_merge <- df_merge %>% 
  arrange(desc(abs(log2FoldChange))) %>% 
  distinct(ENTREZID, .keep_all = TRUE)

gene_list <- df_merge$log2FoldChange
names(gene_list) <- df_merge$ENTREZID

print(df_merge[df_merge$SYMBOL == "Tnf", c("SYMBOL", "log2FoldChange")])


gene_list <- sort(gene_list, decreasing = TRUE)


set.seed(1234)

options(timeout = 600)
Sys.setenv("download.file.method" = "libcurl")

gsea_res <- gseKEGG(geneList     = gene_list,
                    organism     = "mmu",
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,
                    verbose      = FALSE)

target_idx <- grep("MAPK signaling pathway", gsea_res$Description)

if (length(target_idx) > 0) {
  
  pathway_id <- gsea_res$ID[target_idx[1]]
  pathway_name <- gsea_res$Description[target_idx[1]]
  
  nes_score <- gsea_res[pathway_id, "NES"]
  padj_val  <- gsea_res[pathway_id, "p.adjust"]
  
  my_title <- paste0(pathway_name, ": ZnCA_vs_DSS\n",
                     "(NES = ", round(nes_score, 2), 
                     ", P.adj = ", format(padj_val, digits = 3, scientific = TRUE), ")")
  

  p <- gseaplot2(gsea_res, 
                 geneSetID = pathway_id, 
                 title = my_title,
                 pvalue_table = FALSE,
                 base_size = 15,
                 ES_geom = "line")
  
  print(p)

  png("MAPK_GSEA_ZnCA_vs_DSS.png",
      width = 8,
      height = 6,
      units = "in",
      res = 600)
  print(p)
  dev.off()

  pdf("MAPK_GSEA_ZnCA_vs_DSS.pdf",
      width = 8,
      height = 6,
      useDingbats = FALSE)
  print(p)
  dev.off()
  

  
} else {
  print(head(gsea_res))
}



# NF-kappa B signaling pathway

target_idx <- grep("NF-kappa B signaling pathway", gsea_res$Description)

if (length(target_idx) > 0) {
  
  pathway_id <- gsea_res$ID[target_idx[1]]
  pathway_name <- gsea_res$Description[target_idx[1]]
  
  nes_score <- gsea_res[pathway_id, "NES"]
  padj_val  <- gsea_res[pathway_id, "p.adjust"]
  
  my_title <- paste0(pathway_name, ": ZnCA_vs_DSS\n",
                     "(NES = ", round(nes_score, 2), 
                     ", P.adj = ", format(padj_val, digits = 3, scientific = TRUE), ")")
  
  print(paste("Plotting:", my_title))
  
  p <- gseaplot2(gsea_res, 
                 geneSetID = pathway_id, 
                 title = my_title,
                 pvalue_table = FALSE,
                 base_size = 15,
                 ES_geom = "line")
  
  print(p)

  png("NFkB_GSEA_ZnCA_vs_DSS.png",
      width = 8,
      height = 6,
      units = "in",
      res = 600)
  print(p)
  dev.off()
  
  pdf("NFkB_GSEA_ZnCA_vs_DSS.pdf",
      width = 8,
      height = 6,
      useDingbats = FALSE)
  print(p)
  dev.off()
  

  
} else {
  print(head(gsea_res$Description, 20))
}


# JAK-STAT signaling pathway

target_idx <- grep("JAK-STAT signaling pathway", gsea_res$Description)

if (length(target_idx) > 0) {
  
  pathway_id <- gsea_res$ID[target_idx[1]]
  pathway_name <- gsea_res$Description[target_idx[1]]
  
  nes_score <- gsea_res[pathway_id, "NES"]
  padj_val  <- gsea_res[pathway_id, "p.adjust"]
  
  my_title <- paste0(pathway_name, ": ZnCA_vs_DSS\n",
                     "(NES = ", round(nes_score, 2),
                     ", P.adj = ", format(padj_val, digits = 3, scientific = TRUE), ")")
  
  print(paste("Plotting:", my_title))
  
  p <- gseaplot2(gsea_res,
                 geneSetID = pathway_id,
                 title = my_title,
                 pvalue_table = FALSE,
                 base_size = 15,
                 ES_geom = "line")
  
  print(p)

  png("JAK_STAT_GSEA_ZnCA_vs_DSS.png",
      width = 8,
      height = 6,
      units = "in",
      res = 600)
  print(p)
  dev.off()

  pdf("JAK_STAT_GSEA_ZnCA_vs_DSS.pdf",
      width = 8,
      height = 6,
      useDingbats = FALSE)
  print(p)
  dev.off()
  
 
  
} else {
  print(head(gsea_res$Description, 20))
}



# TNF signaling pathway 


target_idx <- grep("TNF signaling pathway", gsea_res$Description)

if (length(target_idx) > 0) {
  
  pathway_id <- gsea_res$ID[target_idx[1]]
  pathway_name <- gsea_res$Description[target_idx[1]]
  
  nes_score <- gsea_res[pathway_id, "NES"]
  padj_val  <- gsea_res[pathway_id, "p.adjust"]
  
  my_title <- paste0(pathway_name, ": ZnCA_vs_DSS\n",
                     "(NES = ", round(nes_score, 2),
                     ", P.adj = ", format(padj_val, digits = 3, scientific = TRUE), ")")
  
  print(paste("Plotting:", my_title))
  
  p <- gseaplot2(
    gsea_res,
    geneSetID = pathway_id,
    title = my_title,
    pvalue_table = FALSE,
    base_size = 15,
    ES_geom = "line"
  )
  
  print(p)

  png("TNF_GSEA_ZnCA_vs_DSS.png",
      width = 8,
      height = 6,
      units = "in",
      res = 600)
  print(p)
  dev.off()

  pdf("TNF_GSEA_ZnCA_vs_DSS.pdf",
      width = 8,
      height = 6,
      useDingbats = FALSE)
  print(p)
  dev.off()
  
  
} else {
  
  print(head(gsea_res$Description, 20))
  
}