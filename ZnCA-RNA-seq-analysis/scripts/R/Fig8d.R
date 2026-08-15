
library(VennDiagram)
library(grid)

setwd("D:/")

Up_ZnCA_vs_DSS <- read.csv(
  "Up_ZnCA_vs_DSS.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

Up_DSS_vs_Control <- read.csv(
  "Up_DSS_vs_Control.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

Down_ZnCA_vs_DSS <- read.csv(
  "Down_ZnCA_vs_DSS.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

Down_DSS_vs_Control <- read.csv(
  "Down_DSS_vs_Control.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

gene_up_ZnCA_DSS <- unique(na.omit(Up_ZnCA_vs_DSS[[1]]))
gene_up_DSS_Control <- unique(na.omit(Up_DSS_vs_Control[[1]]))

gene_down_ZnCA_DSS <- unique(na.omit(Down_ZnCA_vs_DSS[[1]]))
gene_down_DSS_Control <- unique(na.omit(Down_DSS_vs_Control[[1]]))

cat("Up_ZnCA_vs_DSS:", length(gene_up_ZnCA_DSS), "\n")
cat("Up_DSS_vs_Control:", length(gene_up_DSS_Control), "\n")
cat("Down_ZnCA_vs_DSS:", length(gene_down_ZnCA_DSS), "\n")
cat("Down_DSS_vs_Control:", length(gene_down_DSS_Control), "\n")

venn1 <- draw.pairwise.venn(
  
  area1 = length(gene_down_DSS_Control),
  area2 = length(gene_up_ZnCA_DSS),
  
  cross.area = length(intersect(
    gene_down_DSS_Control,
    gene_up_ZnCA_DSS
  )),
  
  category = c(
    "Down_DSS_vs_Control",
    "Up_ZnCA_vs_DSS"
  ),

  fill = c("#F6BE9A", "#A9E5CD"),
  
  alpha = c(0.80, 0.80),

  col = c("#888888", "#888888"),
  lwd = 1.5,
  lty = "solid",

  cex = 1.5,
  fontface = "plain",
  fontfamily = "sans",

  cat.cex = 1.25,
  cat.fontface = "plain",
  cat.fontfamily = "sans",

  cat.pos = c(-25, 25),
  cat.dist = c(0.05, 0.05),

  rotation.degree = 0,
  
  scaled = FALSE
)

pdf(
  "Venn_Down_DSS_vs_Control_Up_ZnCA_vs_DSS.pdf",
  width = 6,
  height = 5
)

grid.draw(venn1)

dev.off()

tiff(
  "Venn_Down_DSS_vs_Control_Up_ZnCA_vs_DSS.tiff",
  width = 1800,
  height = 1500,
  res = 300,
  compression = "lzw"
)

grid.draw(venn1)

dev.off()

venn2 <- draw.pairwise.venn(
  
  area1 = length(gene_up_DSS_Control),
  area2 = length(gene_down_ZnCA_DSS),
  
  cross.area = length(intersect(
    gene_up_DSS_Control,
    gene_down_ZnCA_DSS
  )),
  
  category = c(
    "Up_DSS_vs_Control",
    "Down_ZnCA_vs_DSS"
  ),
  
  fill = c("#F6BE9A", "#A9E5CD"),
  
  alpha = c(0.80, 0.80),
  
  col = c("#888888", "#888888"),
  lwd = 1.5,
  lty = "solid",
  
  cex = 1.5,
  fontface = "plain",
  fontfamily = "sans",
  
  cat.cex = 1.25,
  cat.fontface = "plain",
  cat.fontfamily = "sans",
  
  cat.pos = c(-25, 25),
  cat.dist = c(0.05, 0.05),
  
  rotation.degree = 0,
  
  scaled = FALSE
)

pdf(
  "Venn_Up_DSS_vs_Control_Down_ZnCA_vs_DSS.pdf",
  width = 6,
  height = 5
)

grid.draw(venn2)

dev.off()

tiff(
  "Venn_Up_DSS_vs_Control_Down_ZnCA_vs_DSS.tiff",
  width = 1800,
  height = 1500,
  res = 300,
  compression = "lzw"
)

grid.draw(venn2)

dev.off()

reverse_up_genes <- intersect(
  gene_down_DSS_Control,
  gene_up_ZnCA_DSS
)

reverse_down_genes <- intersect(
  gene_up_DSS_Control,
  gene_down_ZnCA_DSS
)


write.csv(
  reverse_up_genes,
  "Overlap_Down_DSS_Up_ZnCA.csv",
  row.names = FALSE
)

write.csv(
  reverse_down_genes,
  "Overlap_Up_DSS_Down_ZnCA.csv",
  row.names = FALSE
)

cat(
  "\nDown_DSS_vs_Control ∩ Up_ZnCA_vs_DSS =",
  length(reverse_up_genes),
  "\n"
)

cat(
  "Up_DSS_vs_Control ∩ Down_ZnCA_vs_DSS =",
  length(reverse_down_genes),
  "\n"
)