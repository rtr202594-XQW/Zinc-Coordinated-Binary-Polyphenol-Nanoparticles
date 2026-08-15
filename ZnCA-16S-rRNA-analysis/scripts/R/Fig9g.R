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