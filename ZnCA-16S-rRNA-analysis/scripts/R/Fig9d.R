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