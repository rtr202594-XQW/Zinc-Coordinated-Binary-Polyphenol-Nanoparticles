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