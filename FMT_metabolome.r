library(dplyr); library(metadeconfoundR); library(ggplot2); library("stringr"); library(vegan)
setwd("/Users/fjw536/Desktop/Anorexia/AN_DA/AN_FMT_mice_plasma_metabolome/")
metabolome <- read.delim('Metabolome_overlap.txt', header = T, row.names = 1) %>% t %>% data.frame
names(metabolome)<-sapply(str_remove_all(colnames(metabolome),"X"),"[")
meta <- read.delim('Meta.txt', header = T, row.names = 1)
# anno <- read.delim('Met_annotation.txt', header = T, row.names = 1)

meta$Group <- as.factor(meta$Group)
meta$Batch <- as.factor(meta$Batch)

# For the whole dataset
res <- matrix(NA, nrow = ncol(metabolome), ncol = 5) %>% data.frame
rownames(res) <- colnames(metabolome)
colnames(res) <-  c('Log2FC_D', 'p_D', 'p_R', 'log2FC_R', 'Color')
for (i in 1:nrow(res)) {
  res[i, 1] <- log2(mean(as.numeric(metabolome[c(1,3,5),i]))/mean(as.numeric(metabolome[c(2,4,6),i])))
  res[i, 2] <- wilcox.test(as.numeric(metabolome[c(1,3,5),i]), as.numeric(metabolome[c(2,4,6),i]))$p.value
  res[i, 3] <- wilcox.test(as.numeric(metabolome[c(7:10, 15:17, 21:23), i]), as.numeric(metabolome[c(11:14, 18:20, 24:26), i]))$p.value
  res[i, 4] <- log2(mean(as.numeric(metabolome[c(7:10, 15:17, 21:23), i]))/mean(as.numeric(metabolome[c(11:14, 18:20, 24:26), i])))
  if (res[i,1] > 0 & res[i, 4] > 0) {
    res[i, 5] <- 'Yes'
  } else if (res[i, 1] < 0 & res[i, 4] < 0) {
    res[i, 5] <- 'Yes'
  } else {
    res[i, 5] <- 'No'
  }
}
table(res$Color)
res <- res[order(-res$Log2FC_D),]


combine_D <- ggplot(res, aes(x=reorder(rownames(res), +Log2FC_D), y=Log2FC_D, fill = as.factor(Color))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('grey', '#3288ff'))+
  coord_flip() + theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 0),
        axis.text = element_text(size = 12,color = 'black'),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  geom_hline(yintercept = 0)


combine_R <- ggplot(res, aes(x=reorder(rownames(res), +Log2FC_D), y=log2FC_R, fill = as.factor(Color))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('grey', '#3288ff'))+
  coord_flip() + theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 0),
        axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 12,color = 'black'),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  geom_hline(yintercept = 0)


library(ggpubr)

p_combine <- ggarrange(combine_D, combine_R, ncol = 2, widths = c(1, 0.55))
p_combine


p_metabolome <- ggarrange(p_b1, p_b3, p_b4, ncol = 3)

p_metabolome

pdf('/Users/fjw536/Desktop/Anorexia/AN_DA/FMT_metabolome.pdf', height = 8, width = 6)
print(p_combine)
dev.off()
