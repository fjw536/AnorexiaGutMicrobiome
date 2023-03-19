################################################################################
########################### Data Visualization #################################
################################################################################

# Load required libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(ade4)
library(reshape2)
library(ggpubr)
library(metadeconfoundR)
library(ggcharts)
library(rcompanion)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(ggalluvial)
library(wesanderson)
library(ggfortify)
library(fossil)
library(gridExtra)
library(tidyverse)
library(VennDiagram)
source('Functions.r')

# Figure 1a - box plot of beta diversity between AN and NC
# Load datasets
MSP <- read.delim("MSP_matrix.txt", header = T, row.names = 1) %>% t %>% as.data.frame()
meta <- read.delim("Pheno.txt", header = T, row.names = 1)
MSP.anno <- read.delim("IGC2.1989MSPs.taxo.tsv", header = T, row.names = 1)

## Check rowname consistency
MSP <- MSP[rownames(meta), ] 
rownames(MSP) == rownames(meta)

## Remove MSP with low prevalence
MSP.filter <- filter_pre(MSP, 0.1, 'row')

## Define subgroups
an <- meta %>% filter(Category == 'AN') %>% rownames()
nc <- meta %>% filter(Category == 'NC') %>% rownames()
IN <- meta %>% filter(in.out.patient == '1') %>% rownames()
out <- meta %>% filter(in.out.patient == '2') %>% rownames()

## Compute the distance
final.dis <- as.matrix(vegdist(MSP.filter, method = "canberra"))
final.dis.mds<-cmdscale(final.dis, k=5, eig = T)
msp_pcoa <- data.frame(final.dis.mds$points)

## Beta diversity of MSP matrix - box plot
m2 <- melt(final.dis)[melt(upper.tri(final.dis))$value,]
names(m2) <- c("c1", "c2", "distance")
an <- meta %>% filter(Category == 'AN') %>% rownames()
nc <- meta %>% filter(Category == 'NC') %>% rownames()
IN <- meta %>% filter(in.out.patient == '1') %>% rownames()
out <- meta %>% filter(in.out.patient == '2') %>% rownames()


an.row <- c()
for (i in 1:nrow(m2)) {
  if (m2[i, 1] %in% an & m2[i,2] %in% an) {
    an.row <- c(an.row, i)
  }
}

nc.row <- c()
for (i in 1:nrow(m2)) {
  if (m2[i, 1] %in% nc & m2[i,2] %in% nc) {
    nc.row <- c(nc.row, i)
  }
}


IN.row <- c()
for (i in 1:nrow(m2)) {
  if (m2[i, 1] %in% IN & m2[i,2] %in% IN) {
    IN.row <- c(IN.row, i)
  }
}

out.row <- c()
for (i in 1:nrow(m2)) {
  if (m2[i, 1] %in% out & m2[i,2] %in% out) {
    out.row <- c(out.row, i)
  }
}


wilcox.test(m2[nc.row, 3], m2[an.row, 3])
wilcox.test(m2[bp.row, 3], m2[rs.row, 3])
wilcox.test(m2[IN.row, 3], m2[out.row, 3])


m2[nc.row, 'Group'] <- 'NC'
m2[an.row, 'Group'] <- 'AN'
m2.nona <- na.omit(m2)
m2.nona$Group <- factor(m2.nona$Group, levels = c('NC', 'AN'))
m2 <- m2[, -ncol(m2)]
msp.beta.box <- ggplot(m2.nona, aes(Group, distance, fill = Group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c( "#5569a5", "#be8738"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = 'black'),
        axis.text.x = element_text(size = 16, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.position = 'none',
        axis.line = element_line(color = 'black'))+
  ylab('Canberra dissimilarity')+
  annotate("text", x=1.5, y=1, label= "P < 2.2e-16") 

pdf('msp.beta.nc.an.box.pdf', height = 5, width = 3)
print(msp.beta.box)
dev.off()


# Figure 1b - box plot of beta diversity among NC, AN-RS, and AN-BP
df <- melt(as.matrix(final.dis), varnames = c('row', 'col'))
nc.df <- subset(df, df$row %in% nc) %>% filter(col %in% nc) %>% subset(!(value == 0))
bp.df <- subset(df, df$row %in% bp) %>% filter(col %in% bp) %>% subset(!(value == 0))
rs.df <- subset(df, df$row %in% rs) %>% filter(col %in% rs) %>% subset(!(value == 0))
an.df <- subset(df, df$row %in% rs) %>% filter(col %in% an) %>% subset(!(value == 0))

nc.df$'Group' <- 'NC'
bp.df$'Group' <- 'AN-BP'
rs.df$'Group' <- 'AN-RS'
an.df$'Group' <- 'AN'

beta <- data.frame(rbind(nc.df, rs.df, bp.df))


compare_means(value~Group, data=beta)
my_color <- c("#78a0c1", "#a8ded8", "#fae3b5")
my_comparisons <- list(c("NC", "AN-RS"), c("NC", "AN-BP"), c("AN-RS", "AN-BP"))
p_beta <- ggboxplot(beta, x="Group", y="value", fill = "Group",palette = my_color, lwd = 1)+
  stat_compare_means(comparisons = my_comparisons, size =5)+
  theme(axis.title = element_text(size = 22), legend.position = 'none', axis.text = element_text(size = 18)) +
  labs(x = NULL, y = "Canberra distance")

ggsave(p_beta, filename = 'bac.subtype.beta.pdf', width = 5, height = 5)


# FIgure 1c - Significantly contrasted bacterial species between AN and NC
## Data cleaning for drug deconfounding
MSP.final <- read.delim('Datasets_final/Data/MSP.final.txt', header = T, row.names = 1)
meta <- read.delim("Pheno.txt", header = T, row.names = 1)
drug_use <- read.delim('Drug_use.txt', header = T, row.names = 1)

## Remove MSPs with prevalence lower than 10%
msp.filter <- filter_pre(MSP.final, 0.1, 'row')


## Select drug intake as sub metadata
meta.sub <- meta[, c(10, 2, 28, 29, 33:48)]
meta.sub$Category <- as.factor(meta.sub$Category)
levels(meta.sub$Category) <- c('1', '0')

## Perfomr MSPs drug deconfounding
## Check if the rownames are matched between msp and meta.sub matrix
all(rownames(msp.filter)==rownames(meta.sub))
all(order(rownames(msp.filter)) == order(rownames(meta.sub)))

AN_drugdeconfound_output_msp <- MetaDeconfound(featureMat = msp.filter,
                                               metaMat = meta.sub)

msp.tmp <- AN_drugdeconfound_output_msp$status %>% as.data.frame()
msp.sig <-  msp.tmp %>% filter(Category == 'OK_nc' | Category == 'OK_sd' | Category == 'AD') %>% rownames()
msp.cliff <- matrix(NA, nrow = length(msp.sig), ncol = 6) %>% as.data.frame()
rownames(msp.cliff) <- msp.sig
colnames(msp.cliff) <- c('q_value', 'cliff', 'Prevalence_all', 'Prevalence_NC', 'Prevalence_AN', 'msp')


msp.cliff$q_value <- formatC(AN_drugdeconfound_output_msp$Qs[msp.sig, 'Category'], format = "e", digits = 1)

for (i in msp.sig) {
  msp.cliff[i, 'cliff'] <- cliff.delta(MSP.final[an, i],MSP.final[nc, i], return.dm=TRUE)$estimate
  msp.cliff[i, 'Prevalence_all'] <- paste(round(sum(MSP.final[, i] > 0)/147*100), '%', sep = '')
  msp.cliff[i, 'Prevalence_NC'] <- paste(round(sum(MSP.final[nc, i] > 0)/70*100), '%', sep = '')
  msp.cliff[i, 'Prevalence_AN'] <- paste(round(sum(MSP.final[an, i] > 0)/77*100), '%', sep = '')
}

msp.cliff$msp <- MSP.anno[rownames(msp.cliff), "species"]

msp.cliff$"msp.annotation" <- paste(rownames(msp.cliff), ': ',msp.cliff$msp, ' [', msp.cliff$Prevalence_all, ', ', msp.cliff$Prevalence_NC, ', ', msp.cliff$Prevalence_AN, ', ', msp.cliff$q_value,']',sep = '')


## Visualization
p_bac_cliff <- diverging_lollipop_chart(msp.cliff, 
                                        msp.annotation, 
                                        cliff, 
                                        lollipop_colors = c("#be8738", "#5569a5"), 
                                        line_size = 1, 
                                        point_size = 2, 
                                        text_size = 13, 
                                        text_color = c("#be8738", "#5569a5"))+
  labs(x = NULL, y = "Cliff's Delta")+ theme_light()+
  theme(axis.title = element_text(size = 18), 
        panel.grid = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ylim(-0.6, 0.6)
p_bac_cliff

ggsave(p_bac_cliff, filename = "Contrast.msps.pdf", height = 15, width = 15)

# Figure 1d - heatmap of liner regression outcomes between msp and eating disorder scores in AN
## Perform liner regression analysis
cov <- c('Age','BMI', 'Smoking', 'SSRI', 'Antipsychotics', 'Benzodiazepines')
msp_ED_lm <- lm_btw_mats(mat0 = MSP.final[an,], 
                         mat1 = meta[an, EDI], 
                         mat2 = meta[an, ], 
                         covar = cov,
                         direction = c(1,1),
                         y_mat = 0)

msp.ED.lm <- matrix(0, nrow = nrow(msp_ED_lm$p), ncol = ncol(msp_ED_lm$p)) %>% data.frame
rownames(msp.ED.lm) <- rownames(msp_ED_lm$p)
colnames(msp.ED.lm) <- colnames(msp_ED_lm$p)
for (i in 1:nrow(msp.ED.lm)) {
  for (j in 1:ncol(msp.ED.lm)) {
    if (msp_ED_lm$fdr[i, j] < 0.05) {
      msp.ED.lm[i,j] <- msp_ED_lm$beta[i,j]
    }
  }
}


stars= matrix("", nrow = nrow(msp.ED.lm), ncol = ncol(msp.ED.lm))
rownames(stars)=rownames(msp.ED.lm)
colnames(stars)=colnames(msp.ED.lm)
for (i in 1:ncol(stars)) {
  for (j in 1:nrow(stars)) {
    if (msp_ED_lm$fdr [j,i]<0.05) {
      stars[j,i]="+"
    }
  }
}

msp.ED.lm <- t(msp.ED.lm) %>% as.data.frame()
stars <- t(stars) %>% as.data.frame()


msp.ED.plot <- msp.ED.lm[which(rowSums(msp.ED.lm) != 0), ]
stars.plot <- stars[rownames(msp.ED.plot), ]

## Compute cliff delta values for each of msps that correlate with eating disorder scores

meta$Category <- factor(meta$Category, levels = c('NC', 'AN'))
cliff.input <- data.frame(cbind(MSP.final, meta$Category))
for (i in 1:nrow(msp.ED.plot)) {
  form <- formula(paste(rownames(msp.ED.plot)[i], "~meta$Category"))
  msp.ED.plot[i, 'Cliff.delta'] <- cliffDelta(form, data = cliff.input)
  msp.ED.plot[i, 'Abs.cliff'] <- abs(msp.ED.plot[i, 'Cliff.delta'])
}

up <- msp.ED.plot %>% filter(Cliff.delta<0) %>% rownames
down <- msp.ED.plot %>% filter(Cliff.delta>0) %>% rownames
msp.ED.plot[up, 'Trend'] <- 'AN-enriched'
msp.ED.plot[down, 'Trend'] <- 'NC-enriched'
msp.ED.plot[5, 'Trend'] <- 'AN-enriched'


## Make annotation file
row_anno <- matrix(NA, nrow = nrow(msp.ED.plot), ncol = 3) %>% as.data.frame()
rownames(row_anno) <- rownames(msp.ED.plot)
colnames(row_anno) <- c('Species', 'Prevalence_AN', 'Annotation')

for (i in 1:nrow(row_anno)) {
  row_anno[i, 'Species'] <- MSP.anno[rownames(row_anno)[i], 1]
  row_anno[i, 'Prevalence_AN'] <- paste(round(sum(MSP.final[an, rownames(row_anno)[i]] > 0)/77*100), '%', sep = '')
}
row_anno$Annotation <- paste(rownames(row_anno), ': ', row_anno$Species, ' [', row_anno$Prevalence_AN, ']', sep = '')


## Visualization
set.seed(56489) 
mycol <- colorRamp2(c(-0.6,  0,  0.6),
                    c("#02CBFB",  "white", "#F20C34"))

## Rename column
colnames(msp.ED.plot) <- c('Drive for thinness',
                           'Bulimia',
                           'Body dissatisfaction',
                           'Low self-esteem',
                           'Personal alienation',
                           'Interpersonal insecurity',
                           'Interpersonal alienation',
                           'Interoceptive deficits',
                           'Emotional disregulation',
                           'Perfectionism',
                           'Ascetism',
                           'Maturity fears',
                           'Cliff.delta',
                           'Abs.cliff',
                           'Trend')


ht.msp.ED.lm <- Heatmap(as.matrix(msp.ED.plot[, 1:12]), 
                        col = mycol, 
                        rect_gp = gpar(col= "black", lwd =0.5),
                        name = "Beta coefficient",
                        right_annotation = rowAnnotation(
                          Trend = msp.ED.plot$Trend,
                          col = list(Trend = c("AN-enriched" = "#be8738", "NC-enriched"="#5569a5")
                          ),
                          na_col = "white",
                          gp = gpar(col = "black")
                        ),
                        row_labels = row_anno$Annotation,
                        show_row_names = T,
                        show_column_names = T,
                        show_row_dend = F,
                        cluster_columns = FALSE,
                        heatmap_width = unit(0.5, 'npc'),
                        row_names_max_width = unit(20, "cm"),
                        heatmap_legend_param = list(legend_direction = c("horizontal"), 
                                                    legend_width = unit(4, 'cm'), 
                                                    title_position = "topcenter",
                                                    at = seq(-0.6, 0.6, 0.3)),
                        border = "black",
                        cluster_rows = T,
                        show_heatmap_legend = T,
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10),
                        cell_fun = function(j, i, x, y, width, height, fill){
                          grid.text(sprintf("%s", stars.plot[i, j]), x, y, gp = gpar(fontsize = 12))
                        })
ht.msp.ED.lm

pdf("hm.msp.ED.lm.new.pdf", height = 12, width = 8)
draw(ht.msp.ED.lm, heatmap_legend_side = "top")
dev.off()


# Figure 2a and 2b - Richness and Shannon diversity of viral gut microbiota between AN and NC
rich <- apply(virus, 1, chao1)
alpha <- diversity(virus, index = "shannon", MARGIN = 1)
com <- cbind(rich, alpha)
rownames(com) <- rownames(virus)
com <- as.data.frame(com)
com$'subgroup' <- input_metadata$Subtype

rich.bac <- apply(MSP.final, 1, chao1)
alpha.bac <- diversity(MSP.final, index = "shannon", MARGIN = 1)
com.bac <- cbind(rich.bac, alpha.bac)
rownames(com.bac) <- rownames(virus)
com.bac <- as.data.frame(com.bac)
com.bac$'subgroup' <- input_metadata$Subtype

an <- meta %>% filter(Category=='AN') %>% rownames()
nc <- meta %>% filter(Category=='NC') %>% rownames()

com[an, 'Group'] <- c('AN')
com[nc, 'Group'] <- c('NC')
com$Group <- factor(com$Group, levels = c('NC', 'AN'))

## Visualization Figure 5A - Richness
compare_means(rich~Group, data=com)
my_color <- c( "#5569a5", "#be8738")

p_a <- ggboxplot(com, x="Group", y="rich", color = "Group",palette = my_color, add = "jitter", lwd = 1)+
  stat_compare_means(size = 0)+
  theme(axis.title = element_text(size = 13), legend.position = 'none', axis.text = element_text(size = 12)) +
  labs(x = NULL, y = "Chao1 richness of virome")+
  annotate("text", x=1.5, y=70, label= " P = 5.7e-05")
## Visualization Figure 5B - Alpha diversity
p_b <- ggboxplot(com, x="Group", y="alpha", color = "Group",palette = my_color, add = "jitter", lwd = 1)+
  stat_compare_means(size = 0)+
  theme(axis.title = element_text(size = 13), legend.position = 'none', axis.text = element_text(size = 12)) +
  labs(x = NULL, y = "Shannon diversity of virome") + ylim(0, 5)+
  annotate("text", x=1.5, y=4, label= " P = 3.6e-05")

# Figure 2C - contrasted viral gut microbiota between AN and NC
## Check if the rownames are matched between virome and meta.sub matrix
all(rownames(virus)==rownames(meta.sub))
all(order(rownames(virus)) == order(rownames(meta.sub)))
meta.sub1 <- meta.sub[, -2] # Without BMI as confounder

virus_drugdeconfound_output <- MetaDeconfound(featureMat = virus,
                                              metaMat = meta.sub1)

virus.tmp <- virus_drugdeconfound_output$status %>% as.data.frame()
virus.sig <-  virus.tmp %>% filter(Category == 'OK_nc' | Category == 'OK_sd' | Category == 'AD') %>% rownames()
virus.cliff <- matrix(NA, nrow = length(virus.sig), ncol = 6) %>% as.data.frame()
rownames(virus.cliff) <- virus.sig
colnames(virus.cliff) <- c('q_value', 'cliff', 'Prevalence_all', 'Prevalence_NC', 'Prevalence_AN', 'virus')


virus.cliff$q_value <- formatC(virus_drugdeconfound_output$Qs[virus.sig, 'Category'], format = "e", digits = 1)

for (i in virus.sig) {
  virus.cliff[i, 'cliff'] <- cliff.delta(virus[an, i],virus[nc, i], return.dm=TRUE)$estimate
  virus.cliff[i, 'Prevalence_all'] <- paste(round(sum(virus[, i] > 0)/147*100), '%', sep = '')
  virus.cliff[i, 'Prevalence_NC'] <- paste(round(sum(virus[nc, i] > 0)/70*100), '%', sep = '')
  virus.cliff[i, 'Prevalence_AN'] <- paste(round(sum(virus[an, i] > 0)/77*100), '%', sep = '')
}

virus.cliff$virus <-gsub(".", " ", rownames(virus.cliff), fixed = TRUE)

virus.cliff$"virus.annotation" <- paste(virus.cliff$virus, ' [', virus.cliff$q_value,']',sep = '')

## Visualization
p_virus_cliff <- diverging_lollipop_chart(virus.cliff, 
                                          virus.annotation, 
                                          cliff, 
                                          lollipop_colors = c("#be8738", "#5569a5"), 
                                          line_size = 1, 
                                          point_size = 2, 
                                          text_size = 13, 
                                          text_color = c("#be8738", "#5569a5"))+
  labs(x = 'Virual gut microbiota', y = "Cliff's Delta")+ theme_light()+
  theme(axis.title = element_text(size = 18), 
        panel.grid = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ylim(-0.4, 0.4)
p_virus_cliff



# Figure 2d - Numbers of SparCC in AN and NC
setwd("/Users/fjw536/Desktop/Anorexia/Anorexia.sparcc/")
AN.cor <- read.delim("AN_median_correlation.tsv", sep = "\t", header = T, row.names = 1)
NC.cor <- read.delim("NC_median_correlation.tsv", sep = "\t", header = T, row.names = 1)
BP.cor <- read.delim("BP_median_correlation.tsv", sep = "\t", header = T, row.names = 1)
RS.cor <- read.delim("RS_median_correlation.tsv", sep = "\t", header = T, row.names = 1)
AN.p <- read.delim("AN.pvalues.tsv", sep = "\t", header = T, row.names = 1)
NC.p <- read.delim("NC.pvalues.tsv", sep = "\t", header = T, row.names = 1)
BP.p <- read.delim("BP.pvalues.tsv", sep = "\t", header = T, row.names = 1)
RS.p <- read.delim("RS.pvalues.tsv", sep = "\t", header = T, row.names = 1)

AN.sparcc <- matrix(0, nrow = nrow(AN.cor), ncol = ncol(AN.cor)) %>% data.frame
rownames(AN.sparcc) <- rownames(AN.cor)
colnames(AN.sparcc) <- colnames(AN.cor)
for (i in 1:nrow(AN.cor)) {
  for (j in 1:ncol(AN.cor)) {
    if (AN.cor[i,j] > 0.4 | AN.p[i, j] < 0.05) {
      AN.sparcc[i,j] <- AN.cor[i,j]
    }
  }
}
colSums(AN.sparcc !=0) %>% data.frame %>% colSums
colSums(AN.sparcc > 0) %>% data.frame %>% colSums
colSums(AN.sparcc < 0) %>% data.frame %>% colSums

NC.sparcc <- matrix(0, nrow = nrow(NC.cor), ncol = ncol(NC.cor)) %>% data.frame
rownames(NC.sparcc) <- rownames(NC.cor)
colnames(NC.sparcc) <- colnames(NC.cor)
for (i in 1:nrow(NC.cor)) {
  for (j in 1:ncol(NC.cor)) {
    if (NC.cor[i,j] > 0.4 | NC.p[i, j] < 0.05) {
      NC.sparcc[i,j] <- NC.cor[i,j]
    }
  }
}

colSums(NC.sparcc !=0) %>% data.frame %>% colSums
colSums(NC.sparcc > 0) %>% data.frame %>% colSums
colSums(NC.sparcc < 0) %>% data.frame %>% colSums

sparcc.res <- matrix(NA, nrow = 4, ncol = 3) %>% as.data.frame()
sparcc.res[1:4, 1] <- c('NC', 'NC', 'AN', 'AN')
sparcc.res[1:4, 3] <- c('Positive correlation', 'Negative correlation', 'Positive correlation', 'Negative correlation')
sparcc.res[1,2] <- colSums(NC.sparcc > 0) %>% data.frame %>% colSums
sparcc.res[2,2] <- colSums(NC.sparcc < 0) %>% data.frame %>% colSums
sparcc.res[3,2] <- colSums(AN.sparcc > 0) %>% data.frame %>% colSums
sparcc.res[4,2] <- colSums(AN.sparcc < 0) %>% data.frame %>% colSums

## Visualziation of the outcomes between AN and NC
sparcc.res$V1 <- factor(sparcc.res$V1, levels = c('NC', 'AN'))
sparcc_AN_NC.p <- ggplot(sparcc.res, aes(x = V1, y= V2, fill = V3, label = V2)) + 
  geom_bar(stat = "identity", color = 'black') +
  scale_fill_manual(values = c("#02CBFB", "#F20C34"))+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  xlab(label = NULL) +
  ylab('Number of SparCC correlations') +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12, color = 'black'),
        legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank())

ggsave(sparcc_AN_NC.p, filename = '/Users/fjw536/Desktop/Anorexia/AN_DA/sparcc.an.nc.plot.pdf', width = 3, height = 5)


BP.sparcc <- matrix(0, nrow = nrow(BP.cor), ncol = ncol(BP.cor)) %>% data.frame
rownames(BP.sparcc) <- rownames(BP.cor)
colnames(BP.sparcc) <- colnames(BP.cor)
for (i in 1:nrow(BP.cor)) {
  for (j in 1:ncol(BP.cor)) {
    if (BP.cor[i,j] > 0.4 | BP.p[i, j] < 0.05) {
      BP.sparcc[i,j] <- BP.cor[i,j]
    }
  }
}

colSums(BP.sparcc !=0) %>% data.frame %>% colSums

RS.sparcc <- matrix(0, nrow = nrow(RS.cor), ncol = ncol(RS.cor)) %>% data.frame
rownames(RS.sparcc) <- rownames(RS.cor)
colnames(RS.sparcc) <- colnames(RS.cor)
for (i in 1:nrow(RS.cor)) {
  for (j in 1:ncol(RS.cor)) {
    if (RS.cor[i,j] > 0.4 | RS.p[i, j] < 0.05) {
      RS.sparcc[i,j] <- RS.cor[i,j]
    }
  }
}

colSums(RS.sparcc !=0) %>% data.frame %>% colSums

sparcc.sub.res <- matrix(NA, nrow = 4, ncol = 3) %>% as.data.frame()
sparcc.sub.res[1:4, 1] <- c('AN-BP', 'AN-BP', 'AN-RS', 'AN-RS')
sparcc.sub.res[1:4, 3] <- c('Positive correlation', 'Negative correlation', 'Positive correlation', 'Negative correlation')
sparcc.sub.res[1,2] <- colSums(BP.sparcc > 0) %>% data.frame %>% colSums
sparcc.sub.res[2,2] <- colSums(BP.sparcc < 0) %>% data.frame %>% colSums
sparcc.sub.res[3,2] <- colSums(RS.sparcc > 0) %>% data.frame %>% colSums
sparcc.sub.res[4,2] <- colSums(RS.sparcc < 0) %>% data.frame %>% colSums

## Visualziation of the outcomes from the subtypes
sparcc.sub.res$V1 <- factor(sparcc.sub.res$V1, levels = c('AN-BP', 'AN-RS'))
sparcc_sub.p <- ggplot(sparcc.sub.res, aes(x = V1, y= V2, fill = V3, label = V2)) + 
  geom_bar(stat = "identity", color = 'black') +
  scale_fill_manual(values = c("#02CBFB", "#F20C34"))+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  xlab(label = NULL) +
  ylab(label = NULL) +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12, color = 'black'),
        legend.position = 'bottom',
        panel.background = element_blank(),
        panel.grid = element_blank())

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


my_legend <- g_legend(sparcc_sub.p)
sparcc.p <- grid.arrange(arrangeGrob(sparcc_AN_NC.p,
                                     sparcc_sub.p + theme(legend.position="none"),
                                     nrow=1),
                         my_legend, nrow=2,heights=c(10, 1))


p_ab_sparcc <- ggarrange(p_ab, sparcc.p, nrow = 2)

## Assembly all elements of Figure 2 and export
f5 <- ggarrange(p_ab_sparcc, p_virus_cliff, ncol = 2, widths = c(1, 1.5))
ggsave(f5, filename = '/Users/fjw536/Desktop/Anorexia/AN_DA/figure5.pdf', width = 12, height = 8)


# Figure 3a - contrasted gbms between AN and NC
## Perform drug deconfounding
gbm_drugdeconfound_output<- MetaDeconfound(featureMat = input.gbm,
                                           metaMat = meta.sub)

table(gbm_drugdeconfound_output$status[, 'Category'])

gbm.tmp <- gbm_drugdeconfound_output$status %>% as.data.frame()
gbm.sig <-  gbm.tmp %>% filter(Category == 'OK_nc' | Category == 'OK_sd' | Category == 'AD') %>% rownames()
gbm.cliff <- matrix(NA, nrow = length(gbm.sig), ncol = 6) %>% as.data.frame()
rownames(gbm.cliff) <- gbm.sig
colnames(gbm.cliff) <- c('q_value', 'cliff', 'Prevalence_all', 'Prevalence_NC', 'Prevalence_AN', 'gbm')


gbm.cliff$q_value <- formatC(gbm_drugdeconfound_output$Qs[gbm.sig, 'Category'], format = "e", digits = 1)

for (i in gbm.sig) {
  gbm.cliff[i, 'cliff'] <- cliff.delta(input.gbm[an, i],input.gbm[nc, i], return.dm=TRUE)$estimate
  gbm.cliff[i, 'Prevalence_all'] <- paste(round(sum(input.gbm[, i] > 0)/147*100), '%', sep = '')
  gbm.cliff[i, 'Prevalence_NC'] <- paste(round(sum(input.gbm[nc, i] > 0)/70*100), '%', sep = '')
  gbm.cliff[i, 'Prevalence_AN'] <- paste(round(sum(input.gbm[an, i] > 0)/77*100), '%', sep = '')
}

gbm.cliff$gbm <- db@module.names[rownames(gbm.cliff), "V2"]

gbm.cliff$"gbm.annotation" <- paste(gbm.cliff$gbm, ' [', gbm.cliff$q_value,']',sep = '')


## Visualization
p_gbm <- diverging_bar_chart(gbm.cliff, gbm.annotation, cliff, bar_colors = c("#be8738", "#5569a5"), text_size = 16, text_color = c("#be8738", "#5569a5"))+
  labs(x = "Gut-brain modules", y = "Cliff's Delta")+theme_light()+
  theme(axis.title = element_text(size = 16), 
        panel.grid = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ylim(-0.5, 0.5)
p_gbm

pdf("contrast.gbms.pdf", width = 12, height = 7)
print(p_gbm)
dev.off()



# Figure 3b - heatmap visualizing the linear regression model between gbms and bioclinical variables
## Association between GBMs and phenotypes
cov.1 <- c('Age', 'Smoking', 'SSRI', 'Antipsychotics', 'Benzodiazepines')
rownames(input.gbm) == rownames(meta)

gbm_pheno.only <- lm_btw_mats(input.gbm, meta[, c(2:9)], meta, cov.1, direction = c(1,1), y_mat = 1)
gbm_pheno.only.sig <- gbm_pheno.only$table %>% filter(fdr.p < 0.05)

gbm.pheno.cor <- reshape::cast(gbm_pheno.only.sig, Phenotype ~ Taxa, mean, value = 'Beta')
gbm.pheno.cor[is.na(gbm.pheno.cor)] <- 0
rownames(gbm.pheno.cor) <- gbm.pheno.cor[,1]
gbm.pheno.cor <- gbm.pheno.cor[, -1]
star <- reshape::cast(gbm_pheno.only.sig, Phenotype ~ Taxa , mean, value = 'fdr.p')
star[is.na(star)] <- ''
rownames(star) <- star[,1]
star <- star[, -1]
star[star > 0] <- '+'


rownames(gbm.pheno.cor) <- db@module.names[rownames(gbm.pheno.cor), 'V2']
rownames(star) <- rownames(gbm.pheno.cor) 

for (i in 1:nrow(gbm.pheno.cor)) {
  if (gbm.plot[which(gbm.plot$GBM == rownames(gbm.pheno.cor)[i]), 1] > 0) {
    gbm.pheno.cor[, 'Trend'] <- 'AN-enriched'
  }
  if (gbm.plot[which(gbm.plot$GBM == rownames(gbm.pheno.cor)[i]), 1] < 0) {
    gbm.pheno.cor[, 'Trend'] <- 'NC-enriched'
  }
}


## Visualization
mycol <- colorRamp2(c(-0.3, 0),
                    c("#02CBFB", "white"))
ht.gbm.pheno.cor <- Heatmap(as.matrix(gbm.pheno.cor[, 1:5]),
                            col = mycol,
                            # right_annotation = rowAnnotation(
                            #   Trend = gbm.pheno.cor$Trend,
                            #   col = list(Trend = c("AN-enriched" = "#be8738", "NC-enriched"="#5569a5")
                            #   ),
                            #   na_col = "white",
                            #   gp = gpar(col = "black")
                            # ),
                            rect_gp = gpar(col= "black", lwd =0.5),
                            name = "Beta coefficient",
                            row_labels = rownames(gbm.pheno.cor),
                            column_labels = colnames(gbm.pheno.cor)[1:5],
                            show_row_names = T,
                            show_column_names = T,
                            cluster_columns = T,
                            show_column_dend = F,
                            cluster_column_slices = FALSE, 
                            #column_split = c(rep('AN', 10), rep('NC', 10)),
                            column_gap = unit(5, 'mm'),
                            heatmap_width = unit(0.5, 'npc'),
                            row_names_max_width = unit(25, "cm"),
                            heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(4, 'cm'), title_position = "topcenter"),
                            border = "black",
                            cluster_rows = T,
                            show_row_dend = F,
                            show_heatmap_legend = T,
                            row_names_gp = gpar(fontsize = 14),
                            column_names_gp = gpar(fontsize = 14),
                            cell_fun = function(j, i, x, y, width, height, fill){
                              grid.text(sprintf("%s", star[i, j]), x, y, gp = gpar(fontsize = 16))
                            })

pdf("gbm.pheno.cor.pdf", width = 8, height = 5)
draw(ht.gbm.pheno.cor, heatmap_legend_side = "top")
dev.off()



# Figure 4a - bar plot of amount of identified vsdv and dsgv in the combined cohort
## Load and clean SV data
vsgv <- read.csv("vsgv.csv", header = T, check.names = F, row.names = 1)
dsgv <- read.csv("dsgv.csv", header = T, check.names = F, row.names = 1)

changeSVname<-function(SVrawid){
  testname     <- SVrawid
  species_name <- as.character(taxonomy$species[match(str_replace_all(testname, "\\..*",""), taxonomy$NCBI.taxonomy.ID)])
  region       <- str_replace_all(testname, ".*\\:","") 
  region_list  <- str_split(region,";") 
  
  region_suf   <- NULL
  i<-1
  for (i in c(1:length(region_list))){
    if(length(region_list[[i]]) <= 2){
      region_suf[i] <- paste(":",region[i],sep = "")
    }else{
      region_suf[i] <- paste(":",region_list[[i]][1]," and ",length(region_list[[i]])-1," segments", sep = "")
    }
    i <- i+1
  }
  paste(species_name,region_suf,sep = "")
}

taxonomy <- read.csv("ICRAdb_annotated.csv", header = T, check.names = F, sep = ';')
colnames(vsgv) <- changeSVname(colnames(vsgv))
colnames(dsgv) <- changeSVname(colnames(dsgv))

## Split dataset
an.vsgv <- vsgv[an,]
an.dsgv <- dsgv[an,]

## Get name conversion table
## Name conversion
organism<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  .[!duplicated(.)]
Short_name<- organism %>% 
  str_replace_all('\\[','') %>%
  str_replace_all('\\]', '') %>%
  str_replace_all(' cf\\.','')

Short_name[grep(' sp\\.', organism, invert = F)] <- Short_name[grep(' sp\\.', organism, invert = F)] %>%
  str_replace_all('sp\\..*','sp')

Fst_letter<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_replace_all(' .*','') %>%
  str_sub(start = 1,end = 1)

Spe_name<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_extract_all(' .*') %>%
  str_replace_all('^ ', '') %>%
  str_replace_all(' .*', '')

Short_name[grep(' sp\\.', organism, invert = T)] <-paste(Fst_letter,'.', Spe_name, sep = '')



## Get SV number per species in the whole cohort
species_dsgv_n<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_dsgv_n)<-c("Species","Deletion SVs number")
species_vsgv_n<-str_replace_all(colnames(vsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_vsgv_n)<-c("Species","Variable SVs number")

species_sgv_n<-full_join(species_dsgv_n, species_vsgv_n, by = "Species")
species_sgv_n[is.na(species_sgv_n)]<-0

species_sgv_n$'SVs_number' <- apply(species_sgv_n[, c(2,3)], 1, sum)
species_sgv_n_order<- species_sgv_n$Species[order(species_sgv_n$SVs_number, decreasing = T)]

SV_plot <- melt(species_sgv_n, id.vars = 'Species') %>% as.data.frame()
SV_plot <- SV_plot[1:112,]

## Visualization
p_sv_whole <- ggplot(SV_plot, aes(fill=variable, y=value, x=Species, label = value)) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(size = 4, position = position_stack(vjust = 0.5), color = 'white')+
  xlab(NULL)+
  ylab("Number of SVs")+
  scale_x_discrete(limits = species_sgv_n_order)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                    labels = c("Deletion SVs              ", "Variable SVs"),
                    values = c("#0069a6", "#d14e1a"))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(0,1,1,3),"cm"),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        axis.title.y = element_text(size = 16), 
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size = 26),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

pdf("sv_whole.pdf", width = 18, height = 10)
print(p_sv_whole)
dev.off()


# Figure 4b - pie plot of compistion of SVs in the combined cohort
sv_n<-data.frame(items = rep("SVs number", 2),
                 categories = c("Deletion SV", "Variable SV"),
                 value = c(5056, 2423))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

p_pie <- ggplot(sv_n, aes(x="", y=value, fill=categories))+
  geom_bar(width = 10, stat = "identity", color = 'white') +
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values=c( "#0069a6", "#d14e1a"))
p_pie
ggsave(p_pie, filename = "piechart.pdf", width = 4, height = 4)

# Figure 4c - beta diversity of SVs in the combined cohort
## Get distance matrix
## All samples
sgv<-cbind(vsgv, dsgv)
shared_sv_dis<-function(inDf){
  #inDf<-all_sv
  inList <- split(as.matrix(inDf), seq(nrow(inDf)))
  
  shared_n_func<-function(inVec,Vec_i){
    #inVec<-inList[[3]]
    Vec_i<-Vec_i
    insvdf<-data.frame(inVec,Vec_i) %>% na.omit
    sv_dis<- vegdist(t(insvdf), method = "canberra")
    #length(na.omit(Vec_i+inVec))
    return(sv_dis)
  }
  marker_n_mat<-matrix(NA, nrow = nrow(inDf), ncol = nrow(inDf))
  for (i in 1:nrow(inDf)) { #nrow(inDf)
    #i<-2
    Vec_i<-inList[[i]]
    shared_n_i<-sapply(inList, shared_n_func,Vec_i = Vec_i)
    marker_n_mat[i,]<-shared_n_i
  }
  rownames(marker_n_mat)<-rownames(inDf)
  colnames(marker_n_mat)<-rownames(inDf)
  
  return(marker_n_mat)
}

## Remove samples with NA in distance matrix
dist_rmna<-function(inDist){
  while(sort(colSums(is.na(inDist)),decreasing = T)[1] > 0){
    rmid<-names(sort(colSums(is.na(inDist)),decreasing = T)[1])
    inDist<-inDist[-match(rmid, rownames(inDist)),-match(rmid, colnames(inDist))]
  }
  return(inDist)
}
all_shared_sv_dis<-shared_sv_dis(sgv)


## Beta diversity box plot
m.sv <- melt(all_shared_sv_dis)[melt(upper.tri(all_shared_sv_dis))$value,]
names(m.sv) <- c("c1", "c2", "distance")
an <- meta %>% filter(Category == 'AN') %>% rownames()
nc <- meta %>% filter(Category == 'NC') %>% rownames()

an.row.sv <- c()
for (i in 1:nrow(m.sv)) {
  if (m.sv[i, 1] %in% an & m.sv[i,2] %in% an) {
    an.row.sv <- c(an.row.sv, i)
  }
}

nc.row.sv <- c()
for (i in 1:nrow(m.sv)) {
  if (m.sv[i, 1] %in% nc & m.sv[i,2] %in% nc) {
    nc.row.sv <- c(nc.row.sv, i)
  }
}

wilcox.test(m.sv[nc.row.sv, 3], m.sv[an.row.sv, 3])

mean(m.sv[nc.row.sv, 3])
mean(m.sv[an.row.sv, 3], na.rm = T)


m.sv[nc.row.sv, 'Group'] <- 'NC'
m.sv[an.row.sv, 'Group'] <- 'AN'
m.sv.plot <- na.omit(m.sv)
m.sv.plot$Group <- factor(m.sv.plot$Group, levels = c('NC', 'AN'))

## Visualization
sv.beta.box <- ggplot(m.sv.plot, aes(Group, distance, fill = Group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c( "#5569a5", "#be8738"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = 'black'),
        axis.text.x = element_text(size = 16, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.position = 'none',
        axis.line = element_line(color = 'black'))+
  ylab('Canberra dissimilarity')+
  annotate("text", x=1.5, y=0.8, label= "P < 2.2e-16") 

pdf('sv.beta.nc.an.box.pdf', height = 5, width = 3)
print(sv.beta.box)
dev.off()

# Figure 4d - chord diagram of the association between eating disorder scores and bacterial SVs
## Load data
vsgv <- read.csv("vsgv.csv", header = T, check.names = F, row.names = 1)
vsgv.an.filter <- filter_pre(vsgv[an, ], 0.5, 'row')
dsgv <- read.csv("dsgv.csv", header = T, check.names = F, row.names = 1)
dsgv.an.filter <- filter_pre(dsgv[an, ], 0.5, 'row')

taxonomy <- read.csv("ICRAdb_annotated.csv", header = T, check.names = F, sep = ';')
colnames(vsgv.an.filter) <- changeSVname(colnames(vsgv.an.filter))
colnames(dsgv.an.filter) <- changeSVname(colnames(dsgv.an.filter))

sgv.4ED <- cbind(vsgv.an.filter, dsgv.an.filter)

## Execute linear regression adjusting for confounders
sgv_ED_cor <- lm_btw_mats(sgv.4ED, meta[an, EDI], meta[an, ], cov, direction = c(1,1), y_mat = 1)

sgv.ED.cor <- matrix(0, nrow = nrow(sgv_ED_cor$beta), ncol = ncol(sgv_ED_cor$beta)) %>% data.frame
rownames(sgv.ED.cor) <- rownames(sgv_ED_cor$beta)
colnames(sgv.ED.cor) <- colnames(sgv_ED_cor$beta)
for (i in 1:nrow(sgv.ED.cor)) {
  for (j in 1:ncol(sgv.ED.cor)) {
    if (sgv_ED_cor$p[i, j] < 0.05) {
      sgv.ED.cor[i,j] <- sgv_ED_cor$beta[i,j]
    }
  }
}

stars= matrix("", nrow = nrow(sgv.ED.cor), ncol = ncol(sgv.ED.cor))
rownames(stars)=rownames(sgv.ED.cor)
colnames(stars)=colnames(sgv.ED.cor)
for (i in 1:ncol(stars)) {
  for (j in 1:nrow(stars)) {
    if (sgv_ED_cor$p [j,i]<0.05) {
      stars[j,i]="+"
    }
    if (sgv_ED_cor$p [j,i]<0.01){
      stars[j,i]="++"
    }
    if (sgv_ED_cor$p [j,i]<0.001) {
      stars[j,i]="+++"
    }
  }
}

sgv.ED.cor <- t(sgv.ED.cor) %>% as.data.frame()
## Classfication the significant sgvs based on the species
for (i in 1:nrow(sgv.ED.cor)) {
  sgv.ED.cor[i, 'Species'] <- paste0(strsplit(rownames(sgv.ED.cor), '.', fixed = T)[[i]][1], ' ',strsplit(rownames(sgv.ED.cor), '.', fixed = T)[[i]][2])
}

sgv.ED.circos <- matrix(0, nrow = nrow(sgv.ED.cor), ncol = ncol(sgv.ED.cor)) %>% data.frame
colnames(sgv.ED.circos) <- colnames(sgv.ED.cor)
rownames(sgv.ED.circos) <- rownames(sgv.ED.cor)

for (i in 1:nrow(sgv.ED.cor)) {
  for (j in 1:ncol(sgv.ED.cor)) {
    if (sgv.ED.cor[i, j] > 0 | sgv.ED.cor[i, j] < 0 ) {
      sgv.ED.circos [i,j] <- 1
    }
  }
}

for (i in 1:nrow(sgv.ED.circos)) {
  sgv.ED.circos[i, 'Species'] <- paste0(strsplit(rownames(sgv.ED.circos), '.', fixed = T)[[i]][1], ' ',strsplit(rownames(sgv.ED.circos), '.', fixed = T)[[i]][2])
}
## Group by species
sgv.ED.circos.species <- sgv.ED.circos %>% group_by(Species) %>% summarize_each(funs(sum)) %>% t %>% data.frame
colnames(sgv.ED.circos.species) <- sgv.ED.circos.species[1,]
sgv.ED.circos.species <- sgv.ED.circos.species[-1,]
colnames(sgv.ED.circos.species)[14] <- '[Eubacterium] rectale'
sgv.ED.circos.species.num <- apply(sgv.ED.circos.species, 2, as.numeric)
rownames(sgv.ED.circos.species.num) <- rownames(sgv.ED.circos.species)

species.order <- colSums(sgv.ED.circos.species.num) %>% data.frame
colnames(species.order) <- 'Freq'
species.order <- arrange(species.order, desc(Freq))
species.order.str <- rownames(species.order)

ED.order <- rowSums(sgv.ED.circos.species.num) %>% data.frame
colnames(ED.order) <- 'Freq'
ED.order <- arrange(ED.order, desc(Freq))
ED.order.str <- rownames(ED.order)

## Visualization
set.seed(1)
grid.col <- c(wes_palette("Darjeeling1", length(species.order.str), type = "continuous"),
              rep('grey',length(ED.order.str)))
pdf("circus.plot.sv.ED.pdf", width = 10, height = 10)
circos.par(start.degree = 90, "clock.wise" = T)
chordDiagram(sgv.ED.circos.species.num[, species.order.str], 
             grid.col = grid.col, 
             annotationTrack = "grid",
             order = c(species.order.str, rev(ED.order.str)),
             column.col = grid.col,
             preAllocateTracks = list(track.height = 0.5))
# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA)
dev.off()
circos.clear()


# Figure 4e - heatmap showing the association between Bacteroides uniformis SVs and eating disorder scores
## Extract correlation matrix between Bacterioides uniformis and eating disorder socres
cor.BU <- sgv.ED.cor[grep("Bacteroides.uniformis.", rownames(sgv.ED.cor)),]
cor.BU <- cor.BU[,-13]
stars <- t(stars) %>% as.data.frame()
stars.BU <- stars[grep("Bacteroides.uniformis.", rownames(sgv.ED.cor)),] %>% data.frame

tmp <- cor.BU == 0
issig=rowSums(tmp, na.rm = T)
rsig=na.omit(names(issig[issig<12]))
issig=colSums(tmp, na.rm = T)
csig=na.omit(names(issig[issig<63]))

cor.BU.plot <- cor.BU[rsig, csig]
stars.BU.plot <- stars.BU[rsig, csig]

## Visualization
col_fun1 <- colorRamp2(c(-0.6, 0, 0.5),
                       c("#02CBFB", "white","#F20C34"))
ht.cor.BU <- Heatmap(as.matrix(cor.BU.plot),
                     col = col_fun1, 
                     rect_gp = gpar(col= "black", lwd =0.5),
                     name = "Beta coefficient",
                     #row_labels = paste(tax_meta[rownames(asv.contrast), 'Genus'], tax_meta[rownames(asv.contrast), 'Species'], sep = ' '),
                     show_row_names = T,
                     show_column_names = F,
                     cluster_columns = F,
                     cluster_column_slices = FALSE, 
                     #column_split = c(rep('AN', 10), rep('NC', 10)),
                     column_gap = unit(5, 'mm'),
                     heatmap_width = unit(0.5, 'npc'),
                     row_names_max_width = unit(25, "cm"),
                     heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(5, 'cm'), title_position = "topcenter"),
                     border = "black",
                     cluster_rows = F,
                     show_heatmap_legend = T,
                     row_names_gp = gpar(fontsize = 10),
                     column_names_gp = gpar(fontsize = 12),
                     cell_fun = function(j, i, x, y, width, height, fill){
                       grid.text(sprintf("%s", stars.BU.plot[i, j]), x, y, gp = gpar(fontsize = 10))
                     })

pdf("/Users/fjw536/Desktop/Anorexia/AN_DA/cor.BU.pdf", width = 6, height = 6)
draw(ht.cor.BU, heatmap_legend_side = "top")
dev.off()


# Figure 4g - Boxplot for individuals with and without the deletion region
box <- cbind(sgv$Bacteroides.uniformis.3215_3222, meta$EDI3_B, meta$EDI3_A) %>% data.frame
colnames(box) <- c('Bacteroides uniformis 3215_3222','Bulimia', 'Ascetism')
rownames(box) <- rownames(sgv)
box <- na.omit(box)
box$`Bacteroides uniformis 3215_3222` <- as.factor(box$`Bacteroides uniformis 3215_3222`)
box$'ID' <- rownames(box)

box.plot <- melt(box, value.name = 'ID') %>% data.frame

## Visualization
p <- ggplot(box.plot,aes(x = variable, y = ID.1, fill = Bacteroides.uniformis.3215_3222)) +
  # boxplot layer
  geom_boxplot(width = .3,show.legend = T,
               alpha = 0.5,
               position = position_dodge(0.4),
               outlier.color = 'grey50') +
  # geom_point(position = position_jitterdodge(0.1))+
  # geom_violin(position = position_dodge(0.4),alpha = 0.5,
  #             width = 0.7,trim = T,
  #             color = NA) +
  # change theme
  theme(axis.text = element_text(color = 'black', size = 14),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 15),
        axis.line = element_line(size = 0.5 ,linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        legend.position = 'none') +
  labs(x = NULL, y = "EDI3-score (arbitrary unit)", fill = NULL)+
  # set the color scheme
  scale_fill_manual(values = c('0'='#398AB9','1'='red'),
                    name = '') +
  # add significance
  stat_compare_means(aes(group=Bacteroides.uniformis.3215_3222),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif",
                     label.y = 22,size = 10) +
  ylim(0,35)+
  annotate("text", x=0.9, y=15, label= " n = 28")+
  annotate("text", x=1.1, y=20, label= " n = 49")+t
annotate("text", x=1.9, y=20, label= " n = 28")+
  annotate("text", x=2.1, y=28.5, label= " n = 49")
p
ggsave(p, filename = 'Contrasted_eating_score_boxplot.pdf', width = 4, height = 3.5)


# Figure 5a - PCA plot of metabolome between AN and NC
MSP[is.na(MSP)] <- 0
df <- apply(MSP, 2, qtrans)
pca_res <- prcomp(df, scale. = F)
sample_site <- pca_res$x[, 1:2] %>% data.frame

meta$Category <- as.factor(meta$Category)
pdf("Metabolome.PCA.pdf", width = 5, height = 5);
plot(sample_site$PC1,sample_site$PC2, 
     xlab=paste('PC1 (16.36%)'),
     ylab=paste('PC2 (9.5%)'),cex=0.01);
s.class(sample_site, meta$Category, col = c("#be8738", "#5569a5"), grid=F, addaxes=F,axesell =F, cellipse = 0, cpoint = 2, add.plot = T);
dev.off()


# Figure 5b - contrasted metabolites between AN and NC after drug deconfounding
mets <- read.delim('Datasets_final/Data/Metabolome.txt', header = T, row.names = 1)
all(rownames(mets)==rownames(meta.sub))

## Perform drug deconfounding
all(order(rownames(mets)) == order(rownames(meta.sub)))

mets_drugdeconfound_output <- MetaDeconfound(featureMat = mets,
                                             metaMat = meta.sub)

mets.tmp <- mets_drugdeconfound_output$status %>% as.data.frame()
mets.sig <-  mets.tmp %>% filter(Category == 'OK_nc' | Category == 'OK_sd' | Category == 'AD') %>% rownames()
mets.cliff <- matrix(NA, nrow = length(mets.sig), ncol = 6) %>% as.data.frame()
rownames(mets.cliff) <- mets.sig
colnames(mets.cliff) <- c('q_value', 'cliff', 'Prevalence_all', 'Prevalence_NC', 'Prevalence_AN', 'mets')


mets.cliff$q_value <- formatC(mets_drugdeconfound_output$Qs[mets.sig, 'Category'], format = "e", digits = 1)

for (i in mets.sig) {
  mets.cliff[i, 'cliff'] <- cliff.delta(mets[an, i],mets[nc, i], return.dm=TRUE)$estimate
  mets.cliff[i, 'Prevalence_all'] <- paste(round(sum(mets[, i] > 0)/147*100), '%', sep = '')
  mets.cliff[i, 'Prevalence_NC'] <- paste(round(sum(mets[nc, i] > 0)/70*100), '%', sep = '')
  mets.cliff[i, 'Prevalence_AN'] <- paste(round(sum(mets[an, i] > 0)/77*100), '%', sep = '')
}

mets.cliff$mets <- rownames(mets.cliff)

## Visualization
p_mets_cliff <- diverging_lollipop_chart(mets.cliff, 
                                         mets, 
                                         cliff, 
                                         lollipop_colors = c("#be8738", "#5569a5"), 
                                         line_size = 1, 
                                         point_size = 2, 
                                         text_size = 13, 
                                         text_color = c("#be8738", "#5569a5"))+
  labs(x = 'Metabolites', y = "Cliff's Delta")+ theme_light()+
  theme(axis.title = element_text(size = 18), 
        panel.grid = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ ylim(-0.5, 0.5)
p_mets_cliff

ggsave(p_mets_cliff, filename = "Contrast.mets.pdf", height = 6, width = 5)

# Figure 5d - sankey plot visualizing the mediation outcomes where microbial features were treatments, metabolites were mediators, and eating disorder scores were outcomes
## Import cleaned dataset for sankey plot
sankey.ed <- read.delim("Eating_disorder_mediation_forsankey.txt",  header = T)

## Visualization
sankey_colors<-c("#0072B2", "#999999","#D55E00", "#E69F00","#009E73",  "#56B4E9",  "#F0E442", "#CC79A7")

p.sankey.ED <- ggplot(data = sankey.ed,
                      aes(axis1 = sankey.ed$Microbial.features, 
                          axis2 = sankey.ed$Metabolites, 
                          axis3 = sankey.ed$EDI.3.scores, 
                          y = sankey.ed$Freq)) +
  scale_x_discrete(limits = c("Microbial features", "Metabolites", "EDI-3 scores")) +
  geom_alluvium(aes(fill = sankey.ed$EDI.3.scores, alpha =1),
                curve_type = "arctangent") +
  scale_color_manual(values = sankey_colors)+
  scale_fill_manual(values = sankey_colors)+
  geom_stratum(alpha = 0,color = adjustcolor( "white", alpha.f = 1),size=1.2, fill = 'white') + 
  geom_text(stat = "stratum",cex=6, aes(label = after_stat(stratum))) +
  theme_void()+
  theme(legend.position="none",
        axis.text = element_text(size = 24),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid=element_blank())
p.sankey.ED
ggsave(p.sankey.ED, filename = "sankey.ED.pdf", width = 20, height = 6)



# Extended Data Figures
# Extended Figure 2d
final.dis <- as.matrix(vegdist(MSP.genus.final, method = "canberra"))
final.dis.mds<-cmdscale(final.dis, k=5, eig = T)
msp_pcoa <- data.frame(final.dis.mds$points)

## Beta diversity of MSP matrix - box plot
m2 <- melt(final.dis)[melt(upper.tri(final.dis))$value,]
names(m2) <- c("c1", "c2", "distance")
an <- meta %>% filter(Category == 'AN') %>% rownames()
nc <- meta %>% filter(Category == 'NC') %>% rownames()

an.row <- c()
for (i in 1:nrow(m2)) {
  if (m2[i, 1] %in% an & m2[i,2] %in% an) {
    an.row <- c(an.row, i)
  }
}

nc.row <- c()
for (i in 1:nrow(m2)) {
  if (m2[i, 1] %in% nc & m2[i,2] %in% nc) {
    nc.row <- c(nc.row, i)
  }
}

wilcox.test(m2[nc.row, 3], m2[an.row, 3])
m2[nc.row, 'Group'] <- 'NC'
m2[an.row, 'Group'] <- 'AN'
m2.nona <- na.omit(m2)
m2.nona$Group <- factor(m2.nona$Group, levels = c('NC', 'AN'))
m2 <- m2[, -ncol(m2)]

## Visualization
msp.genus.beta.box <- ggplot(m2.nona, aes(Group, distance, fill = Group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c( "#5569a5", "#be8738"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = 'black'),
        axis.text.x = element_text(size = 16, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.position = 'none',
        axis.line = element_line(color = 'black'))+
  ylab('Canberra dissimilarity')+
  annotate("text", x=1.5, y=1, label= "P < 2.2e-16") 

pdf('msp.genus.beta.nc.an.box.pdf', height = 4, width = 5)
print(msp.genus.beta.box)
dev.off()

# Extended Data Figure 4 - heatmap between bacterial genera and eating disorder scores after drug deconfounding
MSP.t <- t(MSP.final) %>% data.frame
MSP.t$'Genus' <- MSP.anno[rownames(MSP.t), 2]
MSP.genus <- MSP.t %>% group_by(Genus) %>% summarize_each(funs(sum)) %>% t %>% data.frame
colnames(MSP.genus) <- MSP.genus[1, ]
MSP.genus <- MSP.genus[-1,]
MSP.genus.new <- apply(MSP.genus, 2, as.numeric) %>% data.frame
rownames(MSP.genus.new) <- rownames(MSP.genus)
MSP.genus.new <- filter_pre(MSP.genus.new, 0.1, direction = 'row')
EDI <- grep('EDI', colnames(meta))
msp.genus_EDI <- lm_btw_mats(MSP.genus.new[an,], meta[an, EDI], meta[an,], cov, direction = c(1,1), y_mat = 1)
msp.genus_EDI.sig <- msp.genus_EDI$table %>% filter(p < 0.005 | abs(Beta) > 0.5) 
msp.genus_EDI.sig

msp.genus_EDI.cor <- reshape::cast(msp.genus_EDI.sig, Phenotype ~ Taxa, mean, value = 'Beta')
msp.genus_EDI.cor[is.na(msp.genus_EDI.cor)] <- 0
rownames(msp.genus_EDI.cor) <- msp.genus_EDI.cor[,1]
msp.genus_EDI.cor <- msp.genus_EDI.cor[, -1]
msp.genus_EDI.cor.star <- reshape::cast(msp.genus_EDI.sig, Phenotype ~ Taxa , mean, value = 'fdr.p')
msp.genus_EDI.cor.star[is.na(msp.genus_EDI.cor.star)] <- ''
rownames(msp.genus_EDI.cor.star) <- msp.genus_EDI.cor.star[,1]
msp.genus_EDI.cor.star <- msp.genus_EDI.cor.star[, -1]
msp.genus_EDI.cor.star[msp.genus_EDI.cor.star > 0] <- '+'

cliff.input <- data.frame(cbind(MSP.genus.new, meta$Category))
for (i in 1:nrow(msp.genus_EDI.cor)) {
  form <- formula(paste(rownames(msp.genus_EDI.cor)[i], "~meta$Category"))
  msp.genus_EDI.cor[i, 'Cliff.delta'] <- cliffDelta(form, data = cliff.input)
  msp.genus_EDI.cor[i, 'Abs.cliff'] <- abs(msp.genus_EDI.cor[i, 'Cliff.delta'])
}

up <- msp.genus_EDI.cor %>% filter(Cliff.delta<0) %>% rownames
down <- msp.genus_EDI.cor %>% filter(Cliff.delta>0) %>% rownames
msp.genus_EDI.cor[up, 'Trend'] <- 'AN-enriched'
msp.genus_EDI.cor[down, 'Trend'] <- 'NC-enriched'

## Visualization
col_fun1 <- colorRamp2(c(-1, 0, 1),
                       c("#02CBFB","#FFFFFF","#F20C34"))
rownames(msp.genus_EDI.cor) <- gsub(".", " ", rownames(msp.genus_EDI.cor), fixed = T)

## Rename column
colnames(msp.genus_EDI.cor) <- c('Drive for thinness',
                                 'Bulimia',
                                 'Body dissatisfaction',
                                 'Low self-esteem',
                                 'Personal alienation',
                                 # 'Interpersonal insecurity',
                                 'Interpersonal alienation',
                                 'Interoceptive deficits',
                                 'Emotional disregulation',
                                 'Perfectionism',
                                 'Ascetism',
                                 'Maturity fears',
                                 'Cliff.delta',
                                 'Abs.cliff',
                                 'Trend')

ht.msp.genus.EDI.cor <- Heatmap(as.matrix(msp.genus_EDI.cor[, 1:11]), 
                                col = col_fun1, 
                                # right_annotation = rowAnnotation(
                                #   Trend = msp.genus_EDI.cor$Trend,
                                #   col = list(Trend = c("AN-enriched" = "#be8738", "NC-enriched"="#5569a5")
                                #   ),
                                #   na_col = "white",
                                #   gp = gpar(col = "black")
                                # ),
                                rect_gp = gpar(col= "black", lwd =0.5),
                                name = "Beta coefficient",
                                row_labels = paste0(rownames(msp.genus_EDI.cor)),
                                column_labels = colnames(msp.genus_EDI.cor)[1:11],
                                show_row_names = T,
                                show_column_names = T,
                                cluster_columns = FALSE,
                                heatmap_width = unit(0.5, 'npc'),
                                row_names_max_width = unit(25, "cm"),
                                heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(4, 'cm'), title_position = "topcenter"),
                                border = "black",
                                cluster_rows = F,
                                show_heatmap_legend = T,
                                row_names_gp = gpar(fontsize = 10),
                                column_names_gp = gpar(fontsize = 12),
                                cell_fun = function(j, i, x, y, width, height, fill){
                                  grid.text(sprintf("%s", msp.genus_EDI.cor.star[i, j]), x, y, gp = gpar(fontsize = 12))
                                })

pdf("ht.msp.genus.EDI.cor.pdf", width = 6, height = 10)
draw(ht.msp.genus.EDI.cor, heatmap_legend_side = "top")
dev.off()

# Extended Data Figure 5 - Fasting plasma concentration of Caseinolytic protease B
meta$'log10_ClpB' <- log10(meta$Clbp.pg.ml)
meta$Subtype <- factor(meta$Subtype, levels = c('NC', 'AN-RS', 'AN-BP'))

compare_means(log10_ClpB~Category, data=meta)
my_2color <- c("#668eb4", "#be8738")
my_2comparisons <-c("NC", 'AN')
meta$Category <- factor(meta$Category, levels = c('NC', 'AN'))

## Visualization - ClpB between AN and NC
clpb.2.box <- ggboxplot(meta, x="Category", y="log10_ClpB", fill = "Category",palette = my_2color, lwd = 1)+
  stat_compare_means( size =5)+
  theme(axis.title = element_text(size = 22), legend.position = 'none', axis.text = element_text(size = 18)) +
  labs(x = NULL, y = "Log10(ClpB)")
ggsave(clpb.2.box, filename = 'ClpB.2box.pdf', height = 5, width = 3)

Entero <- cbind(MSP.family[, which(colnames(MSP.family) == 'Enterobacteriaceae')], meta$Category) %>% as.data.frame()
colnames(Entero) <- c('Enterobacteriaceae', 'Category')
Entero[which(Entero$Category == '1'), 2] <- 'NC'
Entero[which(Entero$Category == '2'), 2] <- 'AN'
Entero$Category <- factor(Entero$Category, levels = c('NC', 'AN'))
Entero$Enterobacteriaceae <- as.numeric(Entero$Enterobacteriaceae)

## Visualization - Enterobacteriaceae between AN and NC
Entero.2.box <- ggboxplot(Entero, x="Category", y="Enterobacteriaceae", fill = "Category",palette = my_2color, lwd = 1)+
  stat_compare_means( size =5)+
  theme(axis.title = element_text(size = 22), legend.position = 'none', axis.text = element_text(size = 18)) +
  labs(x = NULL, y = "Enterobacteriaceae")
ggsave(Entero.2.box, filename = 'Enterobacteriaceae.box.pdf', height = 5, width = 5)

## Visualization - ClpB between NC and two AN subtypes
my_color <- c("#668eb4", "#99d8cf", "#f8dda6")
my_comparisons <- list(c("NC", "AN-RS"), c("NC", "AN-BP"), c("AN-RS", "AN-BP"))
clpb.box <- ggboxplot(meta, x="Subtype", y="log10_ClpB", fill = "Subtype",palette = my_color, lwd = 1)+
  stat_compare_means(comparisons = my_comparisons, size =5)+
  theme(axis.title = element_text(size = 22), legend.position = 'none', axis.text = element_text(size = 18)) +
  labs(x = NULL, y = "Log10(ClpB)")
ggsave(clpb.box, filename = 'ClpB.box.pdf', height = 5, width = 5)

# Extended Data Figure 6 - PTR data visualization
ptr <- read.csv("PTR.result.csv", header = T, sep = ",", quote = "\t", row.names = 1, check.names = F)
meta <- read.delim("Pheno.txt", header = T, row.names = 1)

wilcox.test(rowMeans(ptr, na.rm = T)~meta$Category)
medianWithoutNA<-function(x) {
  median(x[which(!is.na(x))])
}
wilcox.test(apply(ptr, 1, medianWithoutNA)~meta$Category)
kruskal.test(apply(ptr, 1, medianWithoutNA)~meta$Subtype)

ptr$'Median' <- apply(ptr, 1, medianWithoutNA)

fit_ptr <- matrix(NA, nrow = ncol(ptr), ncol = 2)
rownames(fit_ptr) <- colnames(ptr)
colnames(fit_ptr) <- c('pvalue', 'cliffdelta')
fit_ptr <- as.data.frame(fit_ptr)
f.matrix <- data.frame(cbind(ptr, meta$Category, meta$Subtype))

compare_means(Median~meta.Subtype, data=f.matrix)
my_color <- c("#78a0c1", "#a8ded8", "#fae3b5")
f.matrix$meta.Subtype <- factor(f.matrix$meta.Subtype, levels = c ("NC", "AN-RS", "AN-BP"))
my_comparisons <- list(c("NC", "AN-RS"), c("NC", "AN-BP"), c("AN-RS", "AN-BP"))
p_median_ptr <- ggboxplot(f.matrix, x="meta.Subtype", y="Median", fill = "meta.Subtype", palette = my_color, lwd = 1)+
  stat_compare_means(comparisons = my_comparisons, size =5)+
  theme(axis.title = element_text(size = 22), legend.position = 'none', axis.text = element_text(size = 18)) +
  labs(x = NULL, y = "Median Peak-to-trough ratio")

f.matrix$meta.Category <- factor(f.matrix$meta.Category, levels = c("NC", "AN"))
p_median_ptr_category <- ggboxplot(f.matrix, x="meta.Category", y="Median", fill = "meta.Category",palette = c("#5a84a9", "#de8f81"), lwd = 1)+
  stat_compare_means(size =5)+
  theme(axis.title = element_text(size = 22), legend.position = 'none', axis.text = element_text(size = 18)) +
  labs(x = NULL, y = "Peak-to-Trough ratio")

boxplot.input[, 1:7] <- apply(boxplot.input[, 1:7], 2, as.numeric)

f.matrix$meta.Category <- as.character(f.matrix$meta.Category)

boxplot.input <- data.frame(cbind(f.matrix$Median,
                                  f.matrix$Akkermansia.muciniphila, 
                                  f.matrix$Alistipes.finegoldii, 
                                  f.matrix$Coprococcus.catus, 
                                  f.matrix$Eubacterium.siraeum,
                                  f.matrix$Odoribacter.splanchnicus,
                                  f.matrix$butyrate.producing.bacterium.SS3.4,
                                  f.matrix$meta.Category))
colnames(boxplot.input) <- c('Median***','Akkermansia muciniphila*', 'Alistipes finegoldii*', 'Coprococcus catus*', 'Eubacterium siraeum*', 
                             'Odoribacter splanchnicus*', 'butyrate-producing bacterium SS3/4*', 'Category')
boxplot <- melt(boxplot.input, id.vars = "Category")
boxplot$value <- as.numeric(boxplot$value)

boxplot$Category <- factor(boxplot$Category, levels = c('NC', 'AN'))


p.contrast.ptr <- ggboxplot(boxplot, x="variable", y="value", fill = "Category",palette = c("#5a84a9", "#de8f81"), lwd = 0.5)+
  #stat_compare_means()+
  theme(axis.title = element_text(size = 16), legend.position = 'none', axis.text = element_text(size = 14)) +
  labs(x = NULL, y = "Peak-to-Trough ratio")+
  #facet_wrap(~variable, scale="free", nrow = 1)+
  theme(strip.text.x = element_blank(),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        axis.text.y = element_text(face = "italic"))+
  coord_flip()


ggsave(p.contrast.ptr, filename = "Contrast_ptr.pdf", width = 10, height = 5)
pdf("Contrast_ptr.pdf")
print(p.contrast.ptr)
dev.off()



# Extended Data Figures 7a, 7b. 7c - A. putredinis and glucose metabolism
cov.1 <- c('Age', 'Smoking', 'SSRI', 'Antipsychotics', 'Benzodiazepines')
cor.res.AP.glucose <- lm_btw_mats(sgv.an[, c('Alistipes putredinis:936_937', 
                                             'Lachnospiraceae bacterium 3_1_46FAA:1211_1214 and 2 segments')], 
                                  meta[an, c('HOMA.IR', 'P.insulin', 'P.glucose')], 
                                  meta[an,], cov.1, direction = c(1,1), y_mat = 1)
cor.res.AP.glucose$table


# Extended Data Figure 8a - PCA plot of metabolome between AN subtypes and NC
meta$Subtype <-  factor(meta$Subtype, levels = c ("NC", "AN-RS", "AN-BP"))
pdf("Metabolome.PCA.subtype.pdf", width = 5, height = 5);
plot(sample_site$PC1,sample_site$PC2, 
     xlab=paste('PC1 (16.36%)'),
     ylab=paste('PC2 (9.5%)'),cex=0.01);
s.class(sample_site, meta$Subtype, col = c("#5569a5", "#f88d7b", "#79b4a7"), grid=F, addaxes=F, axesell =F, cellipse = 0, cpoint = 2, add.plot = T);
dev.off()


#  Extended Data Figure 8b and 8c - Venn diagrams for the visualization of numbers of mediations
## Visualization of Figure 8B
grid.newpage()
venn.plot.micro.EDI <- draw.pairwise.venn(area1=24, area2=22,cross.area=11,
                                          category=c("Direction 1","Direction 2"),
                                          fill=c("Red","Yellow"),
                                          fontfamily = rep('Arial', 3),
                                          cat.fontfamily = 'Arial',
                                          cex = 2,
                                          cat.pos = 0,
                                          cat.cex = 1.5)
## Visualization of Figure 8c
grid.newpage()
draw.pairwise.venn(area1=94, area2=95,cross.area=77,
                   category=c("Direction 1","Direction 2"),
                   fill=c("Red","Yellow"),
                   fontfamily = rep('Arial', 3),
                   cat.fontfamily = 'Arial',
                   cex = 2,
                   cat.pos = c(10, 170),
                   cat.cex = 1.5)


# Extended Data Figure 8d - sankey plot visualizing the mediation analysis outcomes in direction 1 (microbial features --> metabolites --> phenotypes)
pheno.sankey <- read.delim('Phenotype_mediation_for_sankey.txt', header = T, row.names = 1)

sankey_colors<-c("#0072B2", "#999999","#D55E00", "#E69F00","#009E73",  "#56B4E9",  "#F0E442", "#CC79A7")

phenotype.sankey <- ggplot(data = pheno.sankey,
                           aes(axis1 = pheno.sankey$Microbial.features, 
                               axis2 = pheno.sankey$Metabolites, 
                               axis3 = pheno.sankey$Phenotypes, 
                               y = pheno.sankey$Freq)) +
  scale_x_discrete(limits = c("Microbial features", "Metabolites", "Phenotypes")) +
  geom_alluvium(aes(fill = pheno.sankey$Phenotypes, alpha =1),
                curve_type = "arctangent") +
  scale_color_manual(values = sankey_colors)+
  scale_fill_manual(values = sankey_colors)+
  geom_stratum(alpha = 0,color = adjustcolor( "white", alpha.f = 1),size=1.2, fill = 'white') + 
  geom_text(stat = "stratum",cex=8, aes(label = after_stat(stratum))) +
  theme_void()+
  theme(legend.position="none",
        axis.text = element_text(size = 24),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid=element_blank())
phenotype.sankey
ggsave(phenotype.sankey, filename = "sankey.pheno.pdf", width = 20, height = 10)

# Figure 10d - heatmap showing the correlation between transferred ASVs and gene expression in adipose tissue and hypothalamus of GF recipients
gene <- read.delim("/Users/fjw536/Desktop/Anorexia/16S/AN_FMT_gene.txt", header = T, row.names = 1)
asv.cor.ID <- read.delim('/Users/fjw536/Downloads/AN_FMT_16s/contrast_ASV.anno.update.txt', header = T, row.names = 1) %>% rownames()
asv.rel <- apply(ASV_table, 2, function(x){x/sum(x)})
asv.cor <- t(asv.rel[asv.cor.ID,]) %>% as.data.frame

gene_asv_cor <- cor_mat(asv.cor, gene, cor = 'pearson', direction = c(1,1))

gene.asv.cor <- matrix(0, nrow = nrow(gene_asv_cor$rho), ncol = ncol(gene_asv_cor$rho)) %>% data.frame
rownames(gene.asv.cor) <- rownames(gene_asv_cor$rho)
colnames(gene.asv.cor) <- colnames(gene_asv_cor$rho)
for (i in 1:nrow(gene.asv.cor)) {
  for (j in 1:ncol(gene.asv.cor)) {
    if (gene_asv_cor$BH[i, j] < 0.01) {
      gene.asv.cor[i,j] <- gene_asv_cor$rho[i,j]
    }
  }
}

stars= matrix("", nrow = nrow(gene.asv.cor), ncol = ncol(gene.asv.cor))
rownames(stars)=rownames(gene.asv.cor)
colnames(stars)=colnames(gene.asv.cor)
for (i in 1:ncol(stars)) {
  for (j in 1:nrow(stars)) {
    if (gene_asv_cor$BH [j,i]<0.01) {
      stars[j,i]="+"
    }
    if (gene_asv_cor$BH [j,i]<0.001){
      stars[j,i]="++"
    }
    # if (gene_asv_cor$BH [j,i]<0.0001) {
    #   stars[j,i]="+++"
    # }
  }
}

## Visualization
col_fun1 <- colorRamp2(c(-1, 0, 1),
                       c("#02CBFB", "white","#F20C34"))

## Add annotations for the ASVst
tax_anno <- tax_meta[rownames(gene.asv.cor), ] %>% as.data.frame()
tax_anno$'Annotation' <- paste(rownames(tax_anno), ': p_', tax_anno$Phylum, '|f_', tax_anno$Family, '|g_', tax_anno$Genus, sep = '')

gene.asv.cor.plot <- gene.asv.cor
rownames(gene.asv.cor.plot) <- tax_anno$'Annotation'
ht.bac.gene.cor <- Heatmap(as.matrix(gene.asv.cor.plot),
                           col = col_fun1, 
                           rect_gp = gpar(col= "black", lwd =0.5),
                           name = "Pearson coefficient",
                           #row_labels = paste(tax_meta[rownames(asv.contrast), 'Genus'], tax_meta[rownames(asv.contrast), 'Species'], sep = ' '),
                           show_row_names = T,
                           show_column_names = T,
                           cluster_columns = F,
                           show_column_dend = F,
                           cluster_column_slices = FALSE, 
                           #column_split = c(rep('AN', 10), rep('NC', 10)),
                           column_gap = unit(5, 'mm'),
                           heatmap_width = unit(0.5, 'npc'),
                           row_names_max_width = unit(25, "cm"),
                           heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(4, 'cm'), title_position = "topcenter"),
                           border = "black",
                           cluster_rows = T,
                           show_row_dend = F,
                           show_heatmap_legend = T,
                           row_names_gp = gpar(fontsize = 16),
                           column_names_gp = gpar(fontsize = 16),
                           cell_fun = function(j, i, x, y, width, height, fill){
                             grid.text(sprintf("%s", stars[i, j]), x, y, gp = gpar(fontsize = 18))
                           })
pdf("/Users/fjw536/Desktop/Anorexia/16S/gene.asv.cor.pdf", width = 9, height = 6)
draw(ht.bac.gene.cor, heatmap_legend_side = "top")
dev.off()



# Supplementary Figure 1 - box plot of beta diversity between inpatient and outpatient AN
m2[IN.row, 'Group'] <- 'in'
m2[out.row, 'Group'] <- 'out'
m2.nona <- na.omit(m2)
m2.nona$Group <- factor(m2.nona$Group, levels = c('in', 'out'))
msp.beta.in.out.box <- ggplot(m2.nona, aes(Group, distance, fill = Group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c( "#4ac2e6", "#e4728d"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = 'black'),
        axis.text.x = element_text(size = 16, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.position = 'none',
        axis.line = element_line(color = 'black'))+
  ylab('Canberra dissimilarity')+
  annotate("text", x=1.5, y=1, label= "P = 8.2e-11") 

pdf('msp.beta.in.out.box.pdf', height = 5, width = 3)
print(msp.beta.in.out.box)
dev.off()

# Supplementary Figure 2a - heatmap showing the association between bacterial taxa at genus level and bioclinical variables in the combined dataset
## Genus
MSP.t <- t(MSP.final) %>% data.frame
MSP.t$'Genus' <- MSP.anno[rownames(MSP.t), 2]
MSP.genus <- MSP.t %>% group_by(Genus) %>% summarize_each(funs(sum)) %>% t %>% data.frame
colnames(MSP.genus) <- MSP.genus[1, ]
MSP.genus <- MSP.genus[-1,]
MSP.genus.new <- apply(MSP.genus, 2, as.numeric) %>% data.frame
rownames(MSP.genus.new) <- rownames(MSP.genus)
MSP.genus.final <- filter_pre(MSP.genus.new, 0.1, direction = 'row')

## Association with bioclinical variables in combined cohort
MSP.genus.bio_cor <- lm_btw_mats(MSP.genus.final, meta[,2:9], meta, cov.1, y_mat = 1)
MSP.genus.bio.sig <- MSP.genus.bio_cor$table %>% filter(fdr.p < 0.1)
MSP.genus.bio.an_cor <- lm_btw_mats(MSP.genus.final[an,], meta[an,2:9], meta[an,], cov.1, y_mat = 1)
MSP.genus.bio.an.sig <- MSP.genus.bio.an_cor$table %>% filter(fdr.p < 0.1)

## Visualization
msp.genus_bio.cor <- reshape::cast(MSP.genus.bio.sig, Phenotype ~ Taxa, mean, value = 'Beta')
msp.genus_bio.cor[is.na(msp.genus_bio.cor)] <- 0
rownames(msp.genus_bio.cor) <- msp.genus_bio.cor[,1]
msp.genus_bio.cor <- msp.genus_bio.cor[, -1]
msp.genus_bio.cor.star <- reshape::cast(MSP.genus.bio.sig, Phenotype ~ Taxa , mean, value = 'fdr.p')
msp.genus_bio.cor.star[is.na(msp.genus_bio.cor.star)] <- ''
rownames(msp.genus_bio.cor.star) <- msp.genus_bio.cor.star[,1]
msp.genus_bio.cor.star <- msp.genus_bio.cor.star[, -1]
msp.genus_bio.cor.star[msp.genus_bio.cor.star > 0] <- '+'

col_fun1 <- colorRamp2(c(-0.5, 0, 0.5),
                       c("#02CBFB","#FFFFFF","#F20C34"))
rownames(msp.genus_bio.cor) <- gsub(".", " ", rownames(msp.genus_bio.cor), fixed = T)
colnames(msp.genus_bio.cor) <- gsub(".", "-", colnames(msp.genus_bio.cor), fixed = T)
colnames(msp.genus_bio.cor) <- gsub("_", " ", colnames(msp.genus_bio.cor), fixed = T)

## Make heatmap
ht.msp.genus.bio.cor <- Heatmap(as.matrix(msp.genus_bio.cor), 
                                col = col_fun1, 
                                # right_annotation = rowAnnotation(
                                #   Trend = msp.genus_EDI.cor$Trend,
                                #   col = list(Trend = c("AN-enriched" = "#be8738", "NC-enriched"="#5569a5")
                                #   ),
                                #   na_col = "white",
                                #   gp = gpar(col = "black")
                                # ),
                                rect_gp = gpar(col= "black", lwd =0.5),
                                name = "Beta coefficient",
                                row_labels = paste0(rownames(msp.genus_bio.cor)),
                                column_labels = colnames(msp.genus_bio.cor),
                                show_row_names = T,
                                show_column_names = T,
                                cluster_columns = FALSE,
                                heatmap_width = unit(0.5, 'npc'),
                                row_names_max_width = unit(25, "cm"),
                                heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(4, 'cm'), title_position = "topcenter"),
                                border = "black",
                                cluster_rows = F,
                                show_heatmap_legend = T,
                                row_names_gp = gpar(fontsize = 10),
                                column_names_gp = gpar(fontsize = 12),
                                cell_fun = function(j, i, x, y, width, height, fill){
                                  grid.text(sprintf("%s", msp.genus_bio.cor.star[i, j]), x, y, gp = gpar(fontsize = 12))
                                })

pdf("ht.msp.genus.bio.cor.pdf", width = 3, height = 4)
draw(ht.msp.genus.bio.cor, heatmap_legend_side = "top")
dev.off()

# Supplementary Figure 2b - heatmap showing the association between bacterial taxa at species level and bioclinical variables in the combined dataset
## Species
MSP.species.final <- filter_pre(MSP.final, 0.1, direction = 'row')
MSP.species.bio_cor <- lm_btw_mats(MSP.species.final, meta[,2:9], meta, cov.1, y_mat = 1)
MSP.species.bio.sig <- MSP.species.bio_cor$table %>% filter(fdr.p < 0.1)
MSP.species.bio.an_cor <- lm_btw_mats(MSP.species.final[an,], meta[an,2:9], meta[an,], cov.1, y_mat = 1)
MSP.species.bio.an.cor <- MSP.species.bio.an_cor$table %>% filter(fdr.p < 0.1)

msp.species_bio.cor <- reshape::cast(MSP.species.bio.sig, Phenotype ~ Taxa, mean, value = 'Beta')
msp.species_bio.cor[is.na(msp.species_bio.cor)] <- 0
rownames(msp.species_bio.cor) <- msp.species_bio.cor[,1]
msp.species_bio.cor <- msp.species_bio.cor[, -1]
msp.species_bio.cor.star <- reshape::cast(MSP.species.bio.sig, Phenotype ~ Taxa , mean, value = 'fdr.p')
msp.species_bio.cor.star[is.na(msp.species_bio.cor.star)] <- ''
rownames(msp.species_bio.cor.star) <- msp.species_bio.cor.star[,1]
msp.species_bio.cor.star <- msp.species_bio.cor.star[, -1]
msp.species_bio.cor.star[msp.species_bio.cor.star > 0] <- '+'

col_fun1 <- colorRamp2(c(-0.5, 0, 0.5),
                       c("#02CBFB","#FFFFFF","#F20C34"))
rownames(msp.species_bio.cor) <- gsub(".", " ", rownames(msp.species_bio.cor), fixed = T)
colnames(msp.species_bio.cor) <- gsub(".", "-", colnames(msp.species_bio.cor), fixed = T)
colnames(msp.species_bio.cor) <- gsub("_", " ", colnames(msp.species_bio.cor), fixed = T)
msp.species_bio.cor$'Annotation' <- paste(rownames(msp.species_bio.cor), ': ', MSP.anno[rownames(msp.species_bio.cor), 'species'], sep = '')

## Make heatmap
ht.msp.species.bio.cor <- Heatmap(as.matrix(msp.species_bio.cor[,1:7]), 
                                  col = col_fun1, 
                                  # right_annotation = rowAnnotation(
                                  #   Trend = msp.genus_EDI.cor$Trend,
                                  #   col = list(Trend = c("AN-enriched" = "#be8738", "NC-enriched"="#5569a5")
                                  #   ),
                                  #   na_col = "white",
                                  #   gp = gpar(col = "black")
                                  # ),
                                  rect_gp = gpar(col= "black", lwd =0.5),
                                  name = "Beta coefficient",
                                  row_labels = msp.species_bio.cor$Annotation,
                                  column_labels = colnames(msp.species_bio.cor)[1:7],
                                  show_row_names = T,
                                  show_column_names = T,
                                  cluster_columns = FALSE,
                                  heatmap_width = unit(0.5, 'npc'),
                                  row_names_max_width = unit(25, "cm"),
                                  heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(4, 'cm'), title_position = "topcenter"),
                                  border = "black",
                                  cluster_rows = F,
                                  show_heatmap_legend = T,
                                  row_names_gp = gpar(fontsize = 10),
                                  column_names_gp = gpar(fontsize = 12),
                                  cell_fun = function(j, i, x, y, width, height, fill){
                                    grid.text(sprintf("%s", msp.species_bio.cor.star[i, j]), x, y, gp = gpar(fontsize = 12))
                                  })

pdf("ht.msp.species.bio.cor.pdf", width = 6, height = 9)
draw(ht.msp.species.bio.cor, heatmap_legend_side = "top")
dev.off()


# Supplementary Figure 3 - heatmap showing the viral-bacterial transkingdom interactions between AN and NC groups
colnames(NC.sparcc) <- gsub(".", " ", colnames(NC.sparcc), fixed = T)
colnames(AN.sparcc) <- gsub(".", " ", colnames(AN.sparcc), fixed = T)

col_fun1 <- colorRamp2(c(-0.6, 0, 0.6),
                       c("#02CBFB","#FFFFFF","#F20C34"))

ht.nc <- Heatmap(as.matrix(NC.sparcc), 
                 col = col_fun1, 
                 #rect_gp = gpar(col= "black", lwd =0.5),
                 name = "NC",
                 row_labels = paste0(rownames(NC.sparcc), ":", MSP.anno[rownames(NC.sparcc), 1]),
                 show_row_names = F,
                 cluster_columns = FALSE,
                 heatmap_width = unit(0.5, 'npc'),
                 row_names_max_width = unit(25, "cm"),
                 heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(5, 'cm'), title_position = "topcenter", at = seq(-0.6, 0.6, 0.3)),
                 border = "black",
                 cluster_rows = F,
                 show_heatmap_legend = F,
                 row_names_gp = gpar(fontsize = 6.5),
                 column_names_gp = gpar(fontsize = 6.5))


ht.an <- Heatmap(as.matrix(AN.sparcc), 
                 col = col_fun1, 
                 #rect_gp = gpar(col= "white", lwd =1),
                 name = "SparCC Rho",
                 row_labels = paste0(rownames(AN.sparcc), ":", MSP.anno[rownames(AN.sparcc), 1]),
                 cluster_columns = FALSE,
                 heatmap_width = unit(1, 'npc'),
                 heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(5, 'cm'), title_position = "topcenter", at = seq(-0.6, 0.6, 0.3)),
                 border = "black",
                 show_heatmap_legend = T,
                 row_names_max_width = unit(25, "cm"),
                 cluster_rows = F,
                 row_names_gp = gpar(fontsize = 6.5),
                 column_names_gp = gpar(fontsize = 6.5))

ht.main <- draw(ht.nc + ht.an, heatmap_legend_side = "top")
pdf("ht.main.pdf", width = 16, height = 15)
draw(ht.main)
dev.off()

# Supplementary Figure 4 - heatmap showing the viral-bacterial transkingdom interactions between two AN subtypes
col_fun2 <- colorRamp2(c(-1, 0, 1),
                       c("#02CBFB","#FFFFFF","#F20C34"))
colnames(BP.sparcc) <- gsub(".", " ", colnames(BP.sparcc), fixed = T)
colnames(RS.sparcc) <- gsub(".", " ", colnames(RS.sparcc), fixed = T)
ht.bp <- Heatmap(as.matrix(BP.sparcc), 
                 col = col_fun2, 
                 #rect_gp = gpar(col= "black", lwd =0.5),
                 name = "SparCC Rho",
                 row_labels = paste0(rownames(BP.sparcc), ":", MSP.anno[rownames(BP.sparcc), 1]),
                 show_row_names = F,
                 cluster_columns = FALSE,
                 heatmap_width = unit(0.5, 'npc'),
                 row_names_max_width = unit(20, "cm"),
                 heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(5, 'cm'), title_position = "topcenter", at = seq(-1, 1, 0.5)),
                 border = "black",
                 cluster_rows = F,
                 show_heatmap_legend = T,
                 row_names_gp = gpar(fontsize = 6.5),
                 column_names_gp = gpar(fontsize = 6.5))

ht.rs <- Heatmap(as.matrix(RS.sparcc), 
                 col = col_fun2, 
                 #rect_gp = gpar(col= "black", lwd =0.5),
                 name = "RS",
                 row_labels = paste0(rownames(RS.sparcc), ":", MSP.anno[rownames(RS.sparcc), 1]),
                 cluster_columns = FALSE,
                 heatmap_width = unit(1, 'npc'),
                 row_names_max_width = unit(20, "cm"),
                 heatmap_legend_param = list(legend_direction = c("horizontal"), legend_width = unit(5, 'cm'), title_position = "topcenter"),
                 border = "black",
                 cluster_rows = F,
                 show_heatmap_legend = F,
                 row_names_gp = gpar(fontsize = 6.5),
                 column_names_gp = gpar(fontsize = 6.5))
ht.sub <- draw(ht.bp + ht.rs, heatmap_legend_side = "top")

pdf("ht.sub.pdf", width = 16, height = 15)
draw(ht.sub)
dev.off()


# Supplementary Figure 5 - PCoA plot of the Canberra distance showing the stratification of AN-RS from the AN-BP at species level.
## Subset the virome matrix and metadata
virus.sub <- virus[an,]
meta.sub <- meta %>% filter(Category == 'AN')
meta.sub$Subtype <- factor(meta.sub$Subtype, levels = c('AN-BP', 'AN-RS'))

final.dis <- vegdist(virus.sub, method = "canberra")
pcoa<- dudi.pco(final.dis, scan = FALSE, nf=2)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_site <- data.frame({pcoa$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
pcoa_virus_subtype <- ggplot(sample_site, aes(PCoA1, PCoA2, color = meta.sub$Subtype))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 3)+ 
  scale_color_manual(values = c("#a8ded8", "#fae3b5"))+
  stat_ellipse(type = 't', level = 0.95, linetype = 2, lwd = 1.5)+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank(),
        legend.position = 'bottom', legend.text = element_text(size = 20),
        axis.text.x=element_text(colour="black",family="Arial",size=14), 
        axis.text.y=element_text(colour="black", family="Arial",size=14,face="plain"), 
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), 
        axis.title.x=element_text(family="Arial",size = 20,face="plain"), 
        plot.title = element_text(family = "Arial", size=18, vjust = -5),
        plot.caption = element_text(size = 18))+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'), caption = "adonis p = 0.571")
library(ggExtra)
p_vir_sub <- ggMarginal(pcoa_virus_subtype, type = "boxplot", groupColour = T, groupFill = T)


# Supplementary Figure 6 - example code for plotting SV distributions in bacteriome
## TOTAL SVs per sample ----
total <- cbind(dsgv, vsgv)
a <- as.data.frame(t(total)) ; a$bacteria <- names(total)
a <- data.frame('bact'= unlist(strsplit(a$bacteria, ':'))[c(T,F)], 'regions' = unlist(strsplit(a$bacteria, ':'))[c(F,T)], a) ## divide by bacterium
a$bacteria <- NULL
a.bact <- split(a, a$bact)

smp.SVs <- smp.dels <- data.frame('samples.SVs' = unlist(lapply(a.bact, function(x){ ## number of samples with variable regions per bacterium (non 0)! | there are 147 anorexia samples
  return(147-sum(colSums(is.na(x[,3:149]))/nrow(x))) 
})))

smp.SVs <- data.frame(smp.SVs, 'regions.SVs' = unlist(lapply(a.bact, function(y)(nrow(y[rowSums(is.na(y[,3:149]))<147,]))))) ## compute number of regions detected per bacterium
db.1 <- db[db$species%in%row.names(smp.SVs),5:11] ## retrieve bacterial taxonomy in samples
smp.SVs <- data.frame(smp.SVs, db.1[match(row.names(smp.SVs), db.1$species),])
smp.SVs <- smp.SVs[order(smp.SVs$superkingdom, smp.SVs$phylum, smp.SVs$class, smp.SVs$order, smp.SVs$family, smp.SVs$genus, smp.SVs$species),] ## order dataframe by taxonomical levels
smp.SVs1 <- smp.SVs ## store original

smp.SVs1 <- melt(smp.SVs) ## prepare data.frame for heatmap plotting
textcol <- "grey40"
color.scale <- c('samples' = 'darkolivegreen2', 'regions' = 'orchid') 

smp.SVs1.1 <- smp.SVs1 %>%
  group_by(variable) %>%
  mutate(value.alpha = scale(value)) %>%
  ungroup() %>%
  arrange(variable) ## add plotting alpha scales

smp.SVs1.1$species <- factor(smp.SVs1.1$species, levels = unique(smp.SVs1.1$species)) ## keep bacteria order for the plot
p <- ggplot(smp.SVs1.1, aes(variable, species, fill=variable, alpha = value.alpha)) + geom_tile(colour="white",size=0.2)+
  geom_text(aes(label=value), alpha = 1, cex = 2) +
  # guides(fill=guide_legend(title="Cases per\n100,000 people"))+
  # labs(x="",y="",title="Incidence of Measles in the US")+
  scale_y_discrete(limits = rev(levels(smp.SVs1.1$species)), position = "right") +
  # scale_x_discrete(limits = rev(levels(species)))+
  # scale_fill_manual(values = color.scale) +
  scale_alpha(range = c(0.1, 0.9)) +  coord_fixed() + labs(x = NULL, y = NULL) +
  theme(legend.position="none",legend.direction="vertical", legend.title=element_text(colour=textcol),legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),legend.key.height=grid::unit(0.8,"cm"),legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol, vjust = 0.25, angle = 90, face = 'bold'), axis.text.y=element_text(vjust=0.2,colour=textcol, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_blank(),plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) 
p

## DELETIONS: number of deletions per bacterium ----
b <- as.data.frame(t(dsgv)) ## traspose df
b <- data.frame('bact'= unlist(strsplit(row.names(b), ':'))[c(T,F)], 'regions' = unlist(strsplit(row.names(b), ':'))[c(F,T)], b) ## divide by bacterium
b.bact <- split(b, b$bact) ## separate in list per bacterium

smp.dels <- data.frame('samples.deletions' = unlist(lapply(b.bact, function(x){ ## number of samples with variable regions per bacterium (non 0)! | there are 147 anorexia samples
  return(147-sum(colSums(is.na(x[,3:149]))/nrow(x))) 
})))

smp.dels <- data.frame(smp.dels, 'regions.deletions' = unlist(lapply(b.bact, function(y)(nrow(y[rowSums(is.na(y[,3:149]))<147,]))))) ## compute number of regions detected per bacterium
db.1 <- db[db$species%in%row.names(smp.dels),5:11] ## retrieve bacterial taxonomy in samples
smp.dels <- data.frame(smp.dels, db.1[match(row.names(smp.dels), db.1$species),])
smp.dels <- smp.dels[order(smp.dels$superkingdom, smp.dels$phylum, smp.dels$class, smp.dels$order, smp.dels$family, smp.dels$genus, smp.dels$species),] ## order dataframe by taxonomical levels
smp.dels1 <- smp.dels ## store original

smp.dels1 <- melt(smp.dels) ## prepare data.frame for heatmap plotting
textcol <- "grey40"
color.scale <- c('samples' = 'darkolivegreen2', 'regions' = 'orchid') 

smp.dels1.1 <- smp.dels1 %>%
  group_by(variable) %>%
  mutate(value.alpha = scale(value)) %>%
  ungroup() %>%
  arrange(variable) ## add plotting alpha scales

smp.dels1.1$species <- factor(smp.dels1.1$species, levels = unique(smp.dels1.1$species)) ## keep bacteria order for the plot
p <- ggplot(smp.dels1.1, aes(variable, species, fill=variable, alpha = value.alpha)) + geom_tile(colour="white",size=0.2)+
  geom_text(aes(label=value), alpha = 1, cex = 2) +
  # guides(fill=guide_legend(title="Cases per\n100,000 people"))+
  # labs(x="",y="",title="Incidence of Measles in the US")+
  scale_y_discrete(limits = rev(levels(smp.dels1.1$species)), position = "right") +
  # scale_x_discrete(limits = rev(levels(species)))+
  # scale_fill_manual(values = color.scale) +
  scale_alpha(range = c(0.1, 0.9)) +  coord_fixed() + labs(x = NULL, y = NULL) +
  theme(legend.position="none",legend.direction="vertical", legend.title=element_text(colour=textcol),legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),legend.key.height=grid::unit(0.8,"cm"),legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol, vjust = 0.25, angle = 90, face = 'bold'), axis.text.y=element_text(vjust=0.2,colour=textcol, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_blank(),plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) 

smp.dels.leg <- smp.dels1.1 %>% group_by(variable) %>% 
  summarise(min = min(value), max = max(value)) %>%
  mutate(tile.1 = 0.1, tile.3 = 0.3, tile.5 = 0.5, tile.7 = 0.7, tile.9 = 0.9) %>%
  tidyr::gather(tile, value, -variable, -min, -max)
p.legend <- ggplot(smp.dels.leg, aes(x = value, y = variable, fill = variable, alpha = value)) +
  geom_tile() +
  geom_text(aes(label = variable), x = -0.1) +
  geom_text(aes(label = min), x = 0.1) +
  geom_text(aes(label = max), x = 0.9, color = "white") +
  ggtitle("Legend") +
  coord_cartesian(xlim = c(-0.2, 1.1)) +
  scale_fill_manual(values = color.scale) +
  scale_alpha_identity() +
  theme_void() +
  theme(legend.position = "none")
p.legend
p
cowplot::plot_grid(p, p.legend, 
                   ncol = 2, align = "v", 
                   rel_heights = c(5, 1))

## VARIATIONS: number of SVS per bacterium ----
v <- as.data.frame(t(vsgv)) ## traspose data.frame
v <- data.frame('bact'= unlist(strsplit(row.names(v), ':'))[c(T,F)], 'regions' = unlist(strsplit(row.names(v), ':'))[c(F,T)], v) ## bacteria names and regions
v.bact <- split(v, v$bact) ## split by bacterium

smp.vars <- data.frame('samples.variations' = unlist(lapply(v.bact, function(x){ ## number of samples with variable regions per bacterium (non 0)
  return(147-sum(colSums(is.na(x[,3:149]))/nrow(x))) 
})))

smp.vars <- data.frame(smp.vars, 'regions.variations' = unlist(lapply(v.bact,function(y)(nrow(y[rowSums(is.na(y[,3:149]))<147,]))))) ## recount regions present in at least one sample
db.1 <- db[db$species%in%row.names(smp.vars),5:11] ## retreive bacterial taxonomy in our samples
smp.vars <- data.frame(smp.vars, db.1[match(row.names(smp.vars), db.1$species),])
smp.vars <- smp.vars[order(smp.vars$superkingdom, smp.vars$phylum, smp.vars$class, smp.vars$order, smp.vars$family, smp.vars$genus, smp.vars$species),]
smp.vars1 <- smp.vars
smp.vars1[,2:8] <- sapply(smp.vars[,2:8], function(a)(as.numeric(as.factor(a))))

smp.vars1 <- melt(smp.vars)
textcol <- "grey40"
color.scale <- c('samples.variations' = 'darkolivegreen2', 'regions.variations' = 'orchid')

smp.vars1.1 <- smp.vars1 %>%
  group_by(variable) %>%
  mutate(value.alpha = scale(value)) %>%
  ungroup() %>%
  arrange(variable)
smp.vars1.1$species <- factor(smp.vars1.1$species, levels = unique(smp.vars1.1$species))

p <- ggplot(smp.vars1.1, aes(variable, species, fill=variable, alpha = value.alpha)) + geom_tile(colour="white",size=0.2)+
  geom_text(aes(label=value), alpha = 1, cex = 2) +
  # guides(fill=guide_legend(title="Cases per\n100,000 people"))+
  # labs(x="",y="",title="Incidence of Measles in the US")+
  scale_y_discrete(limits = rev(levels(smp.vars1.1$species)), position = "right") +
  # scale_x_discrete(limits = rev(levels(species)))+
  scale_fill_manual(values = color.scale) +
  scale_alpha(range = c(0.1, 0.9)) +  coord_fixed() + labs(x = NULL, y = NULL) +
  theme(legend.position="none",legend.direction="vertical", legend.title=element_text(colour=textcol),legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),legend.key.height=grid::unit(0.8,"cm"),legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol, vjust = 0.25, angle = 90, face = 'bold'), axis.text.y=element_text(vjust=0.2,colour=textcol, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_blank(),plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"))

p
smp.vars.leg <- smp.vars1.1 %>% group_by(variable) %>% 
  summarise(min = min(value), max = max(value)) %>%
  mutate(tile.1 = 0.1, tile.3 = 0.3, tile.5 = 0.5, tile.7 = 0.7, tile.9 = 0.9) %>%
  tidyr::gather(tile, value, -variable, -min, -max)
p.legend <- ggplot(smp.vars.leg, aes(x = value, y = variable, fill = variable, alpha = value)) +
  geom_tile() +
  geom_text(aes(label = variable), x = -0.1) +
  geom_text(aes(label = min), x = 0.1) +
  geom_text(aes(label = max), x = 0.9, color = "white") +
  ggtitle("Legend") +
  coord_cartesian(xlim = c(-0.2, 1.1)) +
  scale_fill_manual(values = color.scale) +
  scale_alpha_identity() +
  theme_void() +
  theme(legend.position = "none")
p.legend

cowplot::plot_grid(p, p.legend, 
                   ncol = 2, align = "v", 
                   rel_heights = c(5, 1))

## merge deletions and variations for plotting ====
smp.SVs1.1$species == smp.dels1.1$species
smp.dels1.1$species == smp.vars1.1$species
full.df <- rbind(smp.SVs1.1, smp.dels1.1, smp.vars1.1)

color.scale <- c('samples.SVs' = 'darkolivegreen3', 'regions.SVs' = 'hotpink2', 'samples.deletions' = 'blue4', 'regions.deletions' = 'darkred', 'samples.variations' = 'darkorange1', 'regions.variations' = 'darkorchid')

p <- ggplot(full.df, aes(variable, species, fill=variable, alpha = value.alpha)) + geom_tile(colour="white",size=0.2)+
  geom_text(aes(label=value), alpha = 1, cex = 2) +
  # guides(fill=guide_legend(title="Cases per\n100,000 people"))+
  # labs(x="",y="",title="Incidence of Measles in the US")+
  scale_y_discrete(limits = rev(levels(full.df$species)), position = "right") +
  # scale_x_discrete(limits = rev(levels(species)))+
  scale_fill_manual(values = color.scale, guide = 'legend') +
  scale_alpha(range = c(0.1, 0.9), guide = 'legend') +  coord_fixed() + labs(x = NULL, y = NULL, title = 'ALL (n = 147)') +
  theme(legend.position="none",legend.direction="vertical", legend.title=element_text(colour=textcol),legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),legend.key.height=grid::unit(0.8,"cm"),legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol, vjust = 0.25, angle = 90, face = 'bold'), axis.text.y=element_text(vjust=0.2,colour=textcol, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_blank(),plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"))
pdf('ALLSVdistributions.pdf', height = 12, width = 6)
p
dev.off()

full.legend <- full.df %>% group_by(variable) %>% 
  summarise(min = min(value), max = max(value)) %>%
  mutate(tile.1 = 0.1, tile.3 = 0.3, tile.5 = 0.5, tile.7 = 0.7, tile.9 = 0.9) %>%
  tidyr::gather(tile, value, -variable, -min, -max)
p.legend <- ggplot(vegans.legend, aes(x = value, y = variable, fill = variable, alpha = value)) +
  geom_tile() +
  geom_text(aes(label = variable), x = -0.1) +
  geom_text(aes(label = min), x = 0.1) +
  geom_text(aes(label = max), x = 0.9, color = "white") +
  ggtitle("Legend") +
  coord_cartesian(xlim = c(-0.2, 1.1)) +
  scale_fill_manual(values = color.scale) +
  scale_alpha_identity() +
  theme_void() +
  theme(legend.position = "none")
pdf('legend.pdf')
p.legend
dev.off()
