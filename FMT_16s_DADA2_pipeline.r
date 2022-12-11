library(dada2)
setwd("/Users/fjw536/Downloads/AN_FMT_16s/")
list.files()
fnFs <- sort(list.files(pattern = 'L001_R1_001.fastq', full.names = T))
fnFs
fnRs <- sort(list.files(pattern = 'L001_R2_001.fastq', full.names = T))
fnRs
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnFs[c(5:6, 9:10)])
plotQualityProfile(fnFs[25:26])

plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnRs[c(5:6, 9:10)])
plotQualityProfile(fnRs[c(25:26, 35:36)])

filtFs <- file.path("Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs
dadaRs

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
head(taxa)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.names <- paste('ASV_', seq(1:1828), sep = '')
rownames(taxa.print) <- taxa.names

head(seqtab.nochim)
dim(seqtab.nochim)

samples.out <- rownames(seqtab.nochim)


seqtab.nochim[1,]

library(dplyr)
metag <- seqtab.nochim
colnames(metag) <- taxa.names
head(metag)


# Output the combined matrix
metag <- cbind(metag, taxa.print) %>% data.frame

write.table(metag, file = '16s_matrix.txt', quote = F, sep = '\t')

# Reload data matrix
library(dplyr)
metag <- read.delim('16S_matrix_final.txt', header = T, row.names = 1) %>% t %>% data.frame
taxa <- read.delim('Taxa_info.txt', header = T, row.names = 1)

metag.new <- apply(metag, 1, function(x){x/sum(x)}) %>% t %>% data.frame

# Filter ASVs with prevalence of 0
filter_pre <- function(dat, prevalence, direction){
  if (direction=="column") {dat = t(dat)}
  dat[is.na(dat)] <- 0
  dat_logic <- dat==0
  colsig <- colSums(dat_logic)
  dat_sig <- names(colsig[colsig < nrow(dat)*(1- as.numeric(prevalence))])
  dat_filter <- as.data.frame(dat[,dat_sig])
  return(dat_filter)
}

# Select dataset without inoculations and run data fitlering
metag.filter <- filter_pre(metag.new[-c(1:4, 33:36), ], 0.05, 'row') # Yield 358 ASVs

metag.logic <- metag.filter > 0

logic_res <- matrix(0, nrow = ncol(metag.filter), ncol = 3) %>% data.frame
rownames(logic_res) <- colnames(metag.filter)
colnames(logic_res) <- c('Human_only', 'Human_mouse', 'Mouse_only')
for (i in 1:ncol(metag.logic)) {
  x <- sum(metag.logic[1:20, i])
  y <- sum(metag.logic[21:28, i])
  if (x == 0 & y > 0) {
    logic_res[i, 1] <- 1
  } else if (x > 0 & y >0){
    logic_res[i, 1] <- 1
    logic_res[i, 2] <- 1
    logic_res[i, 3] <- 1
  } else {
    logic_res[i, 3] <- 1
  }
}

colSums(logic_res)
library(dplyr)
set.human <- logic_res %>% filter(Human_only == 1) %>% rownames
set.mouse <- logic_res %>% filter(Mouse_only == 1) %>% rownames
library(RColorBrewer); library(VennDiagram)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(set.human, set.mouse),
  category.names = c("Human donors" , "Germ-free mice"),
  filename = 'test.tiff',
  output=TRUE,
  imagetype="tiff" ,
  height = 1000, 
  width = 1000, 
  resolution = 600,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c("#B3E2CD", "#FDCDAC"),
  cex = 0.4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)


# Split the dataset into donors and recipients
# Batch 1 AN donors
b1.an <- metag[c(5:8, 1), ]
# Remove columns that have 0 across all entries for both an and nc
b1.an.logic <- b1.an > 0
b1.an <- b1.an[, -which(colSums(b1.an.logic) == 0)] # Remove the ASVs with all entries of 0
b1.an.res <- matrix(0, ncol = 3, nrow = ncol(b1.an)) %>% data.frame()
rownames(b1.an.res) <- colnames(b1.an)
colnames(b1.an.res) <- c('Taxa_in_donor', 'Taxa_in_donor_and_recipient', 'Taxa_in_recipient')
for (i in 1:ncol(b1.an)) {
  if (b1.an[5, i] > 0) {
    b1.an.res[i, 1] <- 1
    if (b1.an[5, i] > 0 & mean(as.numeric(b1.an[1:4, i]))> 0) {
      b1.an.res[i, 2] <- 1
    }
  }
  if (mean(as.numeric(b1.an[1:4, i])) > 0) {
    b1.an.res[i, 3] <- 1
  }
}

table(b1.an.res$Taxa_in_donor) # 80
table(b1.an.res$Taxa_in_donor_and_recipient) # 18
table(b1.an.res$Taxa_in_recipient) # 286

# Batch 1 NC donors
b1.nc <- metag[c(9:12, 33), ]
b1.nc.logic <- b1.nc > 0
b1.nc <- b1.nc[, -which(colSums(b1.nc.logic) == 0)] # Remove the ASVs with all entries of 0
b1.nc.res <- matrix(0, ncol = 3, nrow = ncol(b1.nc)) %>% data.frame()
rownames(b1.nc.res) <- colnames(b1.nc)
colnames(b1.nc.res) <- c('Taxa_in_donor', 'Taxa_in_donor_and_recipient', 'Taxa_in_recipient')
for (i in 1:ncol(b1.nc)) {
  if (b1.nc[5, i] > 0) {
    b1.nc.res[i, 1] <- 1
    if (b1.nc[5, i] > 0 & mean(as.numeric(b1.nc[1:4, i]))> 0) {
      b1.nc.res[i, 2] <- 1
    }
  }
  if (mean(as.numeric(b1.nc[1:4, i])) > 0) {
    b1.nc.res[i, 3] <- 1
  }
}

table(b1.nc.res$Taxa_in_donor) # 92
table(b1.nc.res$Taxa_in_donor_and_recipient) # 24
table(b1.nc.res$Taxa_in_recipient) # 315

# Batch 3 AN donors
b3.an <- metag[c(13:15, 3), ]
# Remove columns that have 0 across all entries for both an and nc
b3.an.logic <- b3.an > 0
b3.an <- b3.an[, -which(colSums(b3.an.logic) == 0)] # Remove the ASVs with all entries of 0
b3.an.res <- matrix(0, ncol = 3, nrow = ncol(b3.an)) %>% data.frame()
rownames(b3.an.res) <- colnames(b3.an)
colnames(b3.an.res) <- c('Taxa_in_donor', 'Taxa_in_donor_and_recipient', 'Taxa_in_recipient')
for (i in 1:ncol(b3.an)) {
  if (b3.an[4, i] > 0) {
    b3.an.res[i, 1] <- 1
    if (b3.an[4, i] > 0 & mean(as.numeric(b3.an[1:3, i]))> 0) {
      b3.an.res[i, 2] <- 1
    }
  }
  if (mean(as.numeric(b3.an[1:3, i])) > 0) {
    b3.an.res[i, 3] <- 1
  }
}

table(b3.an.res$Taxa_in_donor) # 64
table(b3.an.res$Taxa_in_donor_and_recipient) # 15
table(b3.an.res$Taxa_in_recipient) # 187



# Batch 3 NC donors
b3.nc <- metag[c(16:18, 35), ]
b3.nc.logic <- b3.nc > 0
b3.nc <- b3.nc[, -which(colSums(b3.nc.logic) == 0)] # Remove the ASVs with all entries of 0
b3.nc.res <- matrix(0, ncol = 3, nrow = ncol(b3.nc)) %>% data.frame()
rownames(b3.nc.res) <- colnames(b3.nc)
colnames(b3.nc.res) <- c('Taxa_in_donor', 'Taxa_in_donor_and_recipient', 'Taxa_in_recipient')
for (i in 1:ncol(b3.nc)) {
  if (b3.nc[4, i] > 0) {
    b3.nc.res[i, 1] <- 1
    if (b3.nc[4, i] > 0 & mean(as.numeric(b3.nc[1:3, i]))> 0) {
      b3.nc.res[i, 2] <- 1
    }
  }
  if (mean(as.numeric(b3.nc[1:3, i])) > 0) {
    b3.nc.res[i, 3] <- 1
  }
}

table(b3.nc.res$Taxa_in_donor) # 92
table(b3.nc.res$Taxa_in_donor_and_recipient) # 12
table(b3.nc.res$Taxa_in_recipient) # 192

# Batch 4 AN donors
b4.an <- metag[c(19:21, 4), ]
# Remove columns that have 0 across all entries for both an and nc
b4.an.logic <- b4.an == 0
b4.an <- b4.an[, -which(colSums(b4.an.logic) == 0)] # Remove the ASVs with all entries of 0
b4.an.res <- matrix(0, ncol = 3, nrow = ncol(b4.an)) %>% data.frame()
rownames(b4.an.res) <- colnames(b4.an)
colnames(b4.an.res) <- c('Taxa_in_donor', 'Taxa_in_donor_and_recipient', 'Taxa_in_recipient')
for (i in 1:ncol(b4.an)) {
  if (b4.an[4, i] > 0) {
    b4.an.res[i, 1] <- 1
    if (b4.an[4, i] > 0 & mean(as.numeric(b4.an[1:3, i]))> 0) {
      b4.an.res[i, 2] <- 1
    }
  }
  if (mean(as.numeric(b4.an[1:3, i])) > 0) {
    b4.an.res[i, 3] <- 1
  }
}

table(b4.an.res$Taxa_in_donor) # 109
table(b4.an.res$Taxa_in_donor_and_recipient) # 25
table(b4.an.res$Taxa_in_recipient) # 245

# Batch 4 NC donors
b4.nc <- metag[c(22:24, 36), ]
b4.nc.logic <- b4.nc > 0
b4.nc <- b4.nc[, -which(colSums(b4.nc.logic) == 0)] # Remove the ASVs with all entries of 0
b4.nc.res <- matrix(0, ncol = 3, nrow = ncol(b4.nc)) %>% data.frame()
rownames(b4.nc.res) <- colnames(b4.nc)
colnames(b4.nc.res) <- c('Taxa_in_donor', 'Taxa_in_donor_and_recipient', 'Taxa_in_recipient')
for (i in 1:ncol(b4.nc)) {
  if (b4.nc[4, i] > 0) {
    b4.nc.res[i, 1] <- 1
    if (b4.nc[4, i] > 0 & mean(as.numeric(b4.nc[1:3, i]))> 0) {
      b4.nc.res[i, 2] <- 1
    }
  }
  if (mean(as.numeric(b4.nc[1:3, i])) > 0) {
    b4.nc.res[i, 3] <- 1
  }
}

table(b4.nc.res$Taxa_in_donor) # 78
table(b4.nc.res$Taxa_in_donor_and_recipient) # 23
table(b4.nc.res$Taxa_in_recipient) # 215


# Visualization using Venn Diagram
library(dplyr); library(RColorBrewer); library(VennDiagram)
my_venn <- function(input.res, ID){
  set.human <- input.res %>% filter(Taxa_in_donor == 1) %>% rownames
  set.mouse <- input.res %>% filter(Taxa_in_recipient == 1) %>% rownames
  myCol <- brewer.pal(3, "Pastel2")
  venn.diagram(
    x = list(set.mouse, set.human),
    category.names = c( "Germ-free \n recipients", 
                        "Human \n donor"),
    filename = paste0(ID, '.tiff'),
    output=TRUE,
    imagetype="tiff" ,
    height = 800, 
    width = 800, 
    resolution = 800,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = c("#B3E2CD", "#FDCDAC"),
    cex = 0.2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  )
}



# Combined batch as a whole for the comparison
com.b <- metag[c(1, 3:24, 33, 35:36),]
com.b.logic <- com.b > 0
com.b <- com.b[, -which(colSums(com.b.logic) == 0)] # Remove the ASVs with all entries of 0
com.b.rel <- apply(com.b, 1, my_rel) %>% as.data.frame()

com.b.res <- matrix(0, ncol = 3, nrow = ncol(com.b)) %>% data.frame()
rownames(com.b.res) <- colnames(com.b)
colnames(com.b.res) <- c('Taxa_in_donor', 'Taxa_in_donor_and_recipient', 'Taxa_in_recipient')
for (i in 1:ncol(com.b)) {
  if (mean(as.numeric(com.b[1:3, i])) > 0 | mean(as.numeric(com.b[24:26, i])) > 0) {
    com.b.res[i, 1] <- 1
    if (mean(as.numeric(com.b[c(4:7, 12:14, 18:20), i])) > 0 & 
        mean(as.numeric(com.b[c(8:11, 15:17, 21:23), i])) > 0) {
      com.b.res[i, 2] <- 1
    }
  }
  if (mean(as.numeric(com.b[c(4:7, 12:14, 18:20), i])) > 0 & 
      mean(as.numeric(com.b[c(8:11, 15:17, 21:23), i])) > 0) {
    com.b.res[i, 3] <- 1
  }
}

table(com.b.res$Taxa_in_donor) # 161
table(com.b.res$Taxa_in_donor_and_recipient) # 36
table(com.b.res$Taxa_in_recipient) # 92

com.b.asv <- com.b.res %>% filter(Taxa_in_donor_and_recipient == 1) %>% rownames()
com.b.asv.res <- matrix(NA, nrow = length(com.b.asv), ncol = 5) %>% data.frame
rownames(com.b.asv.res) <- com.b.asv
colnames(com.b.asv.res) <- c('CliffDelta_D', 'p_D', 'p_R', 'CliffDelta_R', 'Color')
for (i in com.b.asv) {
  com.b.asv.res[i, 1] <- cliff.delta(as.numeric(com.b.rel[i,1:3]), as.numeric(com.b.rel[i,24:26]))
  com.b.asv.res[i, 2] <- wilcox.test(as.numeric(com.b.rel[i,1:3]), as.numeric(com.b.rel[i,24:26]))$p.value
  com.b.asv.res[i, 3] <- wilcox.test(as.numeric(com.b.rel[i, c(4:7, 12:14, 18:20)]), as.numeric(com.b.rel[i, c(8:11, 15:17, 21:23)]))$p.value
  com.b.asv.res[i, 4] <- cliff.delta(as.numeric(com.b.rel[i, c(4:7, 12:14, 18:20)]), as.numeric(com.b.rel[i, c(8:11, 15:17, 21:23)]))
  if (com.b.asv.res[i,1] > 0 & com.b.asv.res[i, 4] > 0) {
    com.b.asv.res[i, 5] <- 'Yes'
  } else if (com.b.asv.res[i, 1] < 0 & com.b.asv.res[i, 4] < 0) {
    com.b.asv.res[i, 5] <- 'Yes'
  } else {
    com.b.asv.res[i, 5] <- 'No'
  }
}
table(com.b.asv.res$Color)

com.b_D <- ggplot(com.b.asv.res, aes(x=rownames(com.b.asv.res), y=CliffDelta_D, fill = as.factor(Color))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('grey', '#3288ff'))+
  coord_flip() + theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 0),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  geom_hline(yintercept = 0) +
  ylim(-0.7, 0.4) +
  scale_x_discrete(limits = rev(rownames(com.b.asv.res)))


tmp <- cbind(com.b.asv.res, tax_meta[rownames(com.b.asv.res),])
tmp$'Annotation' <- paste(rownames(tmp), ': p_',tmp$Phylum, '|f_', tmp$Family, '|g_', tmp$Genus, sep = '')
com.b.asv.res[, 6] <- tmp$'Annotation'


com.b.asv.res.plot <- com.b.asv.res[order(-com.b.asv.res$CliffDelta_R), ]
com.b.asv.res.plot$V6 <- factor(com.b.asv.res.plot$V6, levels = com.b.asv.res.plot$V6)

com.b_R <- ggplot(com.b.asv.res.plot, aes(x= V6, y=CliffDelta_R, fill = as.factor(Color))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('grey', '#3288ff'))+
  scale_x_discrete(position = "top")+
  coord_flip() + theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 0),
        axis.text.y = element_text(size = 9, color = 'black'),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank())  +
  geom_hline(yintercept = 0)
com.b_R
ggsave(com.b_R, filename = '/Users/fjw536/Desktop/Anorexia/AN_DA/FMT_recipient.bar.pdf', width = 6, height = 14)


asv.order <- rev(rownames(com.b.asv.res.plot))


library(ggpubr)

p_asv_com.b <- ggarrange(com.b_D, com.b_R, ncol = 2, widths = c(0.9, 0.75))
p_asv_com.b
pdf('/Users/fjw536/Desktop/Anorexia/AN_DA/FMT_recipient.bar.pdf', height = 10, width = 2)
print(com.b_R)
dev.off()
write.table(taxa[com.b.asv, ], file = 'contrast_ASV.anno.txt', row.names = T, quote = F, sep = '\t')
p.test <- draw(p_heat_FMT, heatmap_legend_side = "bottom")
getwd()

# Prepare for heatmap
library(pheatmap)
library(paletteer)
donor <- com.b.rel[rownames(com.b.asv.res), c(1:3, 24:26)]
donor.asv <- cbind(donor[, 4:6], donor[, 1:3]) %>% as.matrix()

p_heat_FMT <- pheatmap(donor.asv[asv.order, ], 
                       scale = 'none', 
                       cluster_cols = F, 
                       cluster_rows = F, 
                       color = paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps"), 
                       border_color = 'black', 
                       show_colnames = F,
                       show_rownames = F)
pushViewport(viewport(width = 0.9, height = 0.9))
p_heat_FMT <- Heatmap(donor.asv[asv.order, ], 
                      col = paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps"),
                      heatmap_legend_param = list(title = "Relative abundance", 
                                                  direction = 'horizontal',
                                                  title_position = "topcenter"),
                      na_col = 'grey',
                      cluster_rows = F,
                      cluster_columns = F,
                      show_row_names = F,
                      show_column_names = F,
                      column_gap = unit(3, "mm"),
                      column_split =  c(rep('a', 3), rep('b', 3)),
                      column_title = NULL,
                      border = 'black',
                      rect_gp = gpar(col = "black", lwd = 1))
getwd()
pdf('/Users/fjw536/Desktop/Anorexia/AN_DA/FMT_donor.hp.pdf', width = 2, height = 14)
draw(p_heat_FMT, heatmap_legend_side = "bottom")
dev.off()
