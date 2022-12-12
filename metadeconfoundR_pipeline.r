library(metadeconfoundR)
setwd("/Users/fjw536/Desktop/Anorexia/AN_DA") # Change directory location when needed
MSP <- data.frame(t(read.delim("MSP_matrix.txt", header = T, row.names = 1)))
MSP.anno <- read.delim("IGC2.1989MSPs.taxo.tsv", header = T, row.names = 1)
meta <- read.delim("Pheno.txt", header = T, row.names = 1)
drug_use <- read.delim('Drug_use.txt', header = T, row.names = 1)

# Remove MSPs with prevalence lower than 10%
msp.filter <- filter_pre(MSP.final, 0.1, 'row')

# Combine drug intake data and metadata
drug <- colnames(drug_use)
meta[, drug] <- 0

for (i in 1:nrow(drug_use)) {
  meta[rownames(drug_use)[i], drug] <- drug_use[i, ]
} 

meta <- cbind(meta, drug_use)

# Select drug intake as sub metadata
meta.sub <- meta[, c(10, 2, 28, 29, 33:64)]
meta.sub$Category <- as.factor(meta.sub$Category)
levels(meta.sub$Category) <- c('1', '0')


# MSPs drug deconfounding
# Check if the rownames are matched between msp and meta.sub matrix
all(rownames(msp.filter)==rownames(meta.sub))
all(order(rownames(msp.filter)) == order(rownames(meta.sub)))

AN_drugdeconfound_output_msp <- MetaDeconfound(featureMat = msp.filter,
                                           metaMat = meta.sub)

BuildHeatmap(AN_drugdeconfound_output_msp, cuneiform = T, d_range = 'full')

AN_deconfound <- as.data.frame(AN_drugdeconfound_output_msp$status)
table(AN_deconfound$Category)
AN_deconfound <- AN_deconfound %>% filter(Category %in% c('AD', 'OK_nc', 'OK_sd'))

# Add annotation
AN_deconfound_sig$'MSP_annotation' <- MSP.anno[rownames(AN_deconfound_sig), 'species']
