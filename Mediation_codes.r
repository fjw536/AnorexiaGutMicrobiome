# Codes to perform mediation analysis
micro_ED <- lm_btw_mats(msp.an, meta[an, EDI], meta[an, ], covar = c(), direction = c(1,1), y_mat = 1)
micro_ED.sig <- micro_ED$table %>% filter(fdr.p < 0.1)

micro_metabolome <- lm_btw_mats(msp.an[, unique(micro_ED.sig$Phenotype)], metabolome.an, meta[an, ], covar = c(), direction = c(1,1), y_mat = 1)
micro_metabolome.sig <- micro_metabolome$table %>% filter(fdr.p < 0.1)

metabolome_ED <- lm_btw_mats(metabolome.an[, unique(micro_metabolome.sig$Taxa)], meta[an, EDI],  meta[an, ], covar = c(), direction = c(1,1), y_mat = 1)
metabolome_ED.sig <- metabolome_ED$table %>% filter(fdr.p < 0.1)

set.seed(100)
micro.ol <- intersect(unique(micro_ED.sig$Phenotype), unique(micro_metabolome.sig$Phenotype))
ED.ol <- intersect(unique(micro_ED.sig$Taxa), unique(metabolome_ED.sig$Taxa))
metabolome.ol <- intersect(unique(micro_metabolome.sig$Taxa), unique(metabolome_ED.sig$Phenotype))

micro.med <- c()
ED.med <-c()
metabolome.med <- c()
for (i in micro.ol) {
  for (j in ED.ol) {
    for (k in metabolome.ol) {
      if (j %in% micro_ED.sig[which(micro_ED.sig$Phenotype == i), 'Taxa'] & 
          k %in% micro_metabolome.sig[which(micro_metabolome.sig$Phenotype == i), 'Taxa'] &
          j %in% metabolome_ED.sig[which(metabolome_ED.sig$Phenotype == k), 'Taxa']) {
        micro.med <- c(micro.med, i)
        ED.med <- c(ED.med, j)
        metabolome.med <- c(metabolome.med, k)
      }
    }
  }
}


# Generate the matrix of candidate mediation groups
med <- matrix(NA, nrow = length(micro.med), ncol = 3) %>% data.frame
colnames(med) <-c('Microbial', 'Metabolites', 'ED')
med$Microbial <- micro.med
med$Metabolites <- metabolome.med
med$ED <- ED.med

med$'Annotation' <- MSP.anno[med$Microbial, 'species']

# Data transformation
micro.qtrans <- apply(msp.an, 2, qtrans) %>% data.frame
metabolome.qtrans <- apply(metabolome.an, 2, qtrans) %>% data.frame
ED.qtrans <- apply(meta[an, EDI], 2, qtrans) %>% data.frame
cov.qtrans <- apply(meta[an, cov], 2, qtrans) %>% data.frame

micro.qtrans <- apply(micro.qtrans, 2, as.numeric) %>% data.frame
metabolome.qtrans <- apply(metabolome.qtrans, 2, as.numeric) %>% data.frame
ED.qtrans <- apply(ED.qtrans, 2, as.numeric) %>% data.frame
cov.qtrans <- apply(cov.qtrans, 2, as.numeric) %>% data.frame

set.seed(100)
# This is an example code testing the hypothesis where phenotype impact microbiome through metabolites
# Now start the mediation analysis (direction 1 : phenotype --> metabolites --> msp)
library(mediation)
pre <- read.delim("/Users/fjw536/Desktop/Anorexia/AN_DA/pre.txt", row.names = 1) # Import a matrix of sample IDs
res <- data.frame(Ass=0, ACME_beta=0,ACME_p=0, ADE_beta=0, ADE_p=0, Total_beta=0, Total_p=0, Prop=0, Prop_p=0)
for (i in 1:nrow(med)) {
  dat = cbind(pre[1:77,], ED.qtrans[,med[i,3]], metabolome.qtrans[,med[i,2]], micro.qtrans[,med[i,1]], cov.qtrans) %>% na.omit %>% data.frame
  dat[,2:10] <- apply(dat[,2:10], 2, as.numeric) %>% data.frame
  colnames(dat)[1] = 'ID'
  colnames(dat)[2] = med[i,1]
  colnames(dat)[3] = med[i,2]
  colnames(dat)[4] = med[i,3]
  colnames(dat)[5:9] = cov
  b = lm(dat[,3]~dat[,2]+ dat[, 5] + dat[, 6] + dat[, 7] + dat[, 8] + dat[, 9]+ dat[, 10], dat)  # Compute a linear model of for the associations (cause to the mediator)
  c = lm(dat[,4]~dat[,3] + dat[,2] + dat[, 5] + dat[, 6] + dat[, 7] + dat[, 8] + dat[, 9]+ dat[, 10], dat)  # Another linear model for the association (pheno ~ mediator + cause)
  med_out = summary(mediate(b, c, treat = "dat[, 2]", mediator = "dat[, 3]", sims = 10))
  tmp=c(paste(med[i, 1], med[i, 2], med[i,3], sep = "_"), med_out$d0, med_out$d0.p, med_out$z0, med_out$z0.p, med_out$tau.coef, med_out$tau.p, med_out$n0, med_out$n0.p)
  res <- rbind(res, tmp)
  remove(tmp)
}

res <- res[-1, ]
res$"ACME_p.adj" <- p.adjust(res$ACME_p)
res$"ADE_p.adj" <- p.adjust(res$ADE_p)
res$"Total_p.adj" <- p.adjust(res$Total_p)
res$"Prop_p.adj" <- p.adjust(res$Prop_p)

### Mediation analysis linking microbiome to phenotyes through metabolome (direction 1 : phenotype --> msp --> metabolites)
res.2 <- data.frame(Ass=0, ACME_beta=0,ACME_p=0, ADE_beta=0, ADE_p=0, Total_beta=0, Total_p=0, Prop=0, Prop_p=0)
for (i in 1:nrow(med)) {
  dat = cbind(pre[1:77,], ED.qtrans[,med[i,3]], metabolome.qtrans[,med[i,2]], micro.qtrans[,med[i,1]], cov.qtrans) %>% na.omit %>% data.frame
  dat[,2:10] <- apply(dat[,2:10], 2, as.numeric) %>% data.frame
  colnames(dat)[1] = 'ID'
  colnames(dat)[2] = med[i,1]
  colnames(dat)[3] = med[i,2]
  colnames(dat)[4] = med[i,3]
  colnames(dat)[5:9] = cov
  b = lm(dat[,4]~dat[,2] + dat[, 5] + dat[, 6] + dat[, 7] + dat[, 8] + dat[, 9]+ dat[, 10], dat)  # Compute a linear model of for the associations (cause to the mediator)
  c = lm(dat[,3]~dat[,4] + dat[,2] + dat[, 5] + dat[, 6] + dat[, 7] + dat[, 8] + dat[, 9]+ dat[, 10], dat)  # Another linear model for the association (pheno ~ mediator + cause)
  med_out = summary(mediate(b, c, treat = "dat[, 2]", mediator = "dat[, 4]", sims = 10))
  tmp=c(paste(med[i, 1], med[i, 2], med[i,3], sep = "_"), med_out$d0, med_out$d0.p, med_out$z0, med_out$z0.p, med_out$tau.coef, med_out$tau.p, med_out$n0, med_out$n0.p)
  res.2 <- rbind(res.2, tmp)
  remove(tmp)
}

res.2 <- res.2[-1, ]
res.2$"ACME_p.adj" <- p.adjust(res.2$ACME_p)
res.2$"ADE_p.adj" <- p.adjust(res.2$ADE_p)
res.2$"Total_p.adj" <- p.adjust(res.2$Total_p)
res.2$"Prop_p.adj" <- p.adjust(res.2$Prop_p)

# Number of inferred causal relationships for direction 1 (from SGVs to phenotypes), direction 2 (from SGVs to Metabolites), and both.

med.res <- cbind(res, res.2[, -1])
colnames(med.res) <- c("Ass", "d1.ACME_beta", 
                       "d1.ACME_p", "d1.ADE_beta", 
                       "d1.ADE_p", "d1.Total_beta", 
                       "d1.Total_p", "d1.Prop", 
                       "d1.Prop_p", "d1.ACME_p.adj", 
                       "d1.ADE_p.adj", "d1.Total_p.adj", 
                       "d1.Prop_p.adj", "d2.ACME_beta", 
                       "d2.ACME_p", "d2.ADE_beta", 
                       "d2.ADE_p", "d2.Total_beta", 
                       "d2.Total_p", "d2.Prop", 
                       "d2.Prop_p", "d2.ACME_p.adj", 
                       "d2.ADE_p.adj", "d2.Total_p.adj", "d2.Prop_p.adj")
med.res$d1.ACME_p.adj <- p.adjust(med.res$d1.ACME_p, method = "BH")
med.res$d2.ACME_p.adj <- p.adjust(med.res$d2.ACME_p, method = "BH")
med.res$d1.Prop_p.adj <- p.adjust(med.res$d1.Prop_p, method = "BH")
med.res$d2.Prop_p.adj <- p.adjust(med.res$d2.Prop_p, method = "BH")
med.res$"med.res.dir"[med.res$d1.ACME_p.adj < 0.05 &
                        med.res$d1.Prop_p.adj < 0.05 &
                        med.res$d2.ACME_p.adj < 0.05 &
                        med.res$d2.Prop_p.adj < 0.05] <- "both"
med.res$"med.res.dir"[med.res$d1.ACME_p.adj < 0.05 &
                        med.res$d1.Prop_p.adj < 0.05 &
                        med.res$d2.ACME_p.adj > 0.05] <- "d1"
med.res$"med.res.dir"[med.res$d1.ACME_p.adj > 0.05 &
                        med.res$d2.ACME_p.adj < 0.05 &
                        med.res$d2.Prop_p.adj < 0.05] <- "d2"
table(med.res$med.res.dir)
med.res.plot.d1 <- filter(med.res, med.res$med.res.dir=="d1")
med.res.plot.d2 <- filter(med.res, med.res$med.res.dir=="d2")
