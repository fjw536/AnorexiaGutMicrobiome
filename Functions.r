## The Normal Quantile Transformation
library(dplyr)
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}

# Filter function for prevalence
filter_pre <- function(dat, prevalence, direction){
  if (direction=="column") {dat = t(dat)}
  dat[is.na(dat)] <- 0
  dat_logic <- dat==0
  colsig <- colSums(dat_logic)
  dat_sig <- names(colsig[colsig < nrow(dat)*(1- as.numeric(prevalence))])
  dat_filter <- as.data.frame(dat[,dat_sig])
  return(dat_filter)
}

## linear model
lm_btw_mats<-function(mat0,mat1,mat2,covar, direction = c(1,1),y_mat = 0){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  # direction: 1 means samples in row and variables in column; 2 means samples in column and variables in row
  # y_mat: 0 means mat0 is y in linear model; 1 mean mat1 is y in linear model
  
  require(reshape2)
  
  ## test block
  #mat0<-lld_tmao[,c(1:3)]
  #mat1<-lld_vsv[,c(1:3)]
  #mat2<-lld_basic
  #covar<-covar1
  #direction<-c(1,1)
  #y_mat <- 0
  ## test block
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  
  col_var<-mat0
  row_var<-mat1
  
  my_lm<-function(a,b){
    # a<-col_var[,1]
    # b<-row_var[,1]
    
    beta    <- NA
    p.value <- NA
    N       <- NA
    uniq_N  <- NA
    
    lm_input<-data.frame(phen = a,
                         mbio = b,
                         mat2[,covar])
    lm_input <- sapply(lm_input, as.numeric)
    lm_input <- na.omit(lm_input)
    N        <- nrow(lm_input)
    uniq_N   <- length(unique(lm_input[,1]))
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame()
    
    
    if(y_mat==0){
      try(lm_res <- summary(lm(phen~.,data = lm_input)), silent = T)
      indv<-'mbio'
    }else{
      try(lm_res <- summary(lm(mbio~.,data = lm_input)), silent = T)
      indv<-'phen'
    }
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    try(return(list(beta = beta, se = se, p.value = p.value, N = N, uniq_N = uniq_N)), silent = T)
  }
  
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )
  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol = 5, byrow = T)
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # se matrix
  col_var_row_var.se<-matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.se)<-colnames(col_var)
  rownames(col_var_row_var.se)<-colnames(row_var)
  #col_var_row_var.se[is.na(col_var_row_var.se)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,3],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  
  
  # N matrix
  col_var_row_var.N <- matrix(col_var_row_var.unlist[,4],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.N)<-colnames(col_var)
  rownames(col_var_row_var.N)<-colnames(row_var)
  col_var_row_var.N[is.na(col_var_row_var.N)]<-0
  
  # uniq_N matrix
  col_var_row_var.uN <- matrix(col_var_row_var.unlist[,5],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.uN)<-colnames(col_var)
  rownames(col_var_row_var.uN)<-colnames(row_var)
  col_var_row_var.uN[is.na(col_var_row_var.uN)]<-0
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_se     <- melt(col_var_row_var.se)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  col_var_row_var_N          <- melt(col_var_row_var.N)
  col_var_row_var_uN         <- melt(col_var_row_var.uN)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,
                              col_var_row_var_edge_se,
                              col_var_row_var_edge_p,
                              col_var_row_var_N,
                              col_var_row_var_uN)[,-c(4,5,7,8,10,11,13,14)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "Beta","SE", "p","N","uniq_N")
  
  
  # add p adjust
  col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
                                   fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
                                   bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$fdr.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  # bonferroni matrix
  col_var_row_var.bon<-matrix(data  = col_var_row_var_edge$bonferroni.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  
  colnames(col_var_row_var.bon)<-colnames(col_var)
  rownames(col_var_row_var.bon)<-colnames(row_var)
  col_var_row_var.bon[is.na(col_var_row_var.bon)]<-1
  
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  return(list(table      = col_var_row_var_edge,
              beta       = col_var_row_var.beta,
              se         = col_var_row_var.se,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr,
              bonferroni = col_var_row_var.bon,
              N          = col_var_row_var.N,
              uniq_N     = col_var_row_var.uN))
}


# Run PCA
library(ade4)
library(ggplot2)
library(RColorBrewer)
run_pca <- function(metag, meta, Group){
  pca<- dudi.pca(metag, scal = FALSE, scan = FALSE)
  pca_eig <- (pca$eig)[1:2] / sum(pca$eig)
  sample_site <- data.frame({pca$li})[1:2]
  sample_site$names <- rownames(sample_site)
  names(sample_site)[1:2] <- c('PCA1', 'PCA2')
  
  sample_site$level<-factor(meta[, Group])
  pca_plot <- ggplot(sample_site, aes(PCA1, PCA2,color=level)) +
    theme_classic()+ #Remove background box
    geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
    geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
    geom_point(size = 3)+  # Change transparency of dots
    scale_color_manual(values = brewer.pal(6,"Set2")) + # Color dots
    theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.title=element_blank()
    )+
    labs(x = paste('PCA1: ', round(100 * pca_eig[1], 2), '%'), y = paste('PCA2: ', round(100 * pca_eig[2], 2), '%'))
  
  return(pca_plot)
}

# Correlation between two matrix
cor_mat <- function(mat1, mat2, cor, direction = c(1,1)){
  if(direction[1]==2){mat1 =t(mat1)}
  if(direction[2]==2){mat2 =t(mat2)}
  tmpMat=array(NA, c(length(colnames(mat1)), length(colnames(mat2)),2))
  dimnames(tmpMat) [[1]]=colnames(mat1)
  dimnames(tmpMat) [[2]]=colnames(mat2)
  dimnames(tmpMat) [[3]]=c("estimate", "p.value")
  for (i in 1:length(colnames(mat1))) {
    for (j in 1:length(colnames(mat2))) {
      data <- as.data.frame(cbind(mat1[, i], mat2[, j])) %>% na.omit
      # if (nrow(data) < 3) {
      #   tmpMat[colnames(mat1)[i],colnames(mat2)[j],"estimate"]= NA
      #   tmpMat[colnames(mat1)[i],colnames(mat2)[j],"p.value"]= NA
      #   remove(data)
      # } else {
      tmpMat[colnames(mat1)[i],colnames(mat2)[j],"estimate"]=cor.test(as.numeric(data[, 1]), as.numeric(data[, 2]), method = cor, exact = F)$estimate
      tmpMat[colnames(mat1)[i],colnames(mat2)[j],"p.value"]=cor.test(as.numeric(data[, 1]), as.numeric(data[, 2]), method = cor, exact = F)$p.value
      remove(data)
    }
  }
  cor_traits=tmpMat
  cor_traits_p=cor_traits[,,"p.value"]
  cor_traits_estimate=cor_traits[,,"estimate"]
  cor_traits_p_adj=p.adjust(cor_traits_p, method = "bonferroni")
  dim(cor_traits_p_adj)=dim(cor_traits_p)
  rownames(cor_traits_p_adj)=rownames(cor_traits_p)
  colnames(cor_traits_p_adj)=colnames(cor_traits_p)
  tmp=cor_traits_p<0.01
  issig=rowSums(tmp, na.rm = T)
  rsig=na.omit(names(issig[issig>0]))
  issig=colSums(tmp, na.rm = T)
  csig=na.omit(names(issig[issig>0]))
  try(return(list(rho = cor_traits_estimate[rsig, csig], p.value = cor_traits_p, BH = cor_traits_p[rsig, csig])), silent = T)
}
