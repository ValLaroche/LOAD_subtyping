# Generate EWAS for methylation data
limma_DA_methyl <- function(beta_values, pheno, subtype1, subtype2, cohort, AD = FALSE){
  cl<- makeCluster(8)
  
  if(AD == FALSE){ 
    tmppheno = pheno[pheno$subtype %in% c(subtype1, subtype2),]
  } else {
    tmppheno = pheno
  }
  
  tmppheno = tmppheno[order(tmppheno$subtype),]
  tmpcounts = beta_values[,colnames(beta_values) %in% tmppheno$Basename]
  tmpcounts = tmpcounts[,match(tmppheno$Basename,colnames(tmpcounts))]
    
  if(AD == FALSE){ 
    group = tmppheno$subtype
  } else {
    group = tmppheno$diag
  }
  tmppheno$group = group
  
  tmppheno = calc_sva_methyl(tmpcounts, tmppheno, cohort)

  group[which(group == subtype1)] = 1
  group[which(group == subtype2)] = 2

  results_lm<-t(parApply(cl, tmpcounts, 1, EWAS, group, tmppheno))
  colnames(results_lm)<-c("Estimate.t1","SE.t1","t.t1","P.t1", "Beta", "SE")

  stopCluster(cl)
  return(results_lm)
  
}

# Generate EWAS for FANS data
limma_DA_FANS <- function(beta_values, tmppheno, subtype1, subtype2){
  cl<- makeCluster(8)
  
  tmppheno = tmppheno[order(tmppheno$subtype),]
  tmpcounts = beta_values[,colnames(beta_values) %in% tmppheno$Basename]
  tmpcounts = tmpcounts[,match(tmppheno$Basename,colnames(tmpcounts))]
  tmppheno$group = tmppheno$subtype
  group = tmppheno$group
  tmppheno = calc_sva_FANS(tmpcounts, tmppheno)
  
  group[which(group == subtype1)] = 1
  group[which(group == subtype2)] = 2
  
  results_lm<-t(parApply(cl, tmpcounts, 1, EWAS, group, tmppheno))
  colnames(results_lm)<-c("Estimate.t1","SE.t1","t.t1","P.t1", "Beta", "SE")
  
  stopCluster(cl)
  return(results_lm)
  
}

# Generate EWAS for SNP data data
limma_DA_SNP <- function(snp_genotype, pheno, subtype1, subtype2, cohort, AD = FALSE){
  cl<- makeCluster(8)
  
  if(AD == FALSE){ 
    tmppheno = pheno[pheno$subtype %in% c(subtype1, subtype2),]
  } else {
    tmppheno = pheno
  }
  
  tmppheno = tmppheno[order(tmppheno$subtype),]
  if(cohort == "ROSMAP"){
    tmpcounts = snp_genotype[,colnames(snp_genotype) %in% tmppheno$sample]  
  } else if(cohort == "UKBBN"){
    tmpcounts = snp_genotype[,colnames(snp_genotype) %in% tmppheno$BBNId]  
  } else if(cohort == "PITT"){
    tmpcounts = snp_genotype[,colnames(snp_genotype) %in% tmppheno$Individual_ID]  
  }
  

  if(AD == FALSE){ 
    group = tmppheno$subtype
  } else {
    group = tmppheno$diag
  }
  tmppheno$group = group
  
  # tmppheno = calc_sva_methyl(tmpcounts, tmppheno, cohort)
  
  group[which(group == subtype1)] = 1
  group[which(group == subtype2)] = 2
  
  results_lm<-t(parApply(cl, tmpcounts, 1, EWAS_SNP, group, tmppheno))
  colnames(results_lm)<-c("Estimate.t1","SE.t1","t.t1","P.t1", "Beta", "SE")
  
  stopCluster(cl)
  return(results_lm)
  
}

# Control for inflation of pvalues
inflation_control <- function(limma_res, subtype2){
  gctrl = gap::gcontrol2(limma_res$P.t1)
  print(gctrl$lambda)
  chisq <- qchisq(1-(limma_res$P.t1),1)
  median(chisq)/qchisq(0.5,1)
  
  if(subtype2 == "Control"){
    bacon_tmp = bacon(effectsizes = as.numeric(limma_res$Estimate.t1), 
                      standarderrors = limma_res$SE.t1, niter = 20000L)
  } else if(subtype2 == "blue"){
    bacon_tmp = bacon(effectsizes = as.numeric(limma_res$Estimate.t1), 
                      standarderrors = limma_res$SE.t1, niter = 10000L)
  }
  
  bacon_fit = bacon::fit(bacon_tmp, n = 100)
  head(pval(bacon_tmp))
  inflation(bacon_tmp)
  bias(bacon_tmp)
  print(plot(bacon_tmp, type="qq"))
  
  limma_res$bacon_effectsize = es(bacon_tmp)[,1]
  limma_res$bacon_stderr = se(bacon_tmp)[,1]
  limma_res$bacon_pval = pval(bacon_tmp)[,1]
  
  gctrl = gap::gcontrol2(limma_res$bacon_pval)
  print(gctrl$lambda)
  chisq <- qchisq(1-(limma_res$bacon_pval),1)
  median(chisq)/qchisq(0.5,1)
  
  return(limma_res)
}

# Perform meta analysis 
meta_analysis <- function(UKBBN_res, PITT_res, ROSMAP_res = NULL,
                          data_type, methyl_type = NA){
  if(data_type == "RNA"){
    rownames(UKBBN_res) = gsub("\\.[0-9][0-9]","",rownames(UKBBN_res))
    rownames(UKBBN_res) = gsub("\\.[0-9]","",rownames(UKBBN_res))
    
    rownames(ROSMAP_res) = gsub("\\.[0-9][0-9]","",rownames(ROSMAP_res))
    rownames(ROSMAP_res) = gsub("\\.[0-9]","",rownames(ROSMAP_res))
    
    rownames(PITT_res) = gsub("\\.[0-9][0-9]","",rownames(PITT_res))
    rownames(PITT_res) = gsub("\\.[0-9]","",rownames(PITT_res))
    
    feature_set = Reduce(intersect, list(rownames(ROSMAP_res),
                                         rownames(UKBBN_res),
                                         rownames(PITT_res)))
    
    UKBBN_res = UKBBN_res[rownames(UKBBN_res) %in% feature_set,]
    ROSMAP_res = ROSMAP_res[rownames(ROSMAP_res) %in% feature_set,]
    PITT_res = PITT_res[rownames(PITT_res) %in% feature_set,]
    
    UKBBN_res = UKBBN_res[order(rownames(UKBBN_res)),]
    ROSMAP_res = ROSMAP_res[order(rownames(ROSMAP_res)),]
    PITT_res = PITT_res[order(rownames(PITT_res)),]
    
    estimate_meta<-cbind(UKBBN_res$bacon_effectsize, ROSMAP_res$bacon_effectsize, PITT_res$bacon_effectsize)
    stderr_meta<-cbind(UKBBN_res$bacon_stderr, ROSMAP_res$bacon_stderr, PITT_res$bacon_stderr)
    data_meta = cbind(estimate_meta, stderr_meta)
    rownames(data_meta) = rownames(UKBBN_res)
    
  } else if(data_type == "Methyl"){
    if(methyl_type == "EPIC") {
      feature_set = Reduce(intersect, list(rownames(BDR_res),
                                           rownames(UKBBN_res),
                                           rownames(PITT_res)))
      
      BDR_res = BDR_res[rownames(BDR_res) %in% feature_set,]
      UKBBN_res = UKBBN_res[rownames(UKBBN_res) %in% feature_set,]
      PITT_res = PITT_res[rownames(PITT_res) %in% feature_set,]
      
      BDR_res = BDR_res[order(rownames(BDR_res)),]
      UKBBN_res = UKBBN_res[order(rownames(UKBBN_res)),]
      PITT_res = PITT_res[order(rownames(PITT_res)),]
      
      estimate_meta<-cbind(UKBBN_res$Estimate.t1, PITT_res$Estimate.t1)
      stderr_meta<-cbind(UKBBN_res$SE.t1, PITT_res$SE.t1)
      data_meta = cbind(estimate_meta, stderr_meta)
      rownames(data_meta) = rownames(UKBBN_res)
      
      
    } else if(methyl_type == "450k"){
      feature_set = Reduce(intersect, list(rownames(BDR_res),
                                           rownames(UKBBN_res),
                                           rownames(ROSMAP_res),
                                           rownames(PITT_res)))
      
      
      BDR_res = BDR_res[rownames(BDR_res) %in% feature_set,]
      UKBBN_res = UKBBN_res[rownames(UKBBN_res) %in% feature_set,]
      ROSMAP_res = ROSMAP_res[rownames(ROSMAP_res) %in% feature_set,]
      PITT_res = PITT_res[rownames(PITT_res) %in% feature_set,]
      
      BDR_res = BDR_res[order(rownames(BDR_res)),]
      UKBBN_res = UKBBN_res[order(rownames(UKBBN_res)),]
      ROSMAP_res = ROSMAP_res[order(rownames(ROSMAP_res)),]
      PITT_res = PITT_res[order(rownames(PITT_res)),]
      
      estimate_meta<-cbind(UKBBN_res$Estimate.t1, ROSMAP_res$bacon_effectsize, PITT_res$Estimate.t1)
      stderr_meta<-cbind(UKBBN_res$SE.t1, ROSMAP_res$bacon_stderr, PITT_res$SE.t1)
      data_meta = cbind(estimate_meta, stderr_meta)
      rownames(data_meta) = rownames(UKBBN_res)
    }
  }
  
  
  res.meta<-matrix(data=NA, nrow=nrow(data_meta), ncol=8)
  rownames(res.meta)<-as.character(rownames(data_meta))
  colnames(res.meta)<-c("estimate","se","t","pval","ci.lb","ci.ub", "I2", "QEp")
  
  for(i in 1:nrow(data_meta)){
    res<-try(rma.uni(data_meta[i,c(1:(ncol(data_meta)/2))],
                     sei=data_meta[i,c((1+ncol(data_meta)/2):ncol(data_meta))], test = "z", method = "REML",control=list(maxiter=2000)), silent = T)
    res.meta[i,]<-c(as.numeric(coef(summary(res))[c(1:6)]),
                        res$I2,res$QEp)
  }
  
  res.meta = as.data.frame(res.meta)
  
  p.bonf<-p.adjust(res.meta[,4],"bonferroni")
  FC<-2^(abs(res.meta[,1]))
  res.meta<-cbind(res.meta,p.bonf,FC)
  res.meta<-as.data.frame(res.meta)
  
  return(res.meta)
  
}
