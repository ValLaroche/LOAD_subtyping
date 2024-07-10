
limma_DA <- function(log2.cpm, pheno, subtype1, subtype2, cohort, AD = FALSE){
  if(AD == FALSE){ 
    tmppheno = pheno[pheno$subtype %in% c(subtype1, subtype2),]
  } else {
    tmppheno = pheno
  }
  tmppheno$RIN = as.numeric(gsub(",",".", tmppheno$RIN))
  
  tmppheno$RINscale = NA ### Making RIN into simpler intervals
  tmppheno[tmppheno$RIN < summary(tmppheno$RIN)[2],]$RINscale = 1
  tmppheno[tmppheno$RIN >= summary(tmppheno$RIN)[2] & 
                  tmppheno$RIN < summary(tmppheno$RIN)[3],]$RINscale = 2
  tmppheno[tmppheno$RIN >= summary(tmppheno$RIN)[3] & 
                  tmppheno$RIN < summary(tmppheno$RIN)[5],]$RINscale = 3
  tmppheno[tmppheno$RIN >= summary(tmppheno$RIN)[5],]$RINscale = 4
  
  if(cohort == "UKBBN"){
    tmpcounts = log2.cpm[,colnames(log2.cpm) %in% tmppheno$BBNId]
    if(AD == FALSE){ 
      group = factor(tmppheno$subtype, levels = c(subtype1, subtype2))
    } else {
      group = factor(tmppheno$diag)
    }
    tmppheno = calc_sva(tmpcounts, tmppheno, cohort)
    tmppheno$group = group
    design = model.matrix(~group + Gender + Age + RINscale + Brain.Bank + plate +
                            sv1 + sv2 + sv3 + sv4 + sv5, data = tmppheno)
    
  } else if(cohort == "PITT"){
    tmpcounts = log2.cpm[,colnames(log2.cpm) %in% tmppheno$Individual_ID]
    if(AD == FALSE){ 
      group = factor(tmppheno$subtype, levels = c(subtype1, subtype2))
    } else {
      group = factor(tmppheno$diag)
    }    
    tmppheno$group = group
    tmppheno = calc_sva(tmpcounts, tmppheno, cohort)
    
    if(is.null(tmppheno$sv5)){
    design = model.matrix(~group + Sex + Age + RINscale + Plate +
                            sv1 + sv2 + sv3 + sv4, data = tmppheno)
    } else {
      design = model.matrix(~group + Sex + Age + RINscale + Plate + 
                              sv1 + sv2 + sv3 + sv4 + sv5, data = tmppheno)
    }
    
  } else if(cohort == "ROSMAP"){
    tmpcounts = log2.cpm[,colnames(log2.cpm) %in% tmppheno$specimenID]
    if(AD == FALSE){ 
      group = factor(tmppheno$subtype, levels = c(subtype1, subtype2))
    } else {
      group = factor(tmppheno$diag)
    }
    tmppheno$group = group
        
    tmppheno = calc_sva(tmpcounts, tmppheno, cohort)

    design = model.matrix(~group + msex + age_death + RINscale + Sample_Plate +
                            sv1 + sv2 + sv3 + sv4 + sv5, data = tmppheno)
  }
  
# Step 5: Calculate SVs and add to the model along with age, sex,.. if necessary
# Step 6: uisng Limma for DEA

fit <- lmFit(tmpcounts, design)

coeffit = coef(fit)
coeffit = coeffit[,is.na(colSums(coeffit)) == FALSE]

design = design[,colnames(design) %in% colnames(coeffit)]
fit <- lmFit(tmpcounts, design)

tmp <- eBayes(fit)

top.table <- limma::topTable(fit = tmp, n = Inf)
top.table$SE<-(fit$stdev.unscaled*fit$sigma)[,2]

return(top.table)

}

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

meta_analysis <- function(BDR_res = NULL, UKBBN_res, PITT_res, ROSMAP_res = NULL,
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

venn_gen <- function(set1,set2, set3 = NULL, set4 = NULL, pval_threshold = 0.05, FC = FALSE){
  if(FC == TRUE){
    sig_set1 = set1[set1$P.t1 < pval_threshold & set1$effect_size > 0.2,]
    sig_set2 = set2[set2$P.t1 < pval_threshold & set2$effect_size > 0.2,]
    if(!is.null(set3)) {
      sig_set3 = set3[set3$P.t1 < pval_threshold & set3$effect_size > 0.2,]
    } 
    if(!is.null(set4)) {
      sig_set4 = set4[set4$P.t1 < pval_threshold & set4$effect_size > 0.2,]
    }
  } else {
    sig_set1 = set1[set1$P.t1 < pval_threshold,]
    sig_set2 = set2[set2$P.t1 < pval_threshold,]
    if(!is.null(set3)){
      sig_set3 = set3[set3$P.t1 < pval_threshold & set3$effect_size > 0.2,]
    }
    if(!is.null(set4)) {
      sig_set4 = set4[set4$P.t1 < pval_threshold & set4$effect_size > 0.2,]
    }
  }
  
  if(is.null(set3)){
    print(Reduce(intersect, list(rownames(sig_set1),
                           rownames(sig_set2))))
    tovenn <- list("set1" = rownames(sig_set1),
                   "set2" = rownames(sig_set2))
    ggvenn(tovenn, columns = c("set1", "set2"))
  } else if (!is.null(set3) & is.null(set4)) {
    print(Reduce(intersect, list(rownames(sig_set1),
                                 rownames(sig_set2),
                                 rownames(sig_set3))))
    tovenn <- list("set1" = rownames(sig_set1),
                   "set2" = rownames(sig_set2),
                   "set3" = rownames(sig_set3))
    ggvenn(tovenn, columns = c("set1", "set2", "set3"))
  } else {
    print(Reduce(intersect, list(rownames(sig_set1),
                                 rownames(sig_set2),
                                 rownames(sig_set3),
                                 rownames(sig_set4))))
    tovenn <- list("set1" = rownames(set1),
                   "set2" = rownames(set2),
                   "set3" = rownames(set3),
                   "set4" = rownames(set4))
    ggvenn(tovenn, columns = c("set1", "set2", "set3", "set4"))
  }
}

ref_scRNA <- function(){
  NG_scrna = Seurat::Read10X("RNAseq/ref_cib/", gene.column = 1)
  NG_seurat = CreateSeuratObject(counts = NG_scrna ,
                                 min.cells = round(ncol(NG_scrna ) / 100),
                                 min.features = 200,
                                 project = "Nagy")
  
  meta <- strsplit(colnames(NG_seurat), "\\.")
  meta <- lapply(meta, function(x) {
    x <- c(x[1], strsplit(x[2], "_")[[1]])
    return(x)
  })

  meta <- do.call("rbind", meta)
  meta <- as.data.frame(meta)
  colnames(meta) <- c("Celltype", 
                      "Individual", # inferred from there being 17 CTL and 17 MDD patients in the study, and there are 17 levels that are always Control and 17 always Suicide in $Disorder 
                      "Disorder", # note that suicide is equivalent to MDD in this study
                      "Batch", # likely batch given its correspondence across individuals, and that it has six levels which is the number of reported batches. not used by us, so of no importance
                      "Barcode")
  
  NG_seurat$orig.celltype <- meta$Celltype
  NG_seurat$orig.ident <- paste0("Nagy_", meta$Barcode)
  NG_seurat$Individual <- meta$Individual
  NG_seurat$Disorder <- meta$Disorder
  
  keep1 <- which(NG_seurat$Disorder == "Control")
  keep2 <- grep("Mix_", NG_seurat$orig.celltype, invert = TRUE)
  keep <- intersect(keep1, keep2)
  NG_seurat = subset(NG_seurat, cells = keep)
  
  max.depth <- get.max.depth(NG_seurat)
  NG_seurat <- preprocess.fun(NG_seurat, max.depth = max.depth)
  
  
  NG_seurat <- rename.ct(NG_seurat, "Astro", "Astrocytes")
  NG_seurat <- rename.ct(NG_seurat, "Inhib", "Inhibitory")
  NG_seurat <- rename.ct(NG_seurat, "Ex", "Excitatory")
  NG_seurat <- rename.ct(NG_seurat, "Endo", "Endothelia")
  NG_seurat <- rename.ct(NG_seurat, "Micro", "Microglia")
  NG_seurat <- rename.ct(NG_seurat, "Oligo", "Oligodendrocytes")
  NG_seurat <- rename.ct(NG_seurat, "OPC", "OPCs")
  
  ct.counts <- table(NG_seurat$brain.ct)
  
  NG_counts = create.seurat.signature(NG_seurat)
  
  save(NG_counts, file = "C:/Users/p70077107/Desktop/rpkm_cpm_scrna_cibersort.Rdata")
  
  tocib_PITT = data_RNA_PITT
  rownames(tocib_PITT) = gsub("\\.[0-9][0-9]","", rownames(tocib_PITT))
  rownames(tocib_PITT) = gsub("\\.[0-9]","", rownames(tocib_PITT))
  tocib_PITT = tocib_PITT[rownames(tocib_PITT) %in% rownames(NG_counts$cpm),]
  
  tocib_UKBBN = data_RNA_UKBBN
  rownames(tocib_UKBBN) = gsub("\\.[0-9][0-9]","", rownames(tocib_UKBBN))
  rownames(tocib_UKBBN) = gsub("\\.[0-9]","", rownames(tocib_UKBBN))
  tocib_UKBBN = tocib_UKBBN[rownames(tocib_UKBBN) %in% rownames(NG_counts$cpm),]
  
  tocib_ROSMAP = data_RNA_ROSMAP
  rownames(tocib_ROSMAP) = gsub("\\.[0-9][0-9]","", rownames(tocib_ROSMAP))
  rownames(tocib_ROSMAP) = gsub("\\.[0-9]","", rownames(tocib_ROSMAP))
  tocib_ROSMAP = tocib_ROSMAP[rownames(tocib_ROSMAP) %in% rownames(NG_counts$cpm),]
  
  tocib_PITT_red = data_RNA_PITT[,colnames(data_RNA_PITT) %in% pheno_RNA_PITT[pheno_RNA_PITT$subtype == "red",]$Individual_ID]
  tocib_PITT_blue = data_RNA_PITT[,colnames(data_RNA_PITT) %in% pheno_RNA_PITT[pheno_RNA_PITT$subtype == "blue",]$Individual_ID]
  tocib_PITT_ctrl = data_RNA_PITT[,colnames(data_RNA_PITT) %in% pheno_RNA_PITT[pheno_RNA_PITT$subtype == "Control",]$Individual_ID]
  
  tocib_UKBBN_red = data_RNA_UKBBN[,colnames(data_RNA_UKBBN) %in% pheno_RNA_UKBBN[pheno_RNA_UKBBN$subtype == "red",]$BBNId]
  tocib_UKBBN_blue = data_RNA_UKBBN[,colnames(data_RNA_UKBBN) %in% pheno_RNA_UKBBN[pheno_RNA_UKBBN$subtype == "blue",]$BBNId]
  tocib_UKBBN_ctrl = data_RNA_UKBBN[,colnames(data_RNA_UKBBN) %in% pheno_RNA_UKBBN[pheno_RNA_UKBBN$subtype == "Control",]$BBNId]
  
  tocib_ROSMAP_red = data_RNA_ROSMAP[,colnames(data_RNA_ROSMAP) %in% pheno_RNA_ROSMAP[pheno_RNA_ROSMAP$subtype == "red",]$specimenID]
  tocib_ROSMAP_blue = data_RNA_ROSMAP[,colnames(data_RNA_ROSMAP) %in% pheno_RNA_ROSMAP[pheno_RNA_ROSMAP$subtype == "blue",]$specimenID]
  tocib_ROSMAP_ctrl = data_RNA_ROSMAP[,colnames(data_RNA_ROSMAP) %in% pheno_RNA_ROSMAP[pheno_RNA_ROSMAP$subtype == "Control",]$specimenID]
  
  keep = filterByExpr(tocib_PITT_red, min.count = 10, min.prop = 0.75)
  tocib_PITT_red = tocib_PITT_red[keep,]
  tocib_PITT_red = log2(tocib_PITT_red+1)
  keep = filterByExpr(tocib_PITT_blue, min.count = 10, min.prop = 0.75)
  tocib_PITT_blue = tocib_PITT_blue[keep,]
  tocib_PITT_blue = log2(tocib_PITT_blue+1)
  
  keep = filterByExpr(tocib_PITT_ctrl, min.count = 10, min.prop = 0.75)
  tocib_PITT_ctrl = tocib_PITT_ctrl[keep,]
  tocib_PITT_ctrl = log2(tocib_PITT_ctrl+1)
  
  keep = filterByExpr(tocib_UKBBN_red, min.count = 10, min.prop = 0.75)
  tocib_UKBBN_red = tocib_UKBBN_red[keep,]
  tocib_UKBBN_red = log2(tocib_UKBBN_red+1)
  keep = filterByExpr(tocib_UKBBN_blue, min.count = 10, min.prop = 0.75)
  tocib_UKBBN_blue = tocib_UKBBN_blue[keep,]
  tocib_UKBBN_blue = log2(tocib_UKBBN_blue+1)
  keep = filterByExpr(tocib_UKBBN_ctrl, min.count = 10, min.prop = 0.75)
  tocib_UKBBN_ctrl = tocib_UKBBN_ctrl[keep,]
  tocib_UKBBN_ctrl = log2(tocib_UKBBN_ctrl+1)
  
  keep = filterByExpr(tocib_ROSMAP_red, min.count = 10, min.prop = 0.75)
  tocib_ROSMAP_red = tocib_ROSMAP_red[keep,]
  tocib_ROSMAP_red = log2(tocib_ROSMAP_red+1)
  keep = filterByExpr(tocib_ROSMAP_blue, min.count = 10, min.prop = 0.75)
  tocib_ROSMAP_blue = tocib_ROSMAP_blue[keep,]
  tocib_ROSMAP_blue = log2(tocib_ROSMAP_blue+1)
  keep = filterByExpr(tocib_ROSMAP_ctrl, min.count = 10, min.prop = 0.75)
  tocib_ROSMAP_ctrl = tocib_ROSMAP_ctrl[keep,]
  tocib_ROSMAP_ctrl = log2(tocib_ROSMAP_ctrl+1)
  
  write.table(tocib_PITT_red, file = "RNAseq/ref_cib/tocib_PITT_red.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(tocib_UKBBN_red, file = "RNAseq/ref_cib/tocib_UKBBN_red.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(tocib_ROSMAP_red, file = "RNAseq/ref_cib/tocib_ROSMAP_red.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  
  write.table(tocib_PITT_blue, file = "RNAseq/ref_cib/tocib_PITT_blue.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(tocib_UKBBN_blue, file = "RNAseq/ref_cib/tocib_UKBBN_blue.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(tocib_ROSMAP_blue, file = "RNAseq/ref_cib/tocib_ROSMAP_blue.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  
  write.table(tocib_PITT_ctrl, file = "RNAseq/ref_cib/tocib_PITT_ctrl.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(tocib_UKBBN_ctrl, file = "RNAseq/ref_cib/tocib_UKBBN_ctrl.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(tocib_ROSMAP_ctrl, file = "RNAseq/ref_cib/tocib_ROSMAP_ctrl.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  
  
  NG_counts$cpm = NG_counts$cpm[rownames(NG_counts$cpm) %in% rownames(tocib_PITT),]
  NG_counts$cpm = NG_counts$cpm[order(rownames(tocib_PITT)),]
  write.table(NG_counts$cpm, file = "RNAseq/ref_cib/scRNA_cpm_ref_NG.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  
  res_cibersort = list()
  idx = 1
  for(data_subtype in dir("RNAseq/ref_cib/subtype_data/",full.names=TRUE)){
    cibersort_results = CIBERSORT("RNAseq/ref_cib/reference_Nature.txt", data_subtype, perm = 100)  
    res_cibersort[[idx]] = cibersort_results
    idx = idx + 1
  }
  names(res_cibersort) = dir("RNAseq/ref_cib/subtype_data/",full.names=F)
  save(res_cibersort, file = "RNAseq/res_cibersort_all.Rdata")
  
  cib_res_allmelt = data.frame(matrix(nrow = 0, ncol = 5))
  cib_res_all = data.frame(matrix(nrow = 0, ncol = 14))
  idx = 1
  for(cib_res in res_cibersort){
    
    filename = names(res_cibersort)[idx]
    cohort = strsplit(filename, split = "_")[[1]][2]
    subtype = strsplit(filename, split = "_")[[1]][3]
    subtype = strsplit(subtype, split = "\\.")[[1]][1]
    
    cib_res = as.data.frame(cib_res)
    cib_res$ID = rownames(cib_res)
    cib_res$cohort = cohort
    cib_res$subtype = subtype
    
    cibersort_results_melt = reshape2::melt(cib_res[,-c(9,10,11)], id.vars = c("ID","subtype","cohort"))
    cib_res_allmelt = rbind(cib_res_allmelt, cibersort_results_melt)
    cib_res_all = rbind(cib_res_all, cib_res)
    idx = idx + 1
  }
  
  n <- ncol(cib_res)-3
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  colnames(cib_res_allmelt) = c("Sample","subtype", "cohort","Cell_type","Value")
  cib_res_allmelt$Cell_type = as.factor(cib_res_allmelt$Cell_type)
  
  cib_res_allmelt[cib_res_allmelt$subtype == "ctrl",]$subtype = "Control"
  cib_res_allmelt[cib_res_allmelt$subtype == "blue",]$subtype = "Blue"
  cib_res_allmelt[cib_res_allmelt$subtype == "red",]$subtype = "Red"
  cib_res_allmelt$subtype = factor(cib_res_allmelt$subtype, levels = c("Red","Blue","Control"))
  
  plot(ggplot(cib_res_allmelt, aes(x = Sample, y = Value, fill = Cell_type)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values = col_vector) +
         facet_wrap(cohort~subtype, scales = "free") +
         xlab("") + ylab("") + theme_minimal() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  pheno_RNA_PITT = pheno_RNA_PITT[, -grep("cibRNA_", colnames(pheno_RNA_PITT))]
  cib_res_tmp = cib_res_all[cib_res_all$cohort == "PITT",]
  pheno_RNA_PITT = pheno_RNA_PITT[order(pheno_RNA_PITT$Individual_ID),]
  cib_res_tmp = cib_res_tmp[order(rownames(cib_res_tmp)),]
  colnames(cib_res_tmp) = paste0("cibRNA_", colnames(cib_res_tmp))
  pheno_RNA_PITT = cbind(pheno_RNA_PITT, cib_res_tmp)
  
  pheno_RNA_UKBBN = pheno_RNA_UKBBN[, -grep("cibRNA_", colnames(pheno_RNA_UKBBN))]
  cib_res_tmp = cib_res_all[cib_res_all$cohort == "UKBBN",]
  pheno_RNA_UKBBN = pheno_RNA_UKBBN[order(pheno_RNA_UKBBN$BBNId),]
  cib_res_tmp = cib_res_tmp[order(rownames(cib_res_tmp)),]
  colnames(cib_res_tmp) = paste0("cibRNA_", colnames(cib_res_tmp))
  pheno_RNA_UKBBN = cbind(pheno_RNA_UKBBN, cib_res_tmp)
  
  
  pheno_RNA_ROSMAP = pheno_RNA_ROSMAP[, -grep("cibRNA_", colnames(pheno_RNA_ROSMAP))]
  ROSMAP_names = read.table("RNAseq/ROSMAP_samples_oddcelltypes.txt")
  ROSMAP_names = ROSMAP_names$V1
  cib_res_tmp = cib_res_all[cib_res_all$cohort == "ROSMAP",]
  cib_res_tmp = cib_res_tmp[!rownames(cib_res_tmp) %in% ROSMAP_names,]
  pheno_RNA_ROSMAP = pheno_RNA_ROSMAP[order(pheno_RNA_ROSMAP$specimenID),]
  cib_res_tmp = cib_res_tmp[order(rownames(cib_res_tmp)),]
  colnames(cib_res_tmp) = paste0("cibRNA_", colnames(cib_res_tmp))
  pheno_RNA_ROSMAP = pheno_RNA_ROSMAP[!pheno_RNA_ROSMAP$specimenID %in% ROSMAP_names,]
  pheno_RNA_ROSMAP = cbind(pheno_RNA_ROSMAP, cib_res_tmp)
  data_RNA_ROSMAP = data_RNA_ROSMAP[,colnames(data_RNA_ROSMAP) %in% pheno_RNA_ROSMAP$specimenID]
  
  
  #### For Ehsan ####
  
  Ehsan_UKBBN = read.table("RNAseq/raw/gene_count_salmon_10570.txt", header = T)
  Ehsan_PITT = read.table("RNAseq/raw/count_star_10396_hg19.txt", header = T)
  Ehsan_ROSMAP = read.table("RNAseq/raw/gene_count_salmon_ROSMAP.txt", header = T)
  
  load()
  
  
  Ehsan_UKBBN$EnsemblID = gsub("\\.[0-9][0-9]","", Ehsan_UKBBN$EnsemblID)
  Ehsan_UKBBN$EnsemblID = gsub("\\.[0-9]","", Ehsan_UKBBN$EnsemblID)
  colnames(Ehsan_UKBBN) = gsub("X","", colnames(Ehsan_UKBBN))
  write.table(Ehsan_UKBBN, file = "RNAseq/ref_cib/Ehsan_UKBBN.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  
  Ehsan_PITT$EnsemblID = gsub("\\.[0-9][0-9]","", Ehsan_PITT$EnsemblID)
  Ehsan_PITT$EnsemblID = gsub("\\.[0-9]","", Ehsan_PITT$EnsemblID)
  write.table(Ehsan_PITT, file = "RNAseq/ref_cib/Ehsan_PITT.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  
  Ehsan_ROSMAP$EnsemblID = gsub("\\.[0-9][0-9]","", Ehsan_ROSMAP$EnsemblID)
  Ehsan_ROSMAP$EnsemblID = gsub("\\.[0-9]","", Ehsan_ROSMAP$EnsemblID)
  colnames(Ehsan_ROSMAP) = gsub("X","", colnames(Ehsan_ROSMAP))
  Ehsan_ROSMAP = Ehsan_ROSMAP[complete.cases(Ehsan_ROSMAP),]
  Ehsan_ROSMAP = Ehsan_ROSMAP[complete.cases(Ehsan_ROSMAP),]
  Ehsan_ROSMAP1 = Ehsan_ROSMAP[,-249]
  
  
  
  write.table(Ehsan_ROSMAP1, file = "RNAseq/ref_cib/Ehsan_ROSMAP1.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(Ehsan_ROSMAP2, file = "RNAseq/ref_cib/Ehsan_ROSMAP2.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(Ehsan_ROSMAP3, file = "RNAseq/ref_cib/Ehsan_ROSMAP3.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(Ehsan_ROSMAP4, file = "RNAseq/ref_cib/Ehsan_ROSMAP4.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(Ehsan_ROSMAPtest, file = "RNAseq/ref_cib/Ehsan_ROSMAPtest.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  
  cibersort_results_UKBBNall = CIBERSORT("RNAseq/ref_cib/reference_Nature.txt", "RNAseq/ref_cib/Ehsan_UKBBN.txt", perm = 100)  
  cibersort_results_PITTall = CIBERSORT("RNAseq/ref_cib/reference_Nature.txt", "RNAseq/ref_cib/Ehsan_PITT.txt", perm = 100)  
  cibersort_results_ROSMAP1 = CIBERSORT("RNAseq/ref_cib/reference_Nature.txt", "RNAseq/ref_cib/Ehsan_ROSMAP1.txt", perm = 100)  
  cibersort_results_ROSMAP2 = CIBERSORT("RNAseq/ref_cib/reference_Nature.txt", "RNAseq/ref_cib/Ehsan_ROSMAP2.txt", perm = 100)  
  cibersort_results_ROSMAP3 = CIBERSORT("RNAseq/ref_cib/reference_Nature.txt", "RNAseq/ref_cib/Ehsan_ROSMAP3.txt", perm = 100)  
  cibersort_results_ROSMAP4 = CIBERSORT("RNAseq/ref_cib/reference_Nature.txt", "RNAseq/ref_cib/Ehsan_ROSMAP4.txt", perm = 100)  
  
  load("RNAseq/ref_cib/mRNA_ROSMAP_processed_QC_normalized.Rdata")
  write.table(log2.cpm, file ="RNAseq/ref_cib/Ehsan_ROSMAP_red.txt", sep = "\t", row.names = T, col.names = T, quote = F)
  
  cibersort_results_ROSMAPtest = CIBERSORT(sig_matrix = "RNAseq/ref_cib/reference_Nature.txt", mixture_file = "RNAseq/ref_cib/Ehsan_ROSMAP_red.txt", perm = 100)  
  
  
  cibersort_results_PITTall = as.data.frame(cibersort_results_PITTall)
  cibersort_results_PITTall$cohort = "PITT"
  
  cibersort_results_ROSMAPtest = as.data.frame(cibersort_results_ROSMAPtest)
  cibersort_results_ROSMAPtest$cohort = "ROSMAP"
  
  cibersort_results_UKBBNall = as.data.frame(cibersort_results_UKBBNall)
  cibersort_results_UKBBNall$cohort = "UKBBN"
  
  
  ciball = rbind(cibersort_results_PITTall, cibersort_results_ROSMAPtest, cibersort_results_UKBBNall)
  ciball$sample = rownames(ciball)
  
  
  ciball$`P-value` = NULL
  ciball$RMSE = NULL
  ciball$Correlation = NULL
  
  ciball_melt = melt(ciball, id.vars = c("sample","cohort"))
  
  plot(ggplot(ciball_melt, aes(x = sample, y = value, fill = variable)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values = col_vector) +
         facet_wrap(~cohort, scales = "free") +
         xlab("") + ylab("") + theme_minimal() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  
  #####
  
  pheno_RNA_PITT_celltypes = data.frame(pheno_RNA_UKBBN$Neuronal_Inhibitory,
                                        pheno_RNA_UKBBN$Neuronal_Excitatory,
                                        pheno_RNA_UKBBN$Astrocyte,
                                        pheno_RNA_UKBBN$Oligodendrocyte,
                                        pheno_RNA_UKBBN$Microglia,
                                        pheno_RNA_UKBBN$cibRNA_fpkm_astrocytes,
                                        pheno_RNA_UKBBN$cibRNA_fpkm_endothelial,
                                        pheno_RNA_UKBBN$cibRNA_fpkm_fetal_quiescent,
                                        pheno_RNA_UKBBN$cibRNA_fpkm_microglia,
                                        pheno_RNA_UKBBN$cibRNA_fpkm_neurons,
                                        pheno_RNA_UKBBN$cibRNA_fpkm_oligodendrocytes,
                                        pheno_RNA_UKBBN$cibRNA_fpkm_OPC)
  corrplot(cor(pheno_RNA_PITT_celltypes))
  
  
  
  }
                          
get.max.depth <- function(x, max.depth.percentile = 0.995) {
  max.depth <- quantile(x@meta.data$nCount_RNA, probs = max.depth.percentile)
}

preprocess.fun <- function(x, run.downsample = downsample, SCTransform = use.SCTransform, max.depth = max.depth, min.depth = 1000, max.mito = 5, downsample = FALSE, use.SCTransform = FALSE){
  # quantify mitochondrial reads
  x[["percent.mito"]] <- PercentageFeatureSet(object = x, pattern = "^MT-")
  
  # filter to remove outlier nuclei: 
  
  x <- subset(x = x, subset = (nCount_RNA > min.depth) & (nCount_RNA < max.depth) & (percent.mito < max.mito))
  
  # downsample
  if (run.downsample) { x <- downsample.fun(x) }
  
  # normalise expression levels
  x <- NormalizeData(object = x, normalization.method = "LogNormalize", scale.factor = 10000) # standard parameters for Seurat
  
  # find variable genes (i.e. features)
  x <- FindVariableFeatures(object = x, selection.method = "vst", nfeatures = 2000)
  
  
  # further normalisation
  if (use.SCTransform) {
    x <- SCTransform(object = x, vars.to.regress = c("nCount_RNA", "percent.mito")) 
  }
  
  # output
  return(x)
}

rename.ct <- function(x, orig.label, new.label) {
  # create metadata column if not already present
  if (!("brain.ct" %in% colnames(x@meta.data))) x$brain.ct <- "."
  
  # get indices of cells with the original labl
  g <- grep(orig.label, x$orig.celltype)
  
  # rename these
  x$brain.ct[g] <- new.label
  
  return(x)
}

create.seurat.signature <- function(w) {
  # print
  print(w@project.name)
  
  # get counts
  counts <- as.data.frame(w@assays$RNA@counts)
  
  # convert to EnsID and remove non-coding genes
  counts <- addENSID(counts)
  
  # get CPM of every cell
  # cpm <- apply(counts, 2, function(x) {
  #   lib.size <- 10^6 / sum(x)
  #   x <- x * lib.size
  #   return(x)
  # })
  # below is an alternate way, easier on RAM 
  cpm <- counts
  for (j in 1:ncol(cpm)) {
    if (j %% 1000 == 0) print(paste0("Processing sample ", j))
    lib.size <- 10^6 / sum(cpm[,j])
    cpm[,j] <- cpm[,j] * lib.size
  }
  
  # cpm <- as.data.frame(cpm)
  
  # get RPKM of every cell
  rpkm <- length.correct(cpm)
  
  # signatures are the average normalised expression of every member
  output <- list(rpkm = list(), cpm = list())
  for (j in rownames(ct.counts)) {
    # print((j))
    k <- which(w$brain.ct == j)
    
    output$cpm[[j]] <- rowMeans(cpm[,k]) 
    output$rpkm[[j]] <- rowMeans(rpkm[,k]) 
  }
  
  output <- lapply(output, function(x) as.data.frame(do.call("cbind", x)))
  
  # add neurons, which come from pooling exc and inh cells
  neu <- which(w$brain.ct %in% c("Excitatory", "Inhibitory"))
  
  output$cpm$Neurons <- rowMeans(cpm[,neu])
  output$rpkm$Neurons <- rowMeans(rpkm[,neu])
  
  
  # expression threshold: a gene is kept if > 1 unit in at least 1 cell-type
  output <- lapply(output, function(x) {
    keep <- which(apply(x, 1, max, na.rm = TRUE) > 1)
    x <- x[keep,]
    return(x)
  })
  
  # return
  return(output)
}


rpkm <- function(counts) { 
  counts <- counts[which(rownames(counts) %in% rownames(exonicLength)),]  # note 1: counts must have EnsID as rownames   
  
  # Format the exonicLength matrix/list
  length <- transpose(exonicLength)[[1]] # note 2: the file "exonicLength" must be loaded into the global environment,
  names(length) <- rownames(exonicLength)
  
  # Stats
  m <- match(rownames(counts), names(length))
  length <- length[m]/1000
  libsize <- apply(counts, 2, sum) / 10^6
  
  # Normalise for library size, then length
  cpm <- counts 
  for (j in c(1:ncol(cpm))) cpm[,j] <- counts[,j]/libsize[j]
  
  rpkm <- cpm 
  for (j in c(1:ncol(rpkm))) rpkm[,j] <- cpm[,j]/length
  
  return(rpkm)
}

## Convert counts to TPM. There is one note in the function.
tpm <- function(counts) { 
  counts <- counts[which(rownames(counts) %in% rownames(exonicLength)),]  # note 1: counts must have EnsID as rownames   
  
  # Format the exonicLength matrix/list
  length <- transpose(exonicLength)[[1]] # note 2: the file "exonicLength" must be loaded into the global environment,
  names(length) <- rownames(exonicLength)
  
  # Stats on length
  m <- match(rownames(counts), names(length))
  length <- length[m]/1000
  
  # Normalise for length
  length.corrected <- counts
  for (j in c(1:ncol(length.corrected))) length.corrected[,j] <- counts[,j]/length
  
  # Stats on libsize
  libsize <- apply(length.corrected, 2, sum) / 10^6
  
  # Normalise for libsize
  tpm <- length.corrected
  for (j in c(1:ncol(tpm))) tpm[,j] <- length.corrected[,j]/libsize[j]
  
  return(tpm)
}

## Divide expression in a dataframe by length in kilobases
length.correct <- function(exp) {
  exp <- exp[which(rownames(exp) %in% rownames(exonicLength)),]  
  
  length <- transpose(exonicLength)[[1]] # note 2: the file "exonicLength" must be loaded into the global environment,
  names(length) <- rownames(exonicLength)
  
  m <- match(rownames(exp), names(length))
  length <- length[m]/1000
  
  norm <- apply(exp, 2, function(x) x / length)  
  norm <- as.data.frame(norm)
  return(norm)
}


## Convert rownames from EnsID to GeneSymbol. 1 note.
addSymbol <- function(dataframe) {
  annot <- geneInfo[which(geneInfo$Biotype == "protein_coding"),] # Note 1: selects protein-coding genes only! 
  
  sharedAnnotation <- annot[which(annot$ensID%in%rownames(dataframe)),]
  matches <- match(sharedAnnotation$ensID, rownames(dataframe))
  
  # Add gene.symbol to a new column
  dataframe$Gene.Symbol <- "-"
  dataframe$Gene.Symbol[matches] <- sharedAnnotation$Gene.Symbol 
  dataframe <- dataframe[which(dataframe$Gene.Symbol != "-"),]      
  dataframe <- dataframe[!(duplicated(dataframe$Gene.Symbol) | duplicated(dataframe$Gene.Symbol, fromLast = TRUE)), ] # remove deprecated entries
  
  # put gene.symbol into rownames
  rownames(dataframe) <- dataframe$Gene.Symbol
  dataframe <- dataframe[,-which(colnames(dataframe) %in% "Gene.Symbol")]
  return(dataframe)
}  


## Convert rownames from GeneSymbol to EnsID. 1 note.
addENSID <- function(dataframe, pc = TRUE) {
  if (pc) {
    annot <- geneInfo[which(geneInfo$Biotype == "protein_coding"),] # Note 1: selects protein-coding genes only!  
  } else {
    annot <- geneInfo
  }
  
  # Match symbols in dataframe and annotation file
  sharedAnnotation <- annot[which(annot$Gene.Symbol%in%rownames(dataframe)),]
  matches <- match(sharedAnnotation$Gene.Symbol, rownames(dataframe))
  
  # Add ensID to a new column
  dataframe$ensID <- "-"
  dataframe$ensID[matches] <- sharedAnnotation$ensID # A sanity check was performed, which bound m$Approved.Symbol, and it always matched to rownames!
  dataframe <- dataframe[which(dataframe$ensID != "-"),]
  
  # Put ensID into rownames
  rownames(dataframe) <- dataframe$ensID
  dataframe <- dataframe[,-which(colnames(dataframe) %in% "ensID")]
  return(dataframe)
}


## Write a gene expression table with geneID in rownames as a CIBERSORT compatible file
write.CIB <- function(data, dir) {
  data <- cbind(rownames(data), data)
  colnames(data)[1] <- "EnsID"
  write.table(data, file = dir, sep = "\t", quote = FALSE, row.names = FALSE)
}

## Write a standard "phenotype class" table for CIBERSORT. This follows the most basic form of every cell-type being compared to every other cell-type
write.CIB.pheno <- function(data, dir, basic = TRUE) {
  n <- ncol(data)
  
  if (basic) {
    pheno <- matrix(data = 2, nrow = n, ncol = n)
    diag(pheno) <- 1
    pheno <- as.data.frame(pheno)
    rownames(pheno) <- colnames(data)
    
    write.table(pheno, col.names = FALSE, sep = "\t", quote = FALSE, file = dir)  
  }
  
}
