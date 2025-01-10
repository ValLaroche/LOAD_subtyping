setwd("D:/valentin/main/")

# load("RNAseq/current.Rdata")

sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}
source("./RNASeq/limma_SE.R")

########### Methylation ############

load("eQTLDB/BDR_final_set_26-9-23.Rdata")
load("eQTLDB/UKBBN_final_set_27-9-2023.Rdata")
pheno_UKBBN = pheno
load("eQTLDB/PITT_final_set_15-11-23.Rdata")
pheno_PITT = pheno
remove(BDR_ok_ADCtrl, pheno)
load("eQTLDB/ROSMAP_final_set_27-11-23.Rdata")
load("RNAseq/pheno_down_BDR.Rdata")

##### Preparing data #####

load("./betas/BDR_methyl.rdata")
load("./betas/PITT_methyl.rdata")
load("./betas/ROSMAP_methyl.rdata")
load("./betas/UKBBN_methyl.rdata")

##### Generate EWAS
methyl_table_UKBBN_blueC = limma_DA_methyl(beta_values = betas_UKBBN, pheno = pheno_UKBBN, subtype1 = "blue", subtype2 = "Control", cohort = "UKBBN")
methyl_table_UKBBN_redC = limma_DA_methyl(beta_values = betas_UKBBN, pheno = pheno_UKBBN, subtype1 = "red", subtype2 = "Control", cohort = "UKBBN")
methyl_table_UKBBN_redblue = limma_DA_methyl(beta_values = betas_UKBBN, pheno = pheno_UKBBN, subtype1 = "red", subtype2 = "blue", cohort = "UKBBN")
methyl_table_UKBBN_ADC = limma_DA_methyl(beta_values = betas_UKBBN, pheno = pheno_UKBBN, subtype1 = "AD", subtype2 = "Control", cohort = "UKBBN", AD = TRUE)

methyl_table_ROSMAP_blueC = limma_DA_methyl(beta_values = betas_ROSMAP, pheno = pheno_ROSMAP, subtype1 = "blue", subtype2 = "Control", cohort = "ROSMAP")
methyl_table_ROSMAP_redC = limma_DA_methyl(beta_values = betas_ROSMAP, pheno = pheno_ROSMAP, subtype1 = "red", subtype2 = "Control", cohort = "ROSMAP")
methyl_table_ROSMAP_redblue = limma_DA_methyl(beta_values = betas_ROSMAP, pheno = pheno_ROSMAP, subtype1 = "red", subtype2 = "blue", cohort = "ROSMAP")
methyl_table_ROSMAP_ADC = limma_DA_methyl(beta_values = betas_ROSMAP, pheno = pheno_ROSMAP, subtype1 = "AD", subtype2 = "Control", cohort = "ROSMAP", AD = TRUE)

methyl_table_PITT_blueC = limma_DA_methyl(beta_values = betas_PITT, pheno = pheno_PITT, subtype1 = "blue", subtype2 = "Control", cohort = "PITT")
methyl_table_PITT_redC = limma_DA_methyl(beta_values = betas_PITT, pheno = pheno_PITT, subtype1 = "red", subtype2 = "Control", cohort = "PITT")
methyl_table_PITT_redblue = limma_DA_methyl(beta_values = betas_PITT, pheno = pheno_PITT, subtype1 = "red", subtype2 = "blue", cohort = "PITT")
methyl_table_PITT_ADC = limma_DA_methyl(beta_values = betas_PITT, pheno = pheno_PITT, subtype1 = "AD", subtype2 = "Control", cohort = "PITT", AD = TRUE)

##### Calculate Cohen's D for effect size

cohend_UKBBN_blueC = pbapply::pbapply(betas_UKBBN, 1, function(x){cohenD(cpg = x, pheno = pheno_UKBBN, groups = c("blue", "Control"), AD = FALSE)})
cohend_UKBBN_redC = pbapply::pbapply(betas_UKBBN, 1, function(x){cohenD(cpg = x, pheno = pheno_UKBBN, groups = c("red", "Control"), AD = FALSE)})
cohend_UKBBN_redblueC = pbapply::pbapply(betas_UKBBN, 1, function(x){cohenD(cpg = x, pheno = pheno_UKBBN, groups = c("red", "blue"), AD = FALSE)})
cohend_UKBBN_ADC = pbapply::pbapply(betas_UKBBN, 1, function(x){cohenD(cpg = x, pheno = pheno_UKBBN, groups = c("AD", "Control"), AD = TRUE)})

cohend_ROSMAP_blueC = pbapply::pbapply(betas_ROSMAP, 1, function(x){cohenD(cpg = x, pheno = pheno_ROSMAP, groups = c("blue", "Control"), AD = FALSE)})
cohend_ROSMAP_redC = pbapply::pbapply(betas_ROSMAP, 1, function(x){cohenD(cpg = x, pheno = pheno_ROSMAP, groups = c("red", "Control"), AD = FALSE)})
cohend_ROSMAP_redblueC = pbapply::pbapply(betas_ROSMAP, 1, function(x){cohenD(cpg = x, pheno = pheno_ROSMAP, groups = c("red", "blue"), AD = FALSE)})
cohend_ROSMAP_ADC = pbapply::pbapply(betas_ROSMAP, 1, function(x){cohenD(cpg = x, pheno = pheno_ROSMAP, groups = c("AD", "Control"), AD = TRUE)})

cohend_PITT_blueC = pbapply::pbapply(betas_PITT, 1, function(x){cohenD(cpg = x, pheno = pheno_PITT, groups = c("blue", "Control"), AD = FALSE)})
cohend_PITT_redC = pbapply::pbapply(betas_PITT, 1, function(x){cohenD(cpg = x, pheno = pheno_PITT, groups = c("red", "Control"), AD = FALSE)})
cohend_PITT_redblueC = pbapply::pbapply(betas_PITT, 1, function(x){cohenD(cpg = x, pheno = pheno_PITT, groups = c("red", "blue"), AD = FALSE)})
cohend_PITT_ADC = pbapply::pbapply(betas_PITT, 1, function(x){cohenD(cpg = x, pheno = pheno_PITT, groups = c("AD", "Control"), AD = TRUE)})

head(cohend_UKBBN_ADC)

##### Turn into dataframes

methyl_table_UKBBN_blueC = as.data.frame(methyl_table_UKBBN_blueC)
methyl_table_UKBBN_redC = as.data.frame(methyl_table_UKBBN_redC)
methyl_table_UKBBN_redblue = as.data.frame(methyl_table_UKBBN_redblue)
methyl_table_UKBBN_ADC = as.data.frame(methyl_table_UKBBN_ADC)

methyl_table_ROSMAP_blueC = as.data.frame(methyl_table_ROSMAP_blueC)
methyl_table_ROSMAP_redC = as.data.frame(methyl_table_ROSMAP_redC)
methyl_table_ROSMAP_redblue = as.data.frame(methyl_table_ROSMAP_redblue)
methyl_table_ROSMAP_ADC = as.data.frame(methyl_table_ROSMAP_ADC)

methyl_table_PITT_blueC = as.data.frame(methyl_table_PITT_blueC)
methyl_table_PITT_redC = as.data.frame(methyl_table_PITT_redC)
methyl_table_PITT_redblue = as.data.frame(methyl_table_PITT_redblue)
methyl_table_PITT_ADC = as.data.frame(methyl_table_PITT_ADC)

##### Add Cohen's D to results of EWAS

methyl_table_UKBBN_blueC = methyl_table_UKBBN_blueC[order(rownames(methyl_table_UKBBN_blueC)),]
methyl_table_UKBBN_blueC$effect_size = cohend_UKBBN_blueC
methyl_table_UKBBN_redC = methyl_table_UKBBN_redC[order(rownames(methyl_table_UKBBN_redC)),]
methyl_table_UKBBN_redC$effect_size = cohend_UKBBN_redC
methyl_table_UKBBN_redblue = methyl_table_UKBBN_redblue[order(rownames(methyl_table_UKBBN_redblue)),]
methyl_table_UKBBN_redblue$effect_size = cohend_UKBBN_redblueC
methyl_table_UKBBN_ADC = methyl_table_UKBBN_ADC[order(rownames(methyl_table_UKBBN_ADC)),]
methyl_table_UKBBN_ADC$effect_size = cohend_UKBBN_ADC

methyl_table_ROSMAP_blueC = methyl_table_ROSMAP_blueC[order(rownames(methyl_table_ROSMAP_blueC)),]
methyl_table_ROSMAP_blueC$effect_size = cohend_ROSMAP_blueC
methyl_table_ROSMAP_redC = methyl_table_ROSMAP_redC[order(rownames(methyl_table_ROSMAP_redC)),]
methyl_table_ROSMAP_redC$effect_size = cohend_ROSMAP_redC
methyl_table_ROSMAP_redblue = methyl_table_ROSMAP_redblue[order(rownames(methyl_table_ROSMAP_redblue)),]
methyl_table_ROSMAP_redblue$effect_size = cohend_ROSMAP_redblueC
methyl_table_ROSMAP_ADC = methyl_table_ROSMAP_ADC[order(rownames(methyl_table_ROSMAP_ADC)),]
methyl_table_ROSMAP_ADC$effect_size = cohend_ROSMAP_ADC

methyl_table_PITT_blueC = methyl_table_PITT_blueC[order(rownames(methyl_table_PITT_blueC)),]
methyl_table_PITT_blueC$effect_size = cohend_PITT_blueC
methyl_table_PITT_redC = methyl_table_PITT_redC[order(rownames(methyl_table_PITT_redC)),]
methyl_table_PITT_redC$effect_size = cohend_PITT_redC
methyl_table_PITT_redblue = methyl_table_PITT_redblue[order(rownames(methyl_table_PITT_redblue)),]
methyl_table_PITT_redblue$effect_size = cohend_PITT_redblueC
methyl_table_PITT_ADC = methyl_table_PITT_ADC[order(rownames(methyl_table_PITT_ADC)),]
methyl_table_PITT_ADC$effect_size = cohend_PITT_ADC

##### Test for inflation and calculate new pvalues

methyl_table_UKBBN_blueC = inflation_control(methyl_table_UKBBN_blueC, subtype2 = "Control")
methyl_table_UKBBN_redC = inflation_control(methyl_table_UKBBN_redC, subtype2 = "Control")
methyl_table_UKBBN_redblue = inflation_control(methyl_table_UKBBN_redblue, subtype2 = "blue")
methyl_table_UKBBN_ADC = inflation_control(methyl_table_UKBBN_ADC, subtype2 = "Control")

methyl_table_ROSMAP_blueC = inflation_control(methyl_table_ROSMAP_blueC, subtype2 = "Control")
methyl_table_ROSMAP_redC = inflation_control(methyl_table_ROSMAP_redC, subtype2 = "Control")
methyl_table_ROSMAP_redblue = inflation_control(methyl_table_ROSMAP_redblue, subtype2 = "blue")
methyl_table_ROSMAP_ADC = inflation_control(methyl_table_ROSMAP_ADC, subtype2 = "Control")

methyl_table_PITT_blueC = inflation_control(methyl_table_PITT_blueC, subtype2 = "Control")
methyl_table_PITT_redC = inflation_control(methyl_table_PITT_redC, subtype2 = "Control")
methyl_table_PITT_redblue = inflation_control(methyl_table_PITT_redblue, subtype2 = "blue")
methyl_table_PITT_ADC = inflation_control(methyl_table_PITT_ADC, subtype2 = "Control")

##### Perform meta analysis of the 3 cohorts

meta_methyl_450k_blueC = meta_analysis(UKBBN_res = methyl_table_UKBBN_blueC, 
                                       PITT_res = methyl_table_PITT_blueC, 
                                       ROSMAP_res = methyl_table_ROSMAP_blueC,
                                       data_type = "Methyl", methyl_type = "450k")
meta_methyl_450k_redC = meta_analysis(UKBBN_res = methyl_table_UKBBN_redC, 
                                       PITT_res = methyl_table_PITT_redC, 
                                       ROSMAP_res = methyl_table_ROSMAP_redC,
                                       data_type = "Methyl", methyl_type = "450k")
meta_methyl_450k_redblue = meta_analysis(UKBBN_res = methyl_table_UKBBN_redblue, 
                                       PITT_res = methyl_table_PITT_redblue, 
                                       ROSMAP_res = methyl_table_ROSMAP_redblue,
                                       data_type = "Methyl", methyl_type = "450k")
meta_methyl_450k_ADC = meta_analysis(UKBBN_res = methyl_table_UKBBN_ADC, 
                                       PITT_res = methyl_table_PITT_ADC, 
                                       ROSMAP_res = methyl_table_ROSMAP_ADC,
                                       data_type = "Methyl", methyl_type = "450k")


# Take only significant effect sizes across cohorts
cpgs_intersect = Reduce(intersect, list(rownames(methyl_table_UKBBN_blueC),
                                        rownames(methyl_table_PITT_blueC),
                                        rownames(methyl_table_ROSMAP_blueC)))

cohenD_cumul_450k_blueC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_blueC[rownames(methyl_table_UKBBN_blueC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_PITT" = methyl_table_PITT_blueC[rownames(methyl_table_PITT_blueC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_ROSMAP" = methyl_table_ROSMAP_blueC[rownames(methyl_table_ROSMAP_blueC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_450k_blueC) = cpgs_intersect
cohenD_cumul_450k_blueC$sig_expr = rowSums(abs(cohenD_cumul_450k_blueC) > 0.2) > 0
cohenD_cumul_450k_blueC$dir_expr = (cohenD_cumul_450k_blueC$cohens_D_UKBBN > 0 & cohenD_cumul_450k_blueC$cohens_D_ROSMAP > 0 & cohenD_cumul_450k_blueC$cohens_D_PITT > 0) | 
  (cohenD_cumul_450k_blueC$cohens_D_UKBBN < 0 & cohenD_cumul_450k_blueC$cohens_D_PITT < 0 & cohenD_cumul_450k_blueC$cohens_D_ROSMAP < 0) 

cohenD_cumul_450k_redC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_redC[rownames(methyl_table_UKBBN_redC) %in% cpgs_intersect,]$effect_size,
                                    "cohens_D_PITT" = methyl_table_PITT_redC[rownames(methyl_table_PITT_redC) %in% cpgs_intersect,]$effect_size,
                                    "cohens_D_ROSMAP" = methyl_table_ROSMAP_redC[rownames(methyl_table_ROSMAP_redC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_450k_redC) = cpgs_intersect
cohenD_cumul_450k_redC$sig_expr = rowSums(abs(cohenD_cumul_450k_redC) > 0.2) > 0
cohenD_cumul_450k_redC$dir_expr = (cohenD_cumul_450k_redC$cohens_D_UKBBN > 0 & cohenD_cumul_450k_redC$cohens_D_ROSMAP > 0 & cohenD_cumul_450k_redC$cohens_D_PITT > 0) | 
  (cohenD_cumul_450k_redC$cohens_D_UKBBN < 0 & cohenD_cumul_450k_redC$cohens_D_PITT < 0 & cohenD_cumul_450k_redC$cohens_D_ROSMAP < 0) 

cohenD_cumul_450k_redblue = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_redblue[rownames(methyl_table_UKBBN_redblue) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_PITT" = methyl_table_PITT_redblue[rownames(methyl_table_PITT_redblue) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_ROSMAP" = methyl_table_ROSMAP_redblue[rownames(methyl_table_ROSMAP_redblue) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_450k_redblue) = cpgs_intersect
cohenD_cumul_450k_redblue$sig_expr = rowSums(abs(cohenD_cumul_450k_redblue) > 0.2) > 0
cohenD_cumul_450k_redblue$dir_expr = (cohenD_cumul_450k_redblue$cohens_D_UKBBN > 0 & cohenD_cumul_450k_redblue$cohens_D_ROSMAP > 0 & cohenD_cumul_450k_redblue$cohens_D_PITT > 0) | 
  (cohenD_cumul_450k_redblue$cohens_D_UKBBN < 0 & cohenD_cumul_450k_redblue$cohens_D_PITT < 0 & cohenD_cumul_450k_redblue$cohens_D_ROSMAP < 0) 

cohenD_cumul_450k_ADC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_ADC[rownames(methyl_table_UKBBN_ADC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_PITT" = methyl_table_PITT_ADC[rownames(methyl_table_PITT_ADC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_ROSMAP" = methyl_table_ROSMAP_ADC[rownames(methyl_table_ROSMAP_ADC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_450k_ADC) = cpgs_intersect
cohenD_cumul_450k_ADC$sig_expr = rowSums(abs(cohenD_cumul_450k_ADC) > 0.2) > 0
cohenD_cumul_450k_ADC$dir_expr = (cohenD_cumul_450k_ADC$cohens_D_UKBBN > 0 & cohenD_cumul_450k_ADC$cohens_D_ROSMAP > 0 & cohenD_cumul_450k_ADC$cohens_D_PITT > 0) | 
  (cohenD_cumul_450k_ADC$cohens_D_UKBBN < 0 & cohenD_cumul_450k_ADC$cohens_D_PITT < 0 & cohenD_cumul_450k_ADC$cohens_D_ROSMAP < 0) 


##### Perform meta analysis at the EPIC level (Results not used because non-trusted)

meta_methyl_EPIC_blueC = meta_analysis(UKBBN_res = methyl_table_UKBBN_blueC, 
                                       PITT_res = methyl_table_PITT_blueC, 
                                       data_type = "Methyl", methyl_type = "EPIC")
meta_methyl_EPIC_redC = meta_analysis(UKBBN_res = methyl_table_UKBBN_redC, 
                                      PITT_res = methyl_table_PITT_redC,
                                      data_type = "Methyl", methyl_type = "EPIC")
meta_methyl_EPIC_redblue = meta_analysis(UKBBN_res = methyl_table_UKBBN_redblue, 
                                         PITT_res = methyl_table_PITT_redblue, 
                                         data_type = "Methyl", methyl_type = "EPIC")
meta_methyl_EPIC_ADC = meta_analysis(UKBBN_res = methyl_table_UKBBN_ADC, 
                                     PITT_res = methyl_table_PITT_ADC, 
                                     data_type = "Methyl", methyl_type = "EPIC")


# Take only significant effect sizes across cohorts
cpgs_intersect = Reduce(intersect, list(rownames(methyl_table_UKBBN_blueC),
                                        rownames(methyl_table_PITT_blueC)))
cohenD_cumul_EPIC_blueC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_blueC[rownames(methyl_table_UKBBN_blueC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_PITT" = methyl_table_PITT_blueC[rownames(methyl_table_PITT_blueC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_EPIC_blueC) = cpgs_intersect
cohenD_cumul_EPIC_blueC$sig_expr = rowSums(abs(cohenD_cumul_EPIC_blueC)> 0.2) > 0
cohenD_cumul_EPIC_blueC$dir_expr = (cohenD_cumul_EPIC_blueC$cohens_D_UKBBN > 0 & cohenD_cumul_EPIC_blueC$cohens_D_PITT > 0) | 
  (cohenD_cumul_EPIC_blueC$cohens_D_UKBBN < 0 & cohenD_cumul_EPIC_blueC$cohens_D_PITT < 0) 


cohenD_cumul_EPIC_redC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_redC[rownames(methyl_table_UKBBN_redC) %in% cpgs_intersect,]$effect_size,
                                    "cohens_D_PITT" = methyl_table_PITT_redC[rownames(methyl_table_PITT_redC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_EPIC_redC) = cpgs_intersect
cohenD_cumul_EPIC_redC$sig_expr = rowSums(abs(cohenD_cumul_EPIC_redC)> 0.2) > 0
cohenD_cumul_EPIC_redC$dir_expr = (cohenD_cumul_EPIC_redC$cohens_D_UKBBN > 0 & cohenD_cumul_EPIC_redC$cohens_D_PITT > 0) | 
  (cohenD_cumul_EPIC_redC$cohens_D_UKBBN < 0 & cohenD_cumul_EPIC_redC$cohens_D_PITT < 0) 

cohenD_cumul_EPIC_redblue = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_redblue[rownames(methyl_table_UKBBN_redblue) %in% cpgs_intersect,]$effect_size,
                                       "cohens_D_PITT" = methyl_table_PITT_redblue[rownames(methyl_table_PITT_redblue) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_EPIC_redblue) = cpgs_intersect
cohenD_cumul_EPIC_redblue$sig_expr = rowSums(abs(cohenD_cumul_EPIC_redblue)> 0.2) > 0
cohenD_cumul_EPIC_redblue$dir_expr = (cohenD_cumul_EPIC_redblue$cohens_D_UKBBN > 0 & cohenD_cumul_EPIC_redblue$cohens_D_PITT > 0) | 
  (cohenD_cumul_EPIC_redblue$cohens_D_UKBBN < 0 & cohenD_cumul_EPIC_redblue$cohens_D_PITT < 0) 

cohenD_cumul_EPIC_ADC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_ADC[rownames(methyl_table_UKBBN_ADC) %in% cpgs_intersect,]$effect_size,
                                   "cohens_D_PITT" = methyl_table_PITT_ADC[rownames(methyl_table_PITT_ADC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_EPIC_ADC) = cpgs_intersect
cohenD_cumul_EPIC_ADC$sig_expr = rowSums(abs(cohenD_cumul_EPIC_ADC)> 0.2) > 0
cohenD_cumul_EPIC_ADC$dir_expr = (cohenD_cumul_EPIC_ADC$cohens_D_UKBBN > 0 & cohenD_cumul_EPIC_ADC$cohens_D_PITT > 0) | 
  (cohenD_cumul_EPIC_ADC$cohens_D_UKBBN < 0 & cohenD_cumul_EPIC_ADC$cohens_D_PITT < 0) 

#####

# Filter for pvalue
pvalue_threshold = 1e-3
sig_cohenD = rownames(cohenD_cumul_450k_blueC[cohenD_cumul_450k_blueC$dir_expr == TRUE &
                                                cohenD_cumul_450k_blueC$sig_expr == TRUE,])

sig_methyl_450k_blueC = meta_methyl_450k_blueC[rownames(meta_methyl_450k_blueC) %in% sig_cohenD,]
sig_methyl_450k_blueC = sig_methyl_450k_blueC[sig_methyl_450k_blueC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_450k_blueC[rownames(cohenD_cumul_450k_blueC) %in% rownames(sig_methyl_450k_blueC),]
sig_methyl_450k_blueC = cbind(sig_methyl_450k_blueC, df_cohenD)
sig_methyl_450k_blueC = sig_methyl_450k_blueC[order(sig_methyl_450k_blueC$pval, decreasing = T),]

#
sig_cohenD = rownames(cohenD_cumul_450k_redC[cohenD_cumul_450k_redC$dir_expr == TRUE &
                                               cohenD_cumul_450k_redC$sig_expr == TRUE,])

sig_methyl_450k_redC = meta_methyl_450k_redC[rownames(meta_methyl_450k_redC) %in% sig_cohenD,]
sig_methyl_450k_redC = sig_methyl_450k_redC[sig_methyl_450k_redC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_450k_redC[rownames(cohenD_cumul_450k_redC) %in% rownames(sig_methyl_450k_redC),]
sig_methyl_450k_redC = cbind(sig_methyl_450k_redC, df_cohenD)
sig_methyl_450k_redC = sig_methyl_450k_redC[order(sig_methyl_450k_redC$pval, decreasing = T),]

#
sig_cohenD = rownames(cohenD_cumul_450k_redblue[cohenD_cumul_450k_redblue$dir_expr == TRUE &
                                                  cohenD_cumul_450k_redblue$sig_expr == TRUE,])

sig_methyl_450k_redblue = meta_methyl_450k_redblue[rownames(meta_methyl_450k_redblue) %in% sig_cohenD,]
sig_methyl_450k_redblue = sig_methyl_450k_redblue[sig_methyl_450k_redblue$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_450k_redblue[rownames(cohenD_cumul_450k_redblue) %in% rownames(sig_methyl_450k_redblue),]
sig_methyl_450k_redblue = cbind(sig_methyl_450k_redblue, df_cohenD)
sig_methyl_450k_redblue = sig_methyl_450k_redblue[order(sig_methyl_450k_redblue$pval, decreasing = T),]

#
sig_cohenD = rownames(cohenD_cumul_450k_ADC[cohenD_cumul_450k_ADC$dir_expr == TRUE &
                                              cohenD_cumul_450k_ADC$sig_expr == TRUE,])

sig_methyl_450k_ADC = meta_methyl_450k_ADC[rownames(meta_methyl_450k_ADC) %in% sig_cohenD,]
sig_methyl_450k_ADC = sig_methyl_450k_ADC[sig_methyl_450k_ADC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_450k_ADC[rownames(cohenD_cumul_450k_ADC) %in% rownames(sig_methyl_450k_ADC),]
sig_methyl_450k_ADC = cbind(sig_methyl_450k_ADC, df_cohenD)
sig_methyl_450k_ADC = sig_methyl_450k_ADC[order(sig_methyl_450k_ADC$pval, decreasing = T),]


tovenn <- list("BluevCtrl" = rownames(sig_methyl_450k_blueC),
               "RedvCtrl" = rownames(sig_methyl_450k_redC),
               "RedvBlue" = rownames(sig_methyl_450k_redblue))
ggvenn(tovenn, columns = c("BluevCtrl", "RedvCtrl", "RedvBlue")) 

tovenn <- list("BLUE vs Control" = rownames(sig_methyl_450k_blueC),
               "RED vs Control" = rownames(sig_methyl_450k_redC),
               "AD vs Control" = rownames(sig_methyl_450k_ADC))
ggvenn(tovenn, columns = c("BLUE vs Control", "RED vs Control", "AD vs Control")) 

# Full venn diagram for significant CpGs
venn <- Venn(tovenn)
d <- process_data(venn)


colorGroups <- c("BLUE vs Control"  = 'steelblue1', "RED vs Control"  = 'brown2', "AD vs Control"  = 'grey')        
# use colorRampPalette to create function that interpolates colors 
colfunc <- colorRampPalette(colorGroups)
col = colfunc(7)

dregion = venn_regionlabel(d)
dregion$perc = paste0("(",format((dregion$count / sum(dregion$count))*100, digits = 1), " %)")
dregion$Y[5] = -2.29
dregion$Y[6] = -2.29
dregion$percY = dregion$Y - 0.52
dregion$percX = c(-1.65065776, 5.6065745,1.9999951,1.9499941,-0.305876,4.2505961,1.9400074)

dgroups = venn_setlabel(d)
dgroups$X[1] = -2.8
dgroups$X[2] = 6.8
pdf("SUBTYPING-PAPER-TABLEFIG/VENN-ALL-EWAS.pdf", width = 9, height = 9)
ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = id, group = id), 
               data = venn_regionedge(d)) +
  # 2. set edge layer
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
            data = dgroups, size = 7) +
  # 4. region label layer
  geom_text(aes(X, Y, label = count), 
            data = dregion, size = 6) +
  geom_text(aes(percX, percY, label = perc), 
            data = dregion, size = 6) +
  coord_equal() +
  scale_fill_manual(values = col, breaks = c("1","1/3", "1/2", "2", "2/3", "1/2/3", "3")) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

# Filter for hypo/hypermethylation
sig_methyl_450k_blueC_neg = sig_methyl_450k_blueC[sig_methyl_450k_blueC$cohens_D_UKBBN < 0,]
sig_methyl_450k_blueC_pos = sig_methyl_450k_blueC[sig_methyl_450k_blueC$cohens_D_UKBBN >= 0,]

sig_methyl_450k_redC_neg = sig_methyl_450k_redC[sig_methyl_450k_redC$cohens_D_UKBBN < 0,]
sig_methyl_450k_redC_pos = sig_methyl_450k_redC[sig_methyl_450k_redC$cohens_D_UKBBN >= 0,]

sig_methyl_450k_ADC_neg = sig_methyl_450k_ADC[sig_methyl_450k_ADC$cohens_D_UKBBN < 0,]
sig_methyl_450k_ADC_pos = sig_methyl_450k_ADC[sig_methyl_450k_ADC$cohens_D_UKBBN >= 0,]

sig_methyl_450k_redblue_neg = sig_methyl_450k_redblue[sig_methyl_450k_redblue$cohens_D_UKBBN < 0,]
sig_methyl_450k_redblue_pos = sig_methyl_450k_redblue[sig_methyl_450k_redblue$cohens_D_UKBBN >= 0,]

tovenn <- list("Blue vs Control" = rownames(sig_methyl_450k_blueC_neg),
               "Red vs Control" = rownames(sig_methyl_450k_redC_neg),
               "AD vs Control" = rownames(sig_methyl_450k_ADC_neg))
pdf("SUBTYPING-PAPER-TABLEFIG/VENN-ALL-EWASneg.pdf", width = 8, height = 8)
ggvenn(tovenn, columns = c("Blue vs Control", "Red vs Control", "AD vs Control"),
       fill_color = c("blue","red", "grey"),
       stroke_color = "white", stroke_size = 0.00001,
       text_size = 6, set_name_size = 7) 
dev.off()

tovenn <- list("Blue vs Control" = rownames(sig_methyl_450k_blueC_pos),
               "Red vs Control" = rownames(sig_methyl_450k_redC_pos),
               "AD vs Control" = rownames(sig_methyl_450k_ADC_pos))
pdf("SUBTYPING-PAPER-TABLEFIG/VENN-ALL-EWASpos.pdf", width = 8, height = 8)
ggvenn(tovenn, columns = c("Blue vs Control", "Red vs Control", "AD vs Control"),
       fill_color = c("steelblue1","brown2", "grey"),
       stroke_color = "white", stroke_size = 0.00001,
       text_size = 6, set_name_size = 7) 
dev.off()

#####

tovenn <- list("EPIC_blueC" = rownames(sig_methyl_EPIC_blueC),
               "450k_blueC" = rownames(sig_methyl_450k_blueC))
ggvenn(tovenn, columns = c("EPIC_blueC", "450k_blueC")) 

tovenn <- list("EPIC_redC" = rownames(sig_methyl_EPIC_redC),
               "450k_redC" = rownames(sig_methyl_450k_redC))
ggvenn(tovenn, columns = c("EPIC_redC", "450k_redC")) 

tovenn <- list("EPIC_redblue" = rownames(sig_methyl_EPIC_redblue),
               "450k_redblue" = rownames(sig_methyl_450k_redblue))
ggvenn(tovenn, columns = c("EPIC_redblue", "450k_redblue")) 


overlap450 = rownames(sig_methyl_EPIC_ADC)[rownames(sig_methyl_EPIC_ADC) %in% rownames(betas_ROSMAP)]
tovenn <- list("EPIC_ADC" = rownames(sig_methyl_EPIC_ADC),
               "450k_ADC" = rownames(sig_methyl_450k_ADC))
ggvenn(tovenn, columns = c("EPIC_ADC", "450k_ADC")) 


#####

# Look at the intersect with previous meta-analysis
inter_red_blue = intersect(rownames(sig_methyl_450k_blueC), rownames(sig_methyl_450k_redC))
sig_methyl_450k_blueC = sig_methyl_450k_blueC[!rownames(sig_methyl_450k_blueC) %in% inter_red_blue,]
sig_methyl_450k_redC = sig_methyl_450k_redC[!rownames(sig_methyl_450k_redC) %in% inter_red_blue,]
CpGs_previous = read.table("AD_BH_CpGs.txt")

intersect(CpGs_previous$V1, rownames(sig_methyl_450k_ADC))
intersect(CpGs_previous$V1, rownames(sig_methyl_450k_blueC))
intersect(CpGs_previous$V1, rownames(sig_methyl_450k_redblue))
intersect(CpGs_previous$V1, rownames(sig_methyl_450k_redC))

Reduce(intersect, list(CpGs_previous$V1, rownames(sig_methyl_450k_blueC),
                       rownames(sig_methyl_450k_redC)))
Reduce(intersect, list(CpGs_previous$V1, rownames(sig_methyl_450k_redC),
                       rownames(sig_methyl_450k_ADC)))

# Control for inflation after meta-analysis
gctrl = gap::gcontrol2(meta_methyl_450k_ADC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_450k_blueC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_450k_redC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_450k_redblue$pval)
print(gctrl$lambda)

# Check for APOE CpG
meta_methyl_450k_ADC[rownames(meta_methyl_450k_ADC) == "cg05066959",]
meta_methyl_450k_redC[rownames(meta_methyl_450k_redC) == "cg05066959",]
meta_methyl_450k_blueC[rownames(meta_methyl_450k_blueC) == "cg05066959",]

methyl_table_UKBBN_ADC[rownames(methyl_table_UKBBN_ADC) == "cg05066959",]
methyl_table_PITT_ADC[rownames(methyl_table_PITT_ADC) == "cg05066959",]
methyl_table_ROSMAP_ADC[rownames(methyl_table_ROSMAP_ADC) == "cg05066959",]

##### gometh 450k #####
gometh_res_450k_blueC = gometh(sig.cpg = rownames(sig_methyl_450k_blueC),
                                   all.cpg = rownames(meta_methyl_450k_ADC),
                                   collection = c("GO", "KEGG"),
                                   array.type = c("450K", "EPIC"),
                                   sig.genes = FALSE)
gometh_res_450k_redC = gometh(sig.cpg = rownames(sig_methyl_450k_redC),
                                  all.cpg = rownames(meta_methyl_450k_ADC),
                                  collection = c("GO", "KEGG"),
                                  array.type = c("450K", "EPIC"),
                                  sig.genes = FALSE)
gometh_res_450k_redblue = gometh(sig.cpg = rownames(sig_methyl_450k_redblue),
                                     all.cpg = rownames(meta_methyl_450k_ADC),
                                     collection = c("GO", "KEGG"),
                                     array.type = c("450K", "EPIC"),
                                     sig.genes = FALSE)
gometh_res_450k_ADC = gometh(sig.cpg = rownames(sig_methyl_450k_ADC),
                                 all.cpg = rownames(meta_methyl_450k_ADC),
                                 collection = c("GO", "KEGG"),
                                 array.type = c("450K", "EPIC"),
                                 sig.genes = FALSE)
##### gometh 450k pos #####
gometh_res_450k_blueC_pos = gometh(sig.cpg = rownames(sig_methyl_450k_blueC_pos),
                                   all.cpg = rownames(meta_methyl_450k_ADC),
                                   collection = c("GO", "KEGG"),
                                   array.type = c("450K", "EPIC"),
                                   sig.genes = FALSE)
gometh_res_450k_redC_pos = gometh(sig.cpg = rownames(sig_methyl_450k_redC_pos),
                                  all.cpg = rownames(meta_methyl_450k_ADC),
                                  collection = c("GO", "KEGG"),
                                  array.type = c("450K", "EPIC"),
                                  sig.genes = FALSE)
gometh_res_450k_redblue_pos = gometh(sig.cpg = rownames(sig_methyl_450k_redblue_pos),
                                     all.cpg = rownames(meta_methyl_450k_ADC),
                                     collection = c("GO", "KEGG"),
                                     array.type = c("450K", "EPIC"),
                                     sig.genes = FALSE)
gometh_res_450k_ADC_pos = gometh(sig.cpg = rownames(sig_methyl_450k_ADC_pos),
                                 all.cpg = rownames(meta_methyl_450k_ADC),
                                 collection = c("GO", "KEGG"),
                                 array.type = c("450K", "EPIC"),
                                 sig.genes = FALSE)
##### gometh EPIC neg #####
gometh_res_450k_blueC_neg = gometh(sig.cpg = rownames(sig_methyl_450k_blueC_neg),
                               all.cpg = rownames(meta_methyl_450k_ADC),
                               collection = c("GO", "KEGG"),
                               array.type = c("450K", "EPIC"),
                               sig.genes = FALSE)
gometh_res_450k_redC_neg = gometh(sig.cpg = rownames(sig_methyl_450k_redC_neg),
                              all.cpg = rownames(meta_methyl_450k_ADC),
                              collection = c("GO", "KEGG"),
                              array.type = c("450K", "EPIC"),
                              sig.genes = FALSE)
gometh_res_450k_redblue_neg = gometh(sig.cpg = rownames(sig_methyl_450k_redblue_neg),
                                 all.cpg = rownames(meta_methyl_450k_ADC),
                                 collection = c("GO", "KEGG"),
                                 array.type = c("450K", "EPIC"),
                                 sig.genes = FALSE)
gometh_res_450k_ADC_neg = gometh(sig.cpg = rownames(sig_methyl_450k_ADC_neg),
                             all.cpg = rownames(meta_methyl_450k_ADC),
                             collection = c("GO", "KEGG"),
                             array.type = c("450K", "EPIC"),
                             sig.genes = FALSE)
##### GREY STUDY #####

methyl_table_UKBBN_greyC = limma_DA_methyl(beta_values = betas_UKBBN, pheno = pheno_UKBBN, subtype1 = "grey90", subtype2 = "Control", cohort = "UKBBN")
methyl_table_PITT_greyC = limma_DA_methyl(beta_values = betas_PITT, pheno = pheno_PITT, subtype1 = "grey90", subtype2 = "Control", cohort = "PITT")
methyl_table_ROSMAP_greyC = limma_DA_methyl(beta_values = betas_ROSMAP, pheno = pheno_ROSMAP, subtype1 = "grey90", subtype2 = "Control", cohort = "ROSMAP")

cohend_UKBBN_greyC = pbapply::pbapply(betas_UKBBN, 1, function(x){cohenD(cpg = x, pheno = pheno_UKBBN, groups = c("grey90", "Control"), AD = FALSE)})
cohend_PITT_greyC = pbapply::pbapply(betas_PITT, 1, function(x){cohenD(cpg = x, pheno = pheno_PITT, groups = c("grey90", "Control"), AD = FALSE)})
cohend_ROSMAP_greyC = pbapply::pbapply(betas_ROSMAP, 1, function(x){cohenD(cpg = x, pheno = pheno_ROSMAP, groups = c("grey90", "Control"), AD = FALSE)})

methyl_table_UKBBN_greyC = as.data.frame(methyl_table_UKBBN_greyC)
methyl_table_PITT_greyC = as.data.frame(methyl_table_PITT_greyC)
methyl_table_ROSMAP_greyC = as.data.frame(methyl_table_ROSMAP_greyC)

methyl_table_UKBBN_greyC = methyl_table_UKBBN_greyC[order(rownames(methyl_table_UKBBN_greyC)),]
methyl_table_UKBBN_greyC$effect_size = cohend_UKBBN_greyC
methyl_table_PITT_greyC = methyl_table_PITT_greyC[order(rownames(methyl_table_PITT_greyC)),]
methyl_table_PITT_greyC$effect_size = cohend_PITT_greyC
methyl_table_ROSMAP_greyC = methyl_table_ROSMAP_greyC[order(rownames(methyl_table_ROSMAP_greyC)),]
methyl_table_ROSMAP_greyC$effect_size = cohend_ROSMAP_greyC

methyl_table_UKBBN_greyC = inflation_control(methyl_table_UKBBN_greyC, subtype2 = "Control")
methyl_table_PITT_greyC = inflation_control(methyl_table_PITT_greyC, subtype2 = "Control")
methyl_table_ROSMAP_greyC = inflation_control(methyl_table_ROSMAP_greyC, subtype2 = "Control")


meta_methyl_450k_greyC = meta_analysis(BDR_res = methyl_table_BDR_blueC, 
                                       UKBBN_res = methyl_table_UKBBN_greyC, 
                                       PITT_res = methyl_table_PITT_greyC, 
                                       ROSMAP_res = methyl_table_ROSMAP_greyC,
                                       data_type = "Methyl", methyl_type = "450k")

cpgs_intersect = Reduce(intersect, list(rownames(methyl_table_UKBBN_blueC),
                                        rownames(methyl_table_PITT_blueC),
                                        rownames(methyl_table_ROSMAP_blueC)))

cohenD_cumul_450k_greyC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_greyC[rownames(methyl_table_UKBBN_greyC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_PITT" = methyl_table_PITT_greyC[rownames(methyl_table_PITT_greyC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_ROSMAP" = methyl_table_ROSMAP_greyC[rownames(methyl_table_ROSMAP_greyC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_450k_greyC) = cpgs_intersect
cohenD_cumul_450k_greyC$sig_expr = rowSums(abs(cohenD_cumul_450k_greyC)> 0.2) > 0
cohenD_cumul_450k_greyC$dir_expr = (cohenD_cumul_450k_greyC$cohens_D_UKBBN > 0 & cohenD_cumul_450k_greyC$cohens_D_ROSMAP > 0 & cohenD_cumul_450k_greyC$cohens_D_PITT > 0) | 
  (cohenD_cumul_450k_greyC$cohens_D_UKBBN < 0 & cohenD_cumul_450k_greyC$cohens_D_PITT < 0 & cohenD_cumul_450k_greyC$cohens_D_ROSMAP < 0) 


meta_methyl_EPIC_greyC = meta_analysis(BDR_res = methyl_table_BDR_blueC, 
                                       UKBBN_res = methyl_table_UKBBN_greyC, 
                                       PITT_res = methyl_table_PITT_greyC, 
                                       data_type = "Methyl", methyl_type = "EPIC")

cpgs_intersect = Reduce(intersect, list(rownames(methyl_table_UKBBN_greyC),
                                        rownames(methyl_table_PITT_greyC)))

cohenD_cumul_EPIC_greyC = data.frame("cohens_D_UKBBN" = methyl_table_UKBBN_greyC[rownames(methyl_table_UKBBN_greyC) %in% cpgs_intersect,]$effect_size,
                                     "cohens_D_PITT" = methyl_table_PITT_greyC[rownames(methyl_table_PITT_greyC) %in% cpgs_intersect,]$effect_size)
rownames(cohenD_cumul_EPIC_greyC) = cpgs_intersect
cohenD_cumul_EPIC_greyC$sig_expr = rowSums(abs(cohenD_cumul_EPIC_greyC)> 0.2) > 0
cohenD_cumul_EPIC_greyC$dir_expr = (cohenD_cumul_EPIC_greyC$cohens_D_UKBBN > 0 & cohenD_cumul_EPIC_greyC$cohens_D_PITT > 0) | 
  (cohenD_cumul_EPIC_greyC$cohens_D_UKBBN < 0 & cohenD_cumul_EPIC_greyC$cohens_D_PITT < 0) 


#
pvalue_threshold = 1e-3

sig_cohenD = rownames(cohenD_cumul_EPIC_greyC[cohenD_cumul_EPIC_greyC$dir_expr == TRUE,])

sig_methyl_EPIC_greyC = meta_methyl_EPIC_greyC[rownames(meta_methyl_EPIC_greyC) %in% sig_cohenD,]
sig_methyl_EPIC_greyC = sig_methyl_EPIC_greyC[sig_methyl_EPIC_greyC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_EPIC_greyC[rownames(cohenD_cumul_EPIC_greyC) %in% rownames(sig_methyl_EPIC_greyC),]
sig_methyl_EPIC_greyC = cbind(sig_methyl_EPIC_greyC, df_cohenD)
sig_methyl_EPIC_greyC = sig_methyl_EPIC_greyC[order(sig_methyl_EPIC_greyC$pval, decreasing = T),]
#


sig_cohenD = rownames(cohenD_cumul_450k_greyC[cohenD_cumul_450k_greyC$dir_expr == TRUE,])

sig_methyl_450k_greyC = meta_methyl_450k_greyC[rownames(meta_methyl_450k_greyC) %in% sig_cohenD,]
sig_methyl_450k_greyC = sig_methyl_450k_greyC[sig_methyl_450k_greyC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_450k_greyC[rownames(cohenD_cumul_450k_greyC) %in% rownames(sig_methyl_450k_greyC),]
sig_methyl_450k_greyC = cbind(sig_methyl_450k_greyC, df_cohenD)
sig_methyl_450k_greyC = sig_methyl_450k_greyC[order(sig_methyl_450k_greyC$pval, decreasing = T),]

tovenn <- list("BluevsCtrl" = rownames(sig_methyl_450k_blueC),
               "RedvsCtrl" = rownames(sig_methyl_450k_redC),
               "GreyvsCtrl" = rownames(sig_methyl_450k_greyC))
ggvenn(tovenn, columns = c("BluevsCtrl", "RedvsCtrl", "GreyvsCtrl"))

tovenn <- list("BluevsCtrl" = rownames(sig_methyl_450k_blueC),
               "RedvsCtrl" = rownames(sig_methyl_450k_redC),
               "GreyvsCtrl" = rownames(sig_methyl_450k_greyC),
               "ADvsCtrl" =  rownames(sig_methyl_450k_ADC))
ggvenn(tovenn, columns = c("BluevsCtrl", "RedvsCtrl", "GreyvsCtrl", "ADvsCtrl"))

meta_methyl_EPIC_greyC[rownames(meta_methyl_EPIC_greyC) == "cg05066959",]
meta_methyl_450k_greyC[rownames(meta_methyl_450k_greyC) == "cg05066959",]





#### Clinical checks ####
subtypes = c("blue", "red")
pheno_ROSMAP_AD = pheno_ROSMAP[pheno_ROSMAP$subtype %in% subtypes,]
pheno_ROSMAP_AD$subtypeNum = NA
pheno_ROSMAP_AD[pheno_ROSMAP_AD$subtype == "blue",]$subtypeNum = 0
pheno_ROSMAP_AD[pheno_ROSMAP_AD$subtype == "red",]$subtypeNum = 1
pheno_ROSMAP_AD$subtypeNum = as.numeric(pheno_ROSMAP_AD$subtypeNum)

pheno_UKBBN_AD = pheno_UKBBN[pheno_UKBBN$subtype %in% subtypes,]
pheno_UKBBN_AD$subtypeNum = NA
pheno_UKBBN_AD[pheno_UKBBN_AD$subtype == "blue",]$subtypeNum = 0
pheno_UKBBN_AD[pheno_UKBBN_AD$subtype == "red",]$subtypeNum = 1
pheno_UKBBN_AD$subtypeNum = as.numeric(pheno_UKBBN_AD$subtypeNum)

pheno_PITT_AD = pheno_PITT[pheno_PITT$subtype %in% subtypes,]
pheno_PITT_AD$subtypeNum = NA
pheno_PITT_AD[pheno_PITT_AD$subtype == "blue",]$subtypeNum = 0
pheno_PITT_AD[pheno_PITT_AD$subtype == "red",]$subtypeNum = 1
pheno_PITT_AD$subtypeNum = as.numeric(pheno_PITT_AD$subtypeNum)


table(pheno_ROSMAP_AD$apoe_genotype, pheno_ROSMAP_AD$subtypeNum)
chisq.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$apoe_genotype)


pheno_ROSMAP_AD[pheno_ROSMAP_AD$age_first_ad_dx == "90+",]$age_first_ad_dx = 90
pheno_ROSMAP_AD$age_first_ad_dx = as.numeric(pheno_ROSMAP_AD$age_first_ad_dx)

RAge1diagvSubtype = pheno_ROSMAP_AD[!is.na(pheno_ROSMAP_AD$age_first_ad_dx),]
RAge1diagvSubtype$age_first_ad_dx = signif(RAge1diagvSubtype$age_first_ad_dx, digits = 2)
table(RAge1diagvSubtype$age_first_ad_dx, RAge1diagvSubtype$subtypeNum)
cor.test(RAge1diagvSubtype$subtypeNum,RAge1diagvSubtype$age_first_ad_dx, method = "spearman")

pheno_ROSMAP_AD$age_death = signif(pheno_ROSMAP_AD$age_death, digits = 1)
table(pheno_ROSMAP_AD$age_death, pheno_ROSMAP_AD$subtypeNum)
cor.test(pheno_ROSMAP_AD$subtypeNum,pheno_ROSMAP_AD$age_death, method = "spearman")

pheno_UKBBN_AD$Age = signif(pheno_UKBBN_AD$Age, digits = 1)
table(pheno_UKBBN_AD$Age, pheno_UKBBN_AD$subtypeNum)
cor.test(pheno_UKBBN_AD$subtypeNum,pheno_UKBBN_AD$Age, method = "spearman")
pheno_PITT_AD$Age = signif(pheno_PITT_AD$Age, digits = 1)
table(pheno_PITT_AD$Age, pheno_PITT_AD$subtypeNum)
cor.test(pheno_PITT_AD$subtypeNum,pheno_PITT_AD$Age, method = "spearman")

pheno_ROSMAP_AD$cts_mmse30_lv = signif(pheno_ROSMAP_AD$cts_mmse30_lv, digits = 1)
table(pheno_ROSMAP_AD$cts_mmse30_lv,pheno_ROSMAP_AD$subtypeNum)
chisq.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$cts_mmse30_lv)
boxplot(pheno_ROSMAP_AD$cts_mmse30_lv~pheno_ROSMAP_AD$subtypeNum)

pheno_ROSMAP_AD$educ = signif(pheno_ROSMAP_AD$educ, digits = 1)
table(pheno_ROSMAP_AD$educ,pheno_ROSMAP_AD$subtypeNum)
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$educ, method = "spearman")

#### Pathology checks ####

table(pheno_ROSMAP_AD$CERAD, pheno_ROSMAP_AD$subtypeNum)
chisq.test(pheno_ROSMAP_AD$subtype,pheno_ROSMAP_AD$CERAD)

table(pheno_UKBBN_AD$Braak, pheno_UKBBN_AD$subtypeNum)
table(pheno_ROSMAP_AD$Braak, pheno_ROSMAP_AD$subtypeNum)
chisq.test(pheno_ROSMAP_AD$subtype,pheno_ROSMAP_AD$Braak)
chisq.test(pheno_UKBBN_AD$subtype,pheno_UKBBN_AD$Braak)

table(pheno_ROSMAP_AD$ceradsc, pheno_ROSMAP_AD$CERAD)
chisq.test(pheno_ROSMAP_AD$subtype,pheno_ROSMAP_AD$ceradsc)


pheno_ROSMAP_AD_CERADBraak = data.frame(pheno_ROSMAP_AD$CERAD, pheno_ROSMAP_AD$Braak)


#### Technical checks ####
pheno_ROSMAP_AD$pmi = signif(pheno_ROSMAP_AD$pmi, digits = 1)
table(pheno_ROSMAP_AD$pmi, pheno_ROSMAP_AD$subtypeNum)
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$pmi, method = "spearman")

table(pheno_ROSMAP_AD$subtype,pheno_ROSMAP_AD$batch)
chisq.test(pheno_ROSMAP_AD$subtype,pheno_ROSMAP_AD$batch)
table(pheno_ROSMAP_AD$Study,pheno_ROSMAP_AD$subtype)
chisq.test(pheno_ROSMAP_AD$subtype, pheno_ROSMAP_AD$Study)

chisq.test(pheno_UKBBN_AD$Brain.Bank,pheno_UKBBN_AD$subtype)
table(pheno_UKBBN_AD$Brain.Bank, pheno_UKBBN_AD$subtype)

table(pheno_ROSMAP_AD$Sample_Plate,pheno_ROSMAP_AD$subtype)
chisq.test(pheno_ROSMAP_AD$subtype, pheno_ROSMAP_AD$Sample_Plate)
chisq.test(pheno_PITT_AD$subtype, pheno_PITT_AD$Plate)
table(pheno_PITT_AD$subtype, pheno_PITT_AD$Plate)


cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$Neuronal_Inhibitory, method = "spearman")
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$Neuronal_Excitatory, method = "spearman")
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$Astrocyte, method = "spearman")
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$Oligodendrocyte, method = "spearman")
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$Microglia, method = "spearman")

cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$Neuronal_Inhibitory, method = "spearman")
cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$Neuronal_Excitatory, method = "spearman")
cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$Astrocyte, method = "spearman")
cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$Oligodendrocyte, method = "spearman")
cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$Microglia, method = "spearman")

cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$Neuronal_Inhibitory, method = "spearman")
cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$Neuronal_Excitatory, method = "spearman")
cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$Astrocyte, method = "spearman")
cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$Oligodendrocyte, method = "spearman")
cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$Microglia, method = "spearman")

#### PRS Check ####
ROSMAP_PRS = read.csv("SUBTYPING-PAPER-TABLEFIG/ROSMAPSummary-PRS.csv")
UKBBN_PRS = read.csv("SUBTYPING-PAPER-TABLEFIG/UKBBNSummary-PRS.csv")
PITT_PRS = read.csv("SUBTYPING-PAPER-TABLEFIG/NIHSummary-PRS.csv")

ROSMAP_PRS = ROSMAP_PRS[ROSMAP_PRS$projid %in% pheno_ROSMAP_AD$sample,]
ROSMAP_PRS = ROSMAP_PRS[order(ROSMAP_PRS$projid),]
UKBBN_PRS = UKBBN_PRS[UKBBN_PRS$BBNId %in% pheno_UKBBN_AD$BBNId,]
UKBBN_PRS = UKBBN_PRS[order(UKBBN_PRS$BBNId),]
PITT_PRS = PITT_PRS[PITT_PRS$Individual_ID %in% pheno_PITT_AD$Individual_ID,]
PITT_PRS = PITT_PRS[order(PITT_PRS$Individual_ID),]

pheno_ROSMAP_AD = pheno_ROSMAP_AD[order(pheno_ROSMAP_AD$projid),]
pheno_UKBBN_AD = pheno_UKBBN_AD[order(pheno_UKBBN_AD$BBNId),]
pheno_UKBBN_AD$BBNId = as.character(pheno_UKBBN_AD$BBNId)
pheno_PITT_AD = pheno_PITT_AD[order(pheno_PITT_AD$Individual_ID),]

pheno_ROSMAP_AD = cbind(pheno_ROSMAP_AD, ROSMAP_PRS)
pheno_UKBBN_AD = cbind(pheno_UKBBN_AD, UKBBN_PRS)
pheno_PITT_AD = cbind(pheno_PITT_AD, PITT_PRS)

cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$PRS, method = "spearman")
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$PC1, method = "spearman")
cor.test(pheno_ROSMAP_AD$subtypeNum, pheno_ROSMAP_AD$PC2, method = "spearman")

cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$PRS, method = "spearman")
cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$PC1, method = "spearman")
cor.test(pheno_UKBBN_AD$subtypeNum, pheno_UKBBN_AD$PC2, method = "spearman")

cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$PRS, method = "spearman")
cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$PC1, method = "spearman")
cor.test(pheno_PITT_AD$subtypeNum, pheno_PITT_AD$PC2, method = "spearman")



#### RF genotype ####
load("phenos_3cohorts.Rdata")
# red_SNPset = read.table("SNP_data/table_mqtls_450kredC.clumped", header = T)
# blue_SNPset = read.table("SNP_data/table_mqtls_450kblueC.clumped", header = T)

ROSMAP_SNPset_red = read.table("SNP_data/ROSMAP_mQTLset_red.raw", sep = " ", header = T)
ROSMAP_SNPset_blue = read.table("SNP_data/ROSMAP_mQTLset_blue.raw", sep = " ", header = T)
ROSMAP_SNPset_red = ROSMAP_SNPset_red[ROSMAP_SNPset_red$FID %in% pheno_ROSMAP$sample,]
ROSMAP_SNPset_blue = ROSMAP_SNPset_blue[ROSMAP_SNPset_blue$FID %in% pheno_ROSMAP$sample,]
pheno_SNP_ROSMAP = pheno_ROSMAP[pheno_ROSMAP$sample %in% ROSMAP_SNPset_red$FID,]

UKBBN_SNPset_red = read.table("SNP_data/UKBBN_mQTLset_red.raw", sep = " ", header = T)
UKBBN_SNPset_blue = read.table("SNP_data/UKBBN_mQTLset_blue.raw", sep = " ", header = T)
UKBBN_SNPset_red = UKBBN_SNPset_red[UKBBN_SNPset_red$BBNId %in% pheno_UKBBN$BBNId,]
UKBBN_SNPset_blue = UKBBN_SNPset_blue[UKBBN_SNPset_blue$BBNId %in% pheno_UKBBN$BBNId,]
pheno_SNP_UKBBN = pheno_UKBBN[pheno_UKBBN$BBNId %in% UKBBN_SNPset_red$BBNId,]

PITT_SNPset_red = read.table("SNP_data/NIH_mQTLset_red.raw", sep = "\t", header = T)
PITT_SNPset_blue = read.table("SNP_data/NIH_mQTLset_blue.raw", sep = "\t", header = T)
PITT_SNPset_red = PITT_SNPset_red[PITT_SNPset_red$FID %in% pheno_PITT$Individual_ID,]
PITT_SNPset_blue = PITT_SNPset_blue[PITT_SNPset_blue$FID %in% pheno_PITT$Individual_ID,]
pheno_SNP_PITT = pheno_PITT[pheno_PITT$Individual_ID %in% PITT_SNPset_red$FID,]

ROSMAP_SNPset_red = ROSMAP_SNPset_red[order(ROSMAP_SNPset_red$FID),]
ROSMAP_SNPset_blue = ROSMAP_SNPset_blue[order(ROSMAP_SNPset_blue$FID),]
UKBBN_SNPset_red = UKBBN_SNPset_red[order(UKBBN_SNPset_red$BBNId),]
UKBBN_SNPset_blue = UKBBN_SNPset_blue[order(UKBBN_SNPset_blue$BBNId),]
PITT_SNPset_red = PITT_SNPset_red[order(PITT_SNPset_red$FID),]
PITT_SNPset_blue = PITT_SNPset_blue[order(PITT_SNPset_blue$FID),]

pheno_SNP_ROSMAP = pheno_SNP_ROSMAP[order(pheno_SNP_ROSMAP$sample),]
pheno_SNP_UKBBN = pheno_SNP_UKBBN[order(pheno_SNP_UKBBN$BBNId),]
pheno_SNP_PITT = pheno_SNP_PITT[order(pheno_SNP_PITT$Individual_ID),]
identical(pheno_SNP_ROSMAP$sample, ROSMAP_SNPset_red$FID)

rownames(ROSMAP_SNPset_red) = ROSMAP_SNPset_red$FID
rownames(ROSMAP_SNPset_blue) = ROSMAP_SNPset_blue$FID
rownames(UKBBN_SNPset_red) = UKBBN_SNPset_red$BBNId
rownames(UKBBN_SNPset_blue) = UKBBN_SNPset_blue$BBNId
rownames(PITT_SNPset_red) = PITT_SNPset_red$FID
rownames(PITT_SNPset_blue) = PITT_SNPset_blue$FID

ROSMAP_SNPset_red = ROSMAP_SNPset_red[,-c(1:6)]
ROSMAP_SNPset_blue = ROSMAP_SNPset_blue[,-c(1:6)]
UKBBN_SNPset_red = UKBBN_SNPset_red[,-c(1:6)]
UKBBN_SNPset_blue = UKBBN_SNPset_blue[,-c(1:6)]
PITT_SNPset_red = PITT_SNPset_red[,-c(1:6)]
PITT_SNPset_blue = PITT_SNPset_blue[,-c(1:6)]

table(pheno_SNP_ROSMAP$subtype)
table(pheno_SNP_UKBBN$subtype)
table(pheno_SNP_PITT$subtype)

col_ROSMAP_red_orig = colnames(ROSMAP_SNPset_red)
col_ROSMAP_red_rep = gsub("X","", col_ROSMAP_red_orig)
col_ROSMAP_red_rep = gsub(".[A-Z].[A-Z]_[A-Z]", "", col_ROSMAP_red_rep)
grep("[A-Z]", col_ROSMAP_red_rep)
colnames(ROSMAP_SNPset_red) = col_ROSMAP_red_rep
# red_SNPset$SNPname = paste(red_SNPset$CHR, red_SNPset$BP, sep = ".")
# intersect(red_SNPset$SNPname, col_ROSMAP_red_rep)
# ROSMAP_SNPset_red = ROSMAP_SNPset_red[,colnames(ROSMAP_SNPset_red) %in% intersect(red_SNPset$SNPname, col_ROSMAP_red_rep)]


col_ROSMAP_blue_orig = colnames(ROSMAP_SNPset_blue)
col_ROSMAP_blue_rep = gsub("X","", col_ROSMAP_blue_orig)
col_ROSMAP_blue_rep = gsub(".[A-Z].[A-Z]_[A-Z]", "", col_ROSMAP_blue_rep)
grep("[A-Z]", col_ROSMAP_blue_rep)
colnames(ROSMAP_SNPset_blue) = col_ROSMAP_blue_rep
# blue_SNPset$SNPname = paste(blue_SNPset$CHR, blue_SNPset$BP, sep = ".")
# intersect(blue_SNPset$SNPname, col_ROSMAP_blue_rep)
# ROSMAP_SNPset_blue = ROSMAP_SNPset_blue[,colnames(ROSMAP_SNPset_blue) %in% intersect(blue_SNPset$SNPname, col_ROSMAP_blue_rep)]


col_UKBBN_red_orig = colnames(UKBBN_SNPset_red)
col_UKBBN_red_rep = gsub("X","", col_UKBBN_red_orig)
col_UKBBN_red_rep = gsub(".[A-Z].[A-Z]_[A-Z]", "", col_UKBBN_red_rep)
grep("[A-Z]", col_UKBBN_red_rep)
colnames(UKBBN_SNPset_red) = col_UKBBN_red_rep
# UKBBN_SNPset_red = UKBBN_SNPset_red[,colnames(UKBBN_SNPset_red) %in% intersect(red_SNPset$SNPname, col_UKBBN_red_rep)]

col_UKBBN_blue_orig = colnames(UKBBN_SNPset_blue)
col_UKBBN_blue_rep = gsub("X","", col_UKBBN_blue_orig)
col_UKBBN_blue_rep = gsub(".[A-Z].[A-Z]_[A-Z]", "", col_UKBBN_blue_rep)
grep("[A-Z]", col_UKBBN_blue_rep)
colnames(UKBBN_SNPset_blue) = col_UKBBN_blue_rep
# UKBBN_SNPset_blue = UKBBN_SNPset_blue[,colnames(UKBBN_SNPset_blue) %in% intersect(blue_SNPset$SNPname, col_UKBBN_blue_rep)]

col_PITT_red_orig = colnames(PITT_SNPset_red)
col_PITT_red_rep = gsub("X","", col_PITT_red_orig)
col_PITT_red_rep = gsub(".[A-Z].[A-Z]_[A-Z]", "", col_PITT_red_rep)
grep("[A-Z]", col_PITT_red_rep)
colnames(PITT_SNPset_red) = col_PITT_red_rep
# PITT_SNPset_red = PITT_SNPset_red[,colnames(PITT_SNPset_red) %in% intersect(red_SNPset$SNPname, col_PITT_red_rep)]

col_PITT_blue_orig = colnames(PITT_SNPset_blue)
col_PITT_blue_rep = gsub("X","", col_PITT_blue_orig)
col_PITT_blue_rep = gsub(".[A-Z].[A-Z]_[A-Z]", "", col_PITT_blue_rep)
grep("[A-Z]", col_PITT_blue_rep)
colnames(PITT_SNPset_blue) = col_PITT_blue_rep
# PITT_SNPset_blue = PITT_SNPset_blue[,colnames(PITT_SNPset_blue) %in% intersect(blue_SNPset$SNPname, col_PITT_blue_rep)]

#### PC data ####
ROSMAP_PCA = read.csv("SNP_data/ROSMAPSummary.csv")
UKBBN_PCA = read.csv("SNP_data/UKBBNSummary.csv")
PITT_PCA = read.csv("SNP_data/NIHSummary.csv")

pheno_SNP_ROSMAP = pheno_SNP_ROSMAP[pheno_SNP_ROSMAP$sample %in% ROSMAP_PCA$projid,]
ROSMAP_PCA = ROSMAP_PCA[ROSMAP_PCA$projid %in% pheno_SNP_ROSMAP$sample,]
ROSMAP_PCA = ROSMAP_PCA[order(ROSMAP_PCA$projid),]
pheno_SNP_ROSMAP = pheno_SNP_ROSMAP[order(pheno_SNP_ROSMAP$sample),]
pheno_SNP_ROSMAP = cbind(pheno_SNP_ROSMAP, ROSMAP_PCA[,c(3:ncol(ROSMAP_PCA))])

UKBBN_PCA = UKBBN_PCA[UKBBN_PCA$BBNId %in% pheno_SNP_UKBBN$BBNId,]
UKBBN_PCA = UKBBN_PCA[order(UKBBN_PCA$BBNId),]
pheno_SNP_UKBBN = pheno_SNP_UKBBN[order(pheno_SNP_UKBBN$BBNId),]
pheno_SNP_UKBBN$BBNId = as.character(pheno_SNP_UKBBN$BBNId)
identical(pheno_SNP_UKBBN$BBNId, UKBBN_PCA$BBNId)
pheno_SNP_UKBBN = cbind(pheno_SNP_UKBBN, UKBBN_PCA[,c(4:ncol(UKBBN_PCA))])


PITT_PCA = PITT_PCA[PITT_PCA$Individual_ID %in% pheno_SNP_PITT$Individual_ID,]
PITT_PCA = PITT_PCA[order(PITT_PCA$Individual_ID),]
identical(pheno_SNP_PITT$Individual_ID, PITT_PCA$Individual_ID)
pheno_SNP_PITT = cbind(pheno_SNP_PITT, PITT_PCA[,c(3:ncol(PITT_PCA))])



#### RED SNPs ####
red_coloc_snps = read.table( file = "./SNP_data/RED_snps.txt")[,1]
red_coloc_snps = gsub(":",".", red_coloc_snps)
blue_coloc_snps = read.table( file = "./SNP_data/blue_snps.txt")[,1]
blue_coloc_snps = gsub(":",".", blue_coloc_snps)
#### vs ctrl
pheno_SNP_ROSMAP_red = pheno_SNP_ROSMAP[pheno_SNP_ROSMAP$subtype %in% c("red", "Control"),]

ROSMAP_SNPset_redC = ROSMAP_SNPset_red[rownames(ROSMAP_SNPset_red) %in% pheno_SNP_ROSMAP_red$sample,]
ROSMAP_SNPset_redC = ROSMAP_SNPset_redC[,colnames(ROSMAP_SNPset_redC) %in% red_coloc_snps]
ROSMAP_SNPset_redC = as.data.frame(t(ROSMAP_SNPset_redC))

snp_table_ROSMAP_redC = limma_DA_SNP(snp_genotype = ROSMAP_SNPset_redC, pheno = pheno_SNP_ROSMAP_red, subtype1 = "red", subtype2 = "Control", cohort = "ROSMAP")
snp_table_ROSMAP_redC = as.data.frame(snp_table_ROSMAP_redC)

pheno_SNP_ROSMAP_blue = pheno_SNP_ROSMAP[pheno_SNP_ROSMAP$subtype %in% c("blue", "Control"),]

ROSMAP_SNPset_blueC = ROSMAP_SNPset_blue[rownames(ROSMAP_SNPset_blue) %in% pheno_SNP_ROSMAP_blue$sample,]
ROSMAP_SNPset_blueC = ROSMAP_SNPset_blueC[,colnames(ROSMAP_SNPset_blueC) %in% blue_coloc_snps]
ROSMAP_SNPset_blueC = as.data.frame(t(ROSMAP_SNPset_blueC))

snp_table_ROSMAP_blueC = limma_DA_SNP(snp_genotype = ROSMAP_SNPset_blueC, pheno = pheno_SNP_ROSMAP_blue, subtype1 = "blue", subtype2 = "Control", cohort = "ROSMAP")
snp_table_ROSMAP_blueC = as.data.frame(snp_table_ROSMAP_blueC)
#
pheno_SNP_UKBBN_red = pheno_SNP_UKBBN[pheno_SNP_UKBBN$subtype %in% c("red", "Control"),]

UKBBN_SNPset_redC = UKBBN_SNPset_red[rownames(UKBBN_SNPset_red) %in% pheno_SNP_UKBBN_red$BBNId,]
UKBBN_SNPset_redC = UKBBN_SNPset_redC[,colnames(UKBBN_SNPset_redC) %in% red_coloc_snps]
UKBBN_SNPset_redC = as.data.frame(t(UKBBN_SNPset_redC))

snp_table_UKBBN_redC = limma_DA_SNP(snp_genotype = UKBBN_SNPset_redC, pheno = pheno_SNP_UKBBN_red, subtype1 = "red", subtype2 = "Control", cohort = "UKBBN")
snp_table_UKBBN_redC = as.data.frame(snp_table_UKBBN_redC)

pheno_SNP_UKBBN_blue = pheno_SNP_UKBBN[pheno_SNP_UKBBN$subtype %in% c("blue", "Control"),]

UKBBN_SNPset_blueC = UKBBN_SNPset_blue[rownames(UKBBN_SNPset_blue) %in% pheno_SNP_UKBBN_blue$BBNId,]
UKBBN_SNPset_blueC = UKBBN_SNPset_blueC[,colnames(UKBBN_SNPset_blueC) %in% blue_coloc_snps]
UKBBN_SNPset_blueC = as.data.frame(t(UKBBN_SNPset_blueC))

snp_table_UKBBN_blueC = limma_DA_SNP(snp_genotype = UKBBN_SNPset_blueC, pheno = pheno_SNP_UKBBN_blue, subtype1 = "blue", subtype2 = "Control", cohort = "UKBBN")
snp_table_UKBBN_blueC = as.data.frame(snp_table_UKBBN_blueC)
#
pheno_SNP_PITT_red = pheno_SNP_PITT[pheno_SNP_PITT$subtype %in% c("red", "Control"),]

PITT_SNPset_redC = PITT_SNPset_red[rownames(PITT_SNPset_red) %in% pheno_SNP_PITT_red$Individual_ID,]
PITT_SNPset_redC = PITT_SNPset_redC[,colnames(PITT_SNPset_redC) %in% red_coloc_snps]
PITT_SNPset_redC = as.data.frame(t(PITT_SNPset_redC))

snp_table_PITT_redC = limma_DA_SNP(snp_genotype = PITT_SNPset_redC, pheno = pheno_SNP_PITT_red, subtype1 = "red", subtype2 = "Control", cohort = "PITT")
snp_table_PITT_redC = as.data.frame(snp_table_PITT_redC)

pheno_SNP_PITT_blue = pheno_SNP_PITT[pheno_SNP_PITT$subtype %in% c("blue", "Control"),]

PITT_SNPset_blueC = PITT_SNPset_blue[rownames(PITT_SNPset_blue) %in% pheno_SNP_PITT_blue$Individual_ID,]
PITT_SNPset_blueC = PITT_SNPset_blueC[,colnames(PITT_SNPset_blueC) %in% blue_coloc_snps]
PITT_SNPset_blueC = as.data.frame(t(PITT_SNPset_blueC))

snp_table_PITT_blueC = limma_DA_SNP(snp_genotype = PITT_SNPset_blueC, pheno = pheno_SNP_PITT_blue, subtype1 = "blue", subtype2 = "Control", cohort = "PITT")
snp_table_PITT_blueC = as.data.frame(snp_table_PITT_blueC)


snp_table_ROSMAP_redC$cohort = "ROSMAP"
snp_table_ROSMAP_blueC$cohort = "ROSMAP"
snp_table_UKBBN_redC$cohort = "UKBBN"
snp_table_UKBBN_blueC$cohort = "UKBBN"
snp_table_PITT_redC$cohort = "PITT"
snp_table_PITT_blueC$cohort = "PITT"

snp_table_ROSMAP_redC = snp_table_ROSMAP_redC[snp_table_ROSMAP_redC$P.t1 < 0.05,]
snp_table_ROSMAP_blueC = snp_table_ROSMAP_blueC[snp_table_ROSMAP_blueC$P.t1 < 0.05,]
snp_table_UKBBN_redC = snp_table_UKBBN_redC[snp_table_UKBBN_redC$P.t1 < 0.05,]
snp_table_UKBBN_blueC = snp_table_UKBBN_blueC[snp_table_UKBBN_blueC$P.t1 < 0.05,]
snp_table_PITT_redC = snp_table_PITT_redC[snp_table_PITT_redC$P.t1 < 0.05,]
snp_table_PITT_blueC = snp_table_PITT_blueC[snp_table_PITT_blueC$P.t1 < 0.05,]

snp_table_all_redC = rbind(snp_table_ROSMAP_redC,snp_table_UKBBN_redC,snp_table_PITT_redC)
snp_table_all_redC = snp_table_all_redC[order(rownames(snp_table_all_redC)),]
snp_table_all_blueC = rbind(snp_table_ROSMAP_blueC,snp_table_UKBBN_blueC,snp_table_PITT_blueC)
snp_table_all_blueC = snp_table_all_blueC[order(rownames(snp_table_all_blueC)),]


#### vs AD

pheno_SNP_ROSMAP_redAD = pheno_SNP_ROSMAP[pheno_SNP_ROSMAP$subtype != "Control",]
pheno_SNP_ROSMAP_redAD$subtype_reduc = pheno_SNP_ROSMAP_redAD$subtype
pheno_SNP_ROSMAP_redAD[pheno_SNP_ROSMAP_redAD$subtype_reduc != "red",]$subtype_reduc = "other"
ROSMAP_SNPset_redAD = ROSMAP_SNPset_red[rownames(ROSMAP_SNPset_red) %in% pheno_SNP_ROSMAP_redAD$sample,]
ROSMAP_SNPset_redAD = ROSMAP_SNPset_redAD[,colnames(ROSMAP_SNPset_redAD) %in% red_coloc_snps]
ROSMAP_SNPset_redAD = as.data.frame(t(ROSMAP_SNPset_redAD))

snp_table_ROSMAP_redAD = limma_DA_SNP(snp_genotype = ROSMAP_SNPset_redAD, pheno = pheno_SNP_ROSMAP_redAD, subtype1 = "red", subtype2 = "other", cohort = "ROSMAP", AD = FALSE)


####

ROSMAP_SNPset_red_BIN = ROSMAP_SNPset_red[,colnames(ROSMAP_SNPset_red) %in% red_coloc_snps[1:63]]
ROSMAP_SNPset_red_BIN = ROSMAP_SNPset_red_BIN[order(rownames(ROSMAP_SNPset_red_BIN)),]
pheno_SNP_ROSMAP_red$sample = as.character(pheno_SNP_ROSMAP_red$sample)
pheno_SNP_ROSMAP_red = pheno_SNP_ROSMAP_red[order(pheno_SNP_ROSMAP_red$sample),]
ROSMAP_SNPset_red_BIN = ROSMAP_SNPset_red_BIN[rownames(ROSMAP_SNPset_red_BIN) %in% pheno_SNP_ROSMAP_red$sample,]
identical(rownames(ROSMAP_SNPset_red_BIN), pheno_SNP_ROSMAP_red$sample)
ROSMAP_SNPset_red_BIN$diag = pheno_SNP_ROSMAP_red$subtype

UKBBN_SNPset_red_BIN = UKBBN_SNPset_red[,colnames(UKBBN_SNPset_red) %in% red_coloc_snps[1:63]]
UKBBN_SNPset_red_BIN = UKBBN_SNPset_red_BIN[order(rownames(UKBBN_SNPset_red_BIN)),]
pheno_SNP_UKBBN_red$BBNId = as.character(pheno_SNP_UKBBN_red$BBNId)
pheno_SNP_UKBBN_red = pheno_SNP_UKBBN_red[order(pheno_SNP_UKBBN_red$BBNId),]
UKBBN_SNPset_red_BIN = UKBBN_SNPset_red_BIN[rownames(UKBBN_SNPset_red_BIN) %in% pheno_SNP_UKBBN_red$BBNId,]
identical(rownames(UKBBN_SNPset_red_BIN), pheno_SNP_UKBBN_red$BBNId)
UKBBN_SNPset_red_BIN$diag = pheno_SNP_UKBBN_red$subtype

PITT_SNPset_red_BIN = PITT_SNPset_red[,colnames(PITT_SNPset_red) %in% red_coloc_snps[1:63]]
PITT_SNPset_red_BIN = PITT_SNPset_red_BIN[order(rownames(PITT_SNPset_red_BIN)),]
pheno_SNP_PITT_red$Individual_ID = as.character(pheno_SNP_PITT_red$Individual_ID)
pheno_SNP_PITT_red = pheno_SNP_PITT_red[order(pheno_SNP_PITT_red$Individual_ID),]
PITT_SNPset_red_BIN = PITT_SNPset_red_BIN[rownames(PITT_SNPset_red_BIN) %in% pheno_SNP_PITT_red$Individual_ID,]
identical(rownames(PITT_SNPset_red_BIN), pheno_SNP_PITT_red$Individual_ID)
PITT_SNPset_red_BIN$diag = pheno_SNP_PITT_red$subtype

all_SNPset_red_BIN = rbind(ROSMAP_SNPset_red_BIN, UKBBN_SNPset_red_BIN, PITT_SNPset_red_BIN)

for(i in seq(1:37)){
  print(chisq.test(table(all_SNPset_red_BIN[,i], all_SNPset_red_BIN$diag)))
}
chisq.test(table(all_SNPset_red_BIN$`2.127772136`, all_SNPset_red_BIN$diag))
