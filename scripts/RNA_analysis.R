setwd("D:/valentin/main/")

load("RNAseq/current.Rdata")

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

#####
plot.new()

boxplot(betas_BDR, outline = FALSE)
boxplot(betas_ROSMAP, outline = FALSE)
boxplot(betas_PITT, outline = FALSE)
boxplot(betas_UKBBN, outline = FALSE)

plot(density(as.matrix(betas_BDR)))
lines(density(as.matrix(betas_UKBBN)), col = "red")
lines(density(as.matrix(betas_PITT)), col = "blue")

#####
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

#####

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

#####

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

#####
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

#####
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

#####

meta_methyl_450k_blueC = meta_analysis(BDR_res = methyl_table_BDR_blueC, 
                                       UKBBN_res = methyl_table_UKBBN_blueC, 
                                       PITT_res = methyl_table_PITT_blueC, 
                                       ROSMAP_res = methyl_table_ROSMAP_blueC,
                                       data_type = "Methyl", methyl_type = "450k")
meta_methyl_450k_redC = meta_analysis(BDR_res = methyl_table_BDR_redC, 
                                       UKBBN_res = methyl_table_UKBBN_redC, 
                                       PITT_res = methyl_table_PITT_redC, 
                                       ROSMAP_res = methyl_table_ROSMAP_redC,
                                       data_type = "Methyl", methyl_type = "450k")
meta_methyl_450k_redblue = meta_analysis(BDR_res = methyl_table_BDR_redblue, 
                                       UKBBN_res = methyl_table_UKBBN_redblue, 
                                       PITT_res = methyl_table_PITT_redblue, 
                                       ROSMAP_res = methyl_table_ROSMAP_redblue,
                                       data_type = "Methyl", methyl_type = "450k")
meta_methyl_450k_ADC = meta_analysis(BDR_res = methyl_table_BDR_ADC, 
                                       UKBBN_res = methyl_table_UKBBN_ADC, 
                                       PITT_res = methyl_table_PITT_ADC, 
                                       ROSMAP_res = methyl_table_ROSMAP_ADC,
                                       data_type = "Methyl", methyl_type = "450k")

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


#####

meta_methyl_EPIC_blueC = meta_analysis(BDR_res = methyl_table_BDR_blueC, 
                                       UKBBN_res = methyl_table_UKBBN_blueC, 
                                       PITT_res = methyl_table_PITT_blueC, 
                                       data_type = "Methyl", methyl_type = "EPIC")
meta_methyl_EPIC_redC = meta_analysis(BDR_res = methyl_table_BDR_redC, 
                                      UKBBN_res = methyl_table_UKBBN_redC, 
                                      PITT_res = methyl_table_PITT_redC,
                                      data_type = "Methyl", methyl_type = "EPIC")
meta_methyl_EPIC_redblue = meta_analysis(BDR_res = methyl_table_BDR_redblue, 
                                         UKBBN_res = methyl_table_UKBBN_redblue, 
                                         PITT_res = methyl_table_PITT_redblue, 
                                         data_type = "Methyl", methyl_type = "EPIC")
meta_methyl_EPIC_ADC = meta_analysis(BDR_res = methyl_table_BDR_ADC, 
                                     UKBBN_res = methyl_table_UKBBN_ADC, 
                                     PITT_res = methyl_table_PITT_ADC, 
                                     data_type = "Methyl", methyl_type = "EPIC")

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
#
pvalue_threshold = 1e-3

sig_cohenD = rownames(cohenD_cumul_EPIC_blueC[cohenD_cumul_EPIC_blueC$dir_expr == TRUE &
                                                cohenD_cumul_EPIC_blueC$sig_expr == TRUE,])

sig_methyl_EPIC_blueC = meta_methyl_EPIC_blueC[rownames(meta_methyl_EPIC_blueC) %in% sig_cohenD,]
sig_methyl_EPIC_blueC = sig_methyl_EPIC_blueC[sig_methyl_EPIC_blueC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_EPIC_blueC[rownames(cohenD_cumul_EPIC_blueC) %in% rownames(sig_methyl_EPIC_blueC),]
sig_methyl_EPIC_blueC = cbind(sig_methyl_EPIC_blueC, df_cohenD)
sig_methyl_EPIC_blueC = sig_methyl_EPIC_blueC[order(sig_methyl_EPIC_blueC$pval, decreasing = T),]

#
sig_cohenD = rownames(cohenD_cumul_EPIC_redC[cohenD_cumul_EPIC_redC$dir_expr == TRUE &
                                               cohenD_cumul_EPIC_redC$sig_expr == TRUE,])

sig_methyl_EPIC_redC = meta_methyl_EPIC_redC[rownames(meta_methyl_EPIC_redC) %in% sig_cohenD,]
sig_methyl_EPIC_redC = sig_methyl_EPIC_redC[sig_methyl_EPIC_redC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_EPIC_redC[rownames(cohenD_cumul_EPIC_redC) %in% rownames(sig_methyl_EPIC_redC),]
sig_methyl_EPIC_redC = cbind(sig_methyl_EPIC_redC, df_cohenD)
sig_methyl_EPIC_redC = sig_methyl_EPIC_redC[order(sig_methyl_EPIC_redC$pval, decreasing = T),]

#
sig_cohenD = rownames(cohenD_cumul_EPIC_redblue[cohenD_cumul_EPIC_redblue$dir_expr == TRUE &
                                                  cohenD_cumul_EPIC_redblue$sig_expr == TRUE,])

sig_methyl_EPIC_redblue = meta_methyl_EPIC_redblue[rownames(meta_methyl_EPIC_redblue) %in% sig_cohenD,]
sig_methyl_EPIC_redblue = sig_methyl_EPIC_redblue[sig_methyl_EPIC_redblue$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_EPIC_redblue[rownames(cohenD_cumul_EPIC_redblue) %in% rownames(sig_methyl_EPIC_redblue),]
sig_methyl_EPIC_redblue = cbind(sig_methyl_EPIC_redblue, df_cohenD)
sig_methyl_EPIC_redblue = sig_methyl_EPIC_redblue[order(sig_methyl_EPIC_redblue$pval, decreasing = T),]

#
sig_cohenD = rownames(cohenD_cumul_EPIC_ADC[cohenD_cumul_EPIC_ADC$dir_expr == TRUE &
                                              cohenD_cumul_EPIC_ADC$sig_expr == TRUE,])

sig_methyl_EPIC_ADC = meta_methyl_EPIC_ADC[rownames(meta_methyl_EPIC_ADC) %in% sig_cohenD,]
sig_methyl_EPIC_ADC = sig_methyl_EPIC_ADC[sig_methyl_EPIC_ADC$pval< pvalue_threshold,]

df_cohenD = cohenD_cumul_EPIC_ADC[rownames(cohenD_cumul_EPIC_ADC) %in% rownames(sig_methyl_EPIC_ADC),]
sig_methyl_EPIC_ADC = cbind(sig_methyl_EPIC_ADC, df_cohenD)
sig_methyl_EPIC_ADC = sig_methyl_EPIC_ADC[order(sig_methyl_EPIC_ADC$pval, decreasing = T),]

tovenn <- list("BluevCtrl" = rownames(sig_methyl_EPIC_blueC),
               "RedvCtrl" = rownames(sig_methyl_EPIC_redC),
               "RedvBlue" = rownames(sig_methyl_EPIC_redblue))
ggvenn(tovenn, columns = c("BluevCtrl", "RedvCtrl", "RedvBlue")) 

tovenn <- list("BluevCtrl" = rownames(sig_methyl_EPIC_blueC),
               "RedvCtrl" = rownames(sig_methyl_EPIC_redC),
               "ADvCtrl" = rownames(sig_methyl_EPIC_ADC))
ggvenn(tovenn, columns = c("BluevCtrl", "RedvCtrl", "ADvCtrl")) 

sig_methyl_EPIC_blueC_neg = sig_methyl_EPIC_blueC[sig_methyl_EPIC_blueC$cohens_D_UKBBN < 0,]
sig_methyl_EPIC_blueC_pos = sig_methyl_EPIC_blueC[sig_methyl_EPIC_blueC$cohens_D_UKBBN >= 0,]

sig_methyl_EPIC_redC_neg = sig_methyl_EPIC_redC[sig_methyl_EPIC_redC$cohens_D_UKBBN < 0,]
sig_methyl_EPIC_redC_pos = sig_methyl_EPIC_redC[sig_methyl_EPIC_redC$cohens_D_UKBBN >= 0,]

sig_methyl_EPIC_ADC_neg = sig_methyl_EPIC_ADC[sig_methyl_EPIC_ADC$cohens_D_UKBBN < 0,]
sig_methyl_EPIC_ADC_pos = sig_methyl_EPIC_ADC[sig_methyl_EPIC_ADC$cohens_D_UKBBN >= 0,]

sig_methyl_EPIC_redblue_neg = sig_methyl_EPIC_redblue[sig_methyl_EPIC_redblue$cohens_D_UKBBN < 0,]
sig_methyl_EPIC_redblue_pos = sig_methyl_EPIC_redblue[sig_methyl_EPIC_redblue$cohens_D_UKBBN >= 0,]

tovenn <- list("BluevCtrl" = rownames(sig_methyl_EPIC_blueC_neg),
               "RedvCtrl" = rownames(sig_methyl_EPIC_redC_neg),
               "ADvCtrl" = rownames(sig_methyl_EPIC_ADC_neg))
ggvenn(tovenn, columns = c("BluevCtrl", "RedvCtrl", "ADvCtrl"))


tovenn <- list("BluevCtrl" = rownames(sig_methyl_EPIC_blueC_pos),
               "RedvCtrl" = rownames(sig_methyl_EPIC_redC_pos),
               "ADvCtrl" = rownames(sig_methyl_EPIC_ADC_pos))
ggvenn(tovenn, columns = c("BluevCtrl", "RedvCtrl", "ADvCtrl"))

#####

#
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

inter_red_blue = intersect(rownames(sig_methyl_450k_blueC), rownames(sig_methyl_450k_redC))
sig_methyl_450k_blueC = sig_methyl_450k_blueC[!rownames(sig_methyl_450k_blueC) %in% inter_red_blue,]
sig_methyl_450k_redC = sig_methyl_450k_redC[!rownames(sig_methyl_450k_redC) %in% inter_red_blue,]

inter_red_blue = intersect(rownames(sig_methyl_EPIC_blueC), rownames(sig_methyl_EPIC_redC))
sig_methyl_EPIC_blueC = sig_methyl_EPIC_blueC[!rownames(sig_methyl_EPIC_blueC) %in% inter_red_blue,]
sig_methyl_EPIC_redC = sig_methyl_EPIC_redC[!rownames(sig_methyl_EPIC_redC) %in% inter_red_blue,]


CpGs_previous = read.table("AD_BH_CpGs.txt")

intersect(CpGs_previous$V1, rownames(sig_methyl_450k_ADC))
intersect(CpGs_previous$V1, rownames(sig_methyl_450k_blueC))
intersect(CpGs_previous$V1, rownames(sig_methyl_450k_redblue))
intersect(CpGs_previous$V1, rownames(sig_methyl_450k_redC))
intersect(CpGs_previous$V1, rownames(sig_methyl_EPIC_ADC))
intersect(CpGs_previous$V1, rownames(sig_methyl_EPIC_blueC))
intersect(CpGs_previous$V1, rownames(sig_methyl_EPIC_redblue))
intersect(CpGs_previous$V1, rownames(sig_methyl_EPIC_redC))

Reduce(intersect, list(CpGs_previous$V1, rownames(sig_methyl_450k_blueC),
                       rownames(sig_methyl_450k_redC)))
Reduce(intersect, list(CpGs_previous$V1, rownames(sig_methyl_450k_redC),
                       rownames(sig_methyl_450k_ADC)))

gctrl = gap::gcontrol2(meta_methyl_EPIC_ADC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_450k_ADC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_EPIC_blueC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_450k_blueC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_EPIC_redC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_450k_redC$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_EPIC_redblue$pval)
print(gctrl$lambda)
gctrl = gap::gcontrol2(meta_methyl_450k_redblue$pval)
print(gctrl$lambda)

meta_methyl_EPIC_ADC[rownames(meta_methyl_EPIC_ADC) == "cg05066959",]
meta_methyl_450k_ADC[rownames(meta_methyl_450k_ADC) == "cg05066959",]
meta_methyl_EPIC_redC[rownames(meta_methyl_EPIC_redC) == "cg05066959",]
meta_methyl_450k_redC[rownames(meta_methyl_450k_redC) == "cg05066959",]
meta_methyl_EPIC_blueC[rownames(meta_methyl_EPIC_blueC) == "cg05066959",]
meta_methyl_450k_blueC[rownames(meta_methyl_450k_blueC) == "cg05066959",]

methyl_table_UKBBN_ADC[rownames(methyl_table_UKBBN_ADC) == "cg05066959",]
methyl_table_PITT_ADC[rownames(methyl_table_PITT_ADC) == "cg05066959",]
methyl_table_ROSMAP_ADC[rownames(methyl_table_ROSMAP_ADC) == "cg05066959",]

##### gometh EPIC #####

gometh_res_EPIC_blueC = gometh(sig.cpg = rownames(sig_methyl_EPIC_blueC),
                                   all.cpg = rownames(meta_methyl_EPIC_ADC),
                                   collection = c("GO", "KEGG"),
                                   array.type = c("450K", "EPIC"),
                                   sig.genes = FALSE)

gometh_res_EPIC_redC = gometh(sig.cpg = rownames(sig_methyl_EPIC_redC),
                                  all.cpg = rownames(meta_methyl_EPIC_ADC),
                                  collection = c("GO", "KEGG"),
                                  array.type = c("450K", "EPIC"),
                                  sig.genes = FALSE)
gometh_res_EPIC_redblue = gometh(sig.cpg = rownames(sig_methyl_EPIC_redblue),
                                     all.cpg = rownames(meta_methyl_EPIC_ADC),
                                     collection = c("GO", "KEGG"),
                                     array.type = c("450K", "EPIC"),
                                     sig.genes = FALSE)
gometh_res_EPIC_ADC = gometh(sig.cpg = rownames(sig_methyl_EPIC_ADC),
                                 all.cpg = rownames(meta_methyl_EPIC_ADC),
                                 collection = c("GO", "KEGG"),
                                 array.type = c("450K", "EPIC"),
                                 sig.genes = FALSE)

##### gometh EPIC pos #####
gometh_res_EPIC_blueC_pos = gometh(sig.cpg = rownames(sig_methyl_EPIC_blueC_pos),
                               all.cpg = rownames(meta_methyl_EPIC_ADC),
                               collection = c("GO", "KEGG"),
                               array.type = c("450K", "EPIC"),
                               sig.genes = FALSE)

gometh_res_EPIC_redC_pos = gometh(sig.cpg = rownames(sig_methyl_EPIC_redC_pos),
                              all.cpg = rownames(meta_methyl_EPIC_ADC),
                              collection = c("GO", "KEGG"),
                              array.type = c("450K", "EPIC"),
                              sig.genes = FALSE)
gometh_res_EPIC_redblue_pos = gometh(sig.cpg = rownames(sig_methyl_EPIC_redblue_pos),
                                 all.cpg = rownames(meta_methyl_EPIC_ADC),
                                 collection = c("GO", "KEGG"),
                                 array.type = c("450K", "EPIC"),
                                 sig.genes = FALSE)
gometh_res_EPIC_ADC_pos = gometh(sig.cpg = rownames(sig_methyl_EPIC_ADC_pos),
                             all.cpg = rownames(meta_methyl_EPIC_ADC),
                             collection = c("GO", "KEGG"),
                             array.type = c("450K", "EPIC"),
                             sig.genes = FALSE)
##### gometh EPIC neg #####

gometh_res_EPIC_blueC_neg = gometh(sig.cpg = rownames(sig_methyl_EPIC_blueC_neg),
                               all.cpg = rownames(meta_methyl_EPIC_ADC),
                               collection = c("GO", "KEGG"),
                               array.type = c("450K", "EPIC"),
                               sig.genes = FALSE)

gometh_res_EPIC_redC_neg = gometh(sig.cpg = rownames(sig_methyl_EPIC_redC_neg),
                              all.cpg = rownames(meta_methyl_EPIC_ADC),
                              collection = c("GO", "KEGG"),
                              array.type = c("450K", "EPIC"),
                              sig.genes = FALSE)
gometh_res_EPIC_redblue_neg = gometh(sig.cpg = rownames(sig_methyl_EPIC_redblue_neg),
                                 all.cpg = rownames(meta_methyl_EPIC_ADC),
                                 collection = c("GO", "KEGG"),
                                 array.type = c("450K", "EPIC"),
                                 sig.genes = FALSE)
gometh_res_EPIC_ADC_neg = gometh(sig.cpg = rownames(sig_methyl_EPIC_ADC_neg),
                             all.cpg = rownames(meta_methyl_EPIC_ADC),
                             collection = c("GO", "KEGG"),
                             array.type = c("450K", "EPIC"),
                             sig.genes = FALSE)

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
##### celltype enrichment #####

prev_study_res = read.table("41467_2021_23243_MOESM3_ESM.txt", sep = "\t", header = TRUE)
ensg_to_gene = read.table("RNAseq/mart_export.txt", sep = "\t", header = TRUE)


EPIC_manifest = read.csv("./MethylationEPIC_v-1-0_B4.csv", header = T)
# EPIC_manifest = EPIC_manifest[-c(1:6),]
# colnames(EPIC_manifest) = EPIC_manifest[1,]
# EPIC_manifest = EPIC_manifest[-1,]

target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_450k_ADC)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)

unconditional_results_meth_450kADC <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_450kADC$results, 
                mtc_method = "BH") 


target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_EPIC_ADC)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)

unconditional_results_meth_EPICADC <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_EPICADC$results, 
                mtc_method = "BH") 
#
target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_EPIC_redC)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)


target_list_meth = target_list_meth[target_list_meth %in% targetList_names_redC]
unconditional_results_meth_EPICredC <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_EPICredC$results, 
                mtc_method = "BH") 

target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_450k_redC)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)

target_list_meth = target_list_meth[target_list_meth %in% targetList_names_redC]
unconditional_results_meth_450kredC <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_450kredC$results, 
                mtc_method = "BH") 
#

target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_EPIC_blueC)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)


target_list_meth = target_list_meth[target_list_meth %in% targetList_names_blueC]
unconditional_results_meth_EPICblueC <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_EPICblueC$results, 
                mtc_method = "BH") 

target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_450k_blueC)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)

target_list_meth = unique(target_list_meth[target_list_meth %in% targetList_names_blueC])
unconditional_results_meth_450kblueC <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_450kblueC$results, 
                mtc_method = "BH") 
#

target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_EPIC_redblue)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)

unconditional_results_meth_EPICredblue <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_EPICredblue$results, 
                mtc_method = "BH") 

target_list_meth = EPIC_manifest[EPIC_manifest$Name %in% intersect(EPIC_manifest$Name, rownames(sig_methyl_450k_redblue)),]
target_list_meth = unique(target_list_meth$GencodeV41_Name)
target_list_meth = unique(unlist(strsplit(target_list_meth, split = ";")))
target_list_meth = gsub("\\.[0-9][0-9]", "", target_list_meth)
target_list_meth = gsub("\\.[0-9]", "", target_list_meth)
ensg_meth = target_list_meth[grepl("ENSG", target_list_meth)]
ensg_meth = unique(ensg_to_gene[ensg_to_gene$Gene.stable.ID %in% ensg_meth,]$Gene.name)
target_list_meth = c(target_list_meth[!grepl("ENSG", target_list_meth)], ensg_meth)

unconditional_results_meth_450kredblue <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = target_list_meth,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100000,
  annotLevel = 1)
EWCE::ewce_plot(unconditional_results_meth_450kredblue$results, 
                mtc_method = "BH") 

##### EPIC Coloc #####

load("D:/valentin/main/sig_betas_to_mQTL.Rdata")

mQTL_PITT = read.csv("PITT.eQTL.Cis.Pval.0.05.Dist.250000.csv")
mQTL_UKBBN = read.csv("UKBBN.eQTL.Cis.Pval.0.5.Dist.250000.csv")

MAF_PITT = read.table("maf/PITT.Imputed.hg19.frq", header = T)
MAF_UKBBN = read.table("maf/UKBBN.frq", header = T)

mQTL_PITT = mQTL_PITT[mQTL_PITT$pvalue < 1e-5,]
# mQTL_PITT = mQTL_PITT[mQTL_PITT$FDR < 5e-2,]

mQTL_UKBBN = mQTL_UKBBN[mQTL_UKBBN$pvalue < 1e-5,]
# mQTL_UKBBN = mQTL_UKBBN[mQTL_UKBBN$FDR < 5e-2,]

#####

mQTL_PITT$se <- mQTL_PITT$beta/abs(mQTL_PITT$statistic)
mQTL_PITT$snps <- str_split(mQTL_PITT$snps,pattern="_",simplify=T)[,1]
mQTL_PITT$BP <- str_split(mQTL_PITT$snps,pattern=":",simplify=T)[,2]
mQTL_PITT$CHR <- str_split(mQTL_PITT$snps,pattern=":",simplify=T)[,1]

mQTL_UKBBN$se <- mQTL_UKBBN$beta/abs(mQTL_UKBBN$statistic)
mQTL_UKBBN$snps <- str_split(mQTL_UKBBN$snps,pattern="_",simplify=T)[,1]
mQTL_UKBBN$BP <- str_split(mQTL_UKBBN$snps,pattern=":",simplify=T)[,2]
mQTL_UKBBN$CHR <- str_split(mQTL_UKBBN$snps,pattern=":",simplify=T)[,1]

mQTL_UKBBN = convert_loc_to_rs(mQTL_UKBBN, SNPlocs.Hsapiens.dbSNP144.GRCh37)
mQTL_PITT = convert_loc_to_rs(mQTL_PITT, SNPlocs.Hsapiens.dbSNP144.GRCh37)

index <- match(mQTL_UKBBN$snps, MAF_UKBBN$SNP)
mQTL_UKBBN$MAF <- MAF_UKBBN$MAF[index]

index <- match(mQTL_PITT$snps, MAF_PITT$SNP)
mQTL_PITT$MAF <- MAF_PITT$MAF[index]

mQTL_PITT = mQTL_PITT[!is.na(mQTL_PITT$MAF),]
mQTL_UKBBN = mQTL_UKBBN[!is.na(mQTL_UKBBN$MAF),]

mQTL_PITT = mQTL_PITT[mQTL_PITT$gene %in% intersect(mQTL_PITT$gene, mQTL_UKBBN$gene),]
mQTL_UKBBN = mQTL_UKBBN[mQTL_UKBBN$gene %in% intersect(mQTL_PITT$gene, mQTL_UKBBN$gene),]

mQTL_PITT_blueC = mQTL_PITT[mQTL_PITT$gene %in% rownames(sig_methyl_EPIC_blueC),]
mQTL_PITT_blueC$N = nrow(pheno_PITT[pheno_PITT$subtype %in% c("blue", "Control"),])

mQTL_PITT_redC = mQTL_PITT[mQTL_PITT$gene %in% rownames(sig_methyl_EPIC_redC),]
mQTL_PITT_redC$N = nrow(pheno_PITT[pheno_PITT$subtype %in% c("red", "Control"),])

mQTL_PITT_redblue = mQTL_PITT[mQTL_PITT$gene %in% rownames(sig_methyl_EPIC_redblue),]
mQTL_PITT_redblue$N = nrow(pheno_PITT[pheno_PITT$subtype %in% c("blue", "red"),])

mQTL_UKBBN_blueC = mQTL_UKBBN[mQTL_UKBBN$gene %in% rownames(sig_methyl_EPIC_blueC),]
mQTL_UKBBN_blueC$N = nrow(pheno_UKBBN[pheno_UKBBN$subtype %in% c("blue", "Control"),])

mQTL_UKBBN_redC = mQTL_UKBBN[mQTL_UKBBN$gene %in% rownames(sig_methyl_EPIC_redC),]
mQTL_UKBBN_redC$N = nrow(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red", "Control"),])

mQTL_UKBBN_redblue = mQTL_UKBBN[mQTL_UKBBN$gene %in% rownames(sig_methyl_EPIC_redblue),]
mQTL_UKBBN_redblue$N = nrow(pheno_UKBBN[pheno_UKBBN$subtype %in% c("blue", "red"),])


intersect(mQTL_PITT_blueC$snps, mQTL_UKBBN_blueC$snps)

unique(mQTL_PITT_blueC$gene)
unique(mQTL_UKBBN_blueC$gene)

write.csv(mQTL_PITT_blueC, file = "./mQTL_PITT_blueC_EPIC.csv")
write.csv(mQTL_PITT_redC, file = "./mQTL_PITT_redC_EPIC.csv")
write.csv(mQTL_PITT_redblue, file = "./mQTL_PITT_redblue_EPIC.csv")

write.csv(mQTL_UKBBN_blueC, file = "./mQTL_UKBBN_blueC_EPIC.csv")
write.csv(mQTL_UKBBN_redC, file = "./mQTL_UKBBN_redC_EPIC.csv")
write.csv(mQTL_UKBBN_redblue, file = "./mQTL_UKBBN_redblue_EPIC.csv")


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


##### Celltype correlation #####
load("betas/UKBBN_betas_only.Rdata")
betas_cohort_UKBBN = betas_cohort
load("betas/PITT_betas_only.Rdata")
betas_cohort_PITT = betas_cohort
load("betas/ROSMAP_betas_only.Rdata")
betas_cohort_ROSMAP = betas_cohort
remove(betas_cohort)

red_cohort_PITT = betas_cohort_PITT[rownames(betas_cohort_PITT) %in% rownames(sig_methyl_450k_redC),
                                    colnames(betas_cohort_PITT) %in% pheno_PITT[pheno_PITT$subtype == "red",]$Basename]
blue_cohort_PITT = betas_cohort_PITT[rownames(betas_cohort_PITT) %in% rownames(sig_methyl_450k_blueC),
                                    colnames(betas_cohort_PITT) %in% pheno_PITT[pheno_PITT$subtype == "blue",]$Basename]

red_cohort_UKBBN = betas_cohort_UKBBN[rownames(betas_cohort_UKBBN) %in% rownames(sig_methyl_450k_redC),
                                    colnames(betas_cohort_UKBBN) %in% pheno_UKBBN[pheno_UKBBN$subtype == "red",]$Basename]
blue_cohort_UKBBN = betas_cohort_UKBBN[rownames(betas_cohort_UKBBN) %in% rownames(sig_methyl_450k_blueC),
                                    colnames(betas_cohort_UKBBN) %in% pheno_UKBBN[pheno_UKBBN$subtype == "blue",]$Basename]

red_cohort_ROSMAP = betas_cohort_ROSMAP[rownames(betas_cohort_ROSMAP) %in% rownames(sig_methyl_450k_redC),
                                    colnames(betas_cohort_ROSMAP) %in% pheno_ROSMAP[pheno_ROSMAP$subtype == "red",]$Basename]
blue_cohort_ROSMAP = betas_cohort_ROSMAP[rownames(betas_cohort_ROSMAP) %in% rownames(sig_methyl_450k_blueC),
                                    colnames(betas_cohort_ROSMAP) %in% pheno_ROSMAP[pheno_ROSMAP$subtype == "blue",]$Basename]

remove(betas_cohort_PITT, betas_cohort_UKBBN, betas_cohort_ROSMAP)

red_betas = cbind(red_cohort_PITT, red_cohort_ROSMAP, red_cohort_UKBBN)
blue_betas = cbind(blue_cohort_PITT, blue_cohort_ROSMAP, blue_cohort_UKBBN)

cellSet_PITT_red = pheno_PITT[pheno_PITT$subtype == "red",c("Basename","Neuronal_Inhibitory", "Neuronal_Excitatory", "Microglia", "Astrocyte", "Oligodendrocyte")]
cellSet_UKBBN_red = pheno_UKBBN[pheno_UKBBN$subtype == "red",c("Basename","Neuronal_Inhibitory", "Neuronal_Excitatory", "Microglia", "Astrocyte", "Oligodendrocyte")]
cellSet_ROSMAP_red = pheno_ROSMAP[pheno_ROSMAP$subtype == "red",c("Basename","Neuronal_Inhibitory", "Neuronal_Excitatory", "Microglia", "Astrocyte", "Oligodendrocyte")]

cellSet_PITT_blue = pheno_PITT[pheno_PITT$subtype == "blue",c("Basename","Neuronal_Inhibitory", "Neuronal_Excitatory", "Microglia", "Astrocyte", "Oligodendrocyte")]
cellSet_UKBBN_blue = pheno_UKBBN[pheno_UKBBN$subtype == "blue",c("Basename","Neuronal_Inhibitory", "Neuronal_Excitatory", "Microglia", "Astrocyte", "Oligodendrocyte")]
cellSet_ROSMAP_blue = pheno_ROSMAP[pheno_ROSMAP$subtype == "blue",c("Basename","Neuronal_Inhibitory", "Neuronal_Excitatory", "Microglia", "Astrocyte", "Oligodendrocyte")]

#####
#Create cor matrix with number of Cpgs as rows and number of cell types as cols
cellCor_red <- matrix(nrow = nrow(red_betas), ncol = 5)
colnames(cellCor_red) <- colnames(cellSet_red)[-1]
rownames(cellCor_red) <- rownames(red_betas)

x <- rownames(red_betas)[1]
for(x in rownames(red_betas)){
  corNeunI <- cor(red_betas[x,],cellSet_red$Neuronal_Inhibitory)
  corNeunE <- cor(red_betas[x,],cellSet_red$Neuronal_Excitatory)
  corMic <- cor(red_betas[x,],cellSet_red$Microglia)
  corAst <- cor(red_betas[x,],cellSet_red$Astrocyte)
  corOlig <- cor(red_betas[x,],cellSet_red$Oligodendrocyte)
  cellCor_red[x,] <- c(corNeunI,corNeunE,corMic,corAst,corOlig)
}

cellCor_red = abs(cellCor_red)
cellCorLong_red <- melt(cellCor_red)
cellCorLong_red$subtype = "red"

cellCor_blue <- matrix(nrow = nrow(blue_betas), ncol = 5)
colnames(cellCor_blue) <- colnames(cellSet_blue)[-1]
rownames(cellCor_blue) <- rownames(blue_betas)

x <- rownames(blue_betas)[1]
for(x in rownames(blue_betas)){
  corNeunI <- cor(blue_betas[x,],cellSet_blue$Neuronal_Inhibitory)
  corNeunE <- cor(blue_betas[x,],cellSet_blue$Neuronal_Excitatory)
  corMic <- cor(blue_betas[x,],cellSet_blue$Microglia)
  corAst <- cor(blue_betas[x,],cellSet_blue$Astrocyte)
  corOlig <- cor(blue_betas[x,],cellSet_blue$Oligodendrocyte)
  cellCor_blue[x,] <- c(corNeunI,corNeunE,corMic,corAst,corOlig)
}

cellCor_blue = abs(cellCor_blue)
cellCorLong_blue <- melt(cellCor_blue)
cellCorLong_blue$subtype = "blue"


png("./betas/cor_red.png")
heatmap(cellCor_red)
dev.off()

png("./betas/box_red.png")
boxplot(cellCorLong_red$value~as.factor(cellCorLong_red$Var2))
dev.off()

png("./betas/cor_blue.png")
heatmap(cellCor_blue)
dev.off()

png("./betas/box_blue.png")
boxplot(cellCorLong_blue$value~as.factor(cellCorLong_blue$Var2))
dev.off()

cellCorLong_all = rbind(cellCorLong_blue, cellCorLong_red)
cellCorLong_all$cellsubtype = paste(cellCorLong_all$subtype, cellCorLong_all$Var2, sep = "_")

cellCorLong_all$cellsubtype = factor(cellCorLong_all$cellsubtype,
                                     levels = c("blue_Neuronal_Inhibitory",
                                                "red_Neuronal_Inhibitory",
                                                "blue_Neuronal_Excitatory",
                                                "red_Neuronal_Excitatory",
                                                "blue_Microglia",
                                                "red_Microglia",
                                                "blue_Astrocyte",
                                                "red_Astrocyte",
                                                "blue_Oligodendrocyte",
                                                "red_Oligodendrocyte"))
boxplot(cellCorLong_all$value~cellCorLong_all$cellsubtype)

comparisons <- list( c("blue_Neuronal_Inhibitory", "red_Neuronal_Inhibitory"), 
                        c("blue_Neuronal_Excitatory", "red_Neuronal_Excitatory"), 
                        c("blue_Microglia", "red_Microglia"),
                        c("blue_Astrocyte", "red_Astrocyte"),
                        c("blue_Oligodendrocyte", "red_Oligodendrocyte"))

ggplot(cellCorLong_all, aes(x=cellsubtype, y=value)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test", comparisons = comparisons) +
  theme_bw()


boxplot(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Neuronal_Inhibitory~
          pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Neuronal_Inhibitory,
           as.numeric(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Neuronal_Inhibitory~
          pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Neuronal_Inhibitory,
         as.numeric(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Neuronal_Inhibitory~
          pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Neuronal_Inhibitory,
         as.numeric(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Neuronal_Excitatory~
          pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Neuronal_Excitatory,
         as.numeric(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Neuronal_Excitatory~
          pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Neuronal_Excitatory,
         as.numeric(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Neuronal_Excitatory~
          pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Neuronal_Excitatory,
         as.numeric(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Microglia~
          pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Microglia,
         as.numeric(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Microglia~
          pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Microglia,
         as.numeric(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Microglia~
          pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Microglia,
         as.numeric(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Astrocyte~
          pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Astrocyte,
         as.numeric(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Astrocyte~
          pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Astrocyte,
         as.numeric(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Astrocyte~
          pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Astrocyte,
         as.numeric(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Oligodendrocyte~
          pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$Oligodendrocyte,
         as.numeric(pheno_PITT[pheno_PITT$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Oligodendrocyte~
          pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$Oligodendrocyte,
         as.numeric(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")

boxplot(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Oligodendrocyte~
          pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$subtype)
cor.test(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$Oligodendrocyte,
         as.numeric(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red","blue"),]$HC_clusters), method = "spearman")
#####

# Per cohort
# Correlation Celltype vs Subtype - Make sure subtypes and celltype are not similar - OK
# Celltype specificity Subtype vs all samples - Once to check for each subtype
# Random correlation sampling vs true signatures x100 + Add a line for the median of all the rest

betas_cohort_PITT_blue = betas_cohort_PITT[, colnames(betas_cohort_PITT) %in% pheno_PITT[pheno_PITT$subtype == "blue",]$Basename]
betas_cohort_PITT_red = betas_cohort_PITT[, colnames(betas_cohort_PITT) %in% pheno_PITT[pheno_PITT$subtype == "red",]$Basename]

betas_cohort_UKBBN_blue = betas_cohort_UKBBN[, colnames(betas_cohort_UKBBN) %in% pheno_UKBBN[pheno_UKBBN$subtype == "blue",]$Basename]
betas_cohort_UKBBN_red = betas_cohort_UKBBN[, colnames(betas_cohort_UKBBN) %in% pheno_UKBBN[pheno_UKBBN$subtype == "red",]$Basename]

betas_cohort_ROSMAP_blue = betas_cohort_ROSMAP[, colnames(betas_cohort_ROSMAP) %in% pheno_ROSMAP[pheno_ROSMAP$subtype == "blue",]$Basename]
betas_cohort_ROSMAP_red = betas_cohort_ROSMAP[, colnames(betas_cohort_ROSMAP) %in% pheno_ROSMAP[pheno_ROSMAP$subtype == "red",]$Basename]

#
corNeunI = pbapply::pbapply(betas_cohort_PITT_blue, 1, function(x){cor(x, cellSet_PITT_blue$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_PITT_blue, 1, function(x){cor(x, cellSet_PITT_blue$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_PITT_blue, 1, function(x){cor(x, cellSet_PITT_blue$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_PITT_blue, 1, function(x){cor(x, cellSet_PITT_blue$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_PITT_blue, 1, function(x){cor(x, cellSet_PITT_blue$Oligodendrocyte)})

cellCor_Pblue_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_Pblue_all) <- colnames(cellSet_PITT_blue)[-1]
rownames(cellCor_Pblue_all) <- rownames(betas_cohort_PITT_blue)

corNeunI = pbapply::pbapply(betas_cohort_PITT_red, 1, function(x){cor(x, cellSet_PITT_red$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_PITT_red, 1, function(x){cor(x, cellSet_PITT_red$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_PITT_red, 1, function(x){cor(x, cellSet_PITT_red$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_PITT_red, 1, function(x){cor(x, cellSet_PITT_red$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_PITT_red, 1, function(x){cor(x, cellSet_PITT_red$Oligodendrocyte)})

cellCor_Pred_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_Pred_all) <- colnames(cellSet_PITT_red)[-1]
rownames(cellCor_Pred_all) <- rownames(betas_cohort_PITT_red)

#
corNeunI = pbapply::pbapply(betas_cohort_UKBBN_blue, 1, function(x){cor(x, cellSet_UKBBN_blue$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_UKBBN_blue, 1, function(x){cor(x, cellSet_UKBBN_blue$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_UKBBN_blue, 1, function(x){cor(x, cellSet_UKBBN_blue$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_UKBBN_blue, 1, function(x){cor(x, cellSet_UKBBN_blue$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_UKBBN_blue, 1, function(x){cor(x, cellSet_UKBBN_blue$Oligodendrocyte)})

cellCor_Ublue_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_Ublue_all) <- colnames(cellSet_UKBBN_blue)[-1]
rownames(cellCor_Ublue_all) <- rownames(betas_cohort_UKBBN_blue)

corNeunI = pbapply::pbapply(betas_cohort_UKBBN_red, 1, function(x){cor(x, cellSet_UKBBN_red$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_UKBBN_red, 1, function(x){cor(x, cellSet_UKBBN_red$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_UKBBN_red, 1, function(x){cor(x, cellSet_UKBBN_red$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_UKBBN_red, 1, function(x){cor(x, cellSet_UKBBN_red$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_UKBBN_red, 1, function(x){cor(x, cellSet_UKBBN_red$Oligodendrocyte)})

cellCor_Ured_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_Ured_all) <- colnames(cellSet_UKBBN_red)[-1]
rownames(cellCor_Ured_all) <- rownames(betas_cohort_UKBBN_red)


#
corNeunI = pbapply::pbapply(betas_cohort_ROSMAP_blue, 1, function(x){cor(x, cellSet_ROSMAP_blue$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_ROSMAP_blue, 1, function(x){cor(x, cellSet_ROSMAP_blue$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_ROSMAP_blue, 1, function(x){cor(x, cellSet_ROSMAP_blue$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_ROSMAP_blue, 1, function(x){cor(x, cellSet_ROSMAP_blue$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_ROSMAP_blue, 1, function(x){cor(x, cellSet_ROSMAP_blue$Oligodendrocyte)})

cellCor_Rblue_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_Rblue_all) <- colnames(cellSet_ROSMAP_blue)[-1]
rownames(cellCor_Rblue_all) <- rownames(betas_cohort_ROSMAP_blue)

corNeunI = pbapply::pbapply(betas_cohort_ROSMAP_red, 1, function(x){cor(x, cellSet_ROSMAP_red$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_ROSMAP_red, 1, function(x){cor(x, cellSet_ROSMAP_red$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_ROSMAP_red, 1, function(x){cor(x, cellSet_ROSMAP_red$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_ROSMAP_red, 1, function(x){cor(x, cellSet_ROSMAP_red$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_ROSMAP_red, 1, function(x){cor(x, cellSet_ROSMAP_red$Oligodendrocyte)})

cellCor_Rred_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_Rred_all) <- colnames(cellSet_ROSMAP_red)[-1]
rownames(cellCor_Rred_all) <- rownames(betas_cohort_ROSMAP_red)

dc <- detectCores()-1
registerDoParallel(dc)

###### PITT ######
#### Blue ####
cellCor_Pblue_orig = cellCor_Pblue_all[rownames(cellCor_Pblue_all) %in% 
                                         rownames(sig_methyl_450k_blueC),]
cellCor_Pblue_orig = abs(cellCor_Pblue_orig)
cellCor_Pblue_random = foreach(i = seq(1:100), .combine = cbind) %dopar% {
  rows = sample(seq(1:nrow(cellCor_Pblue_all)), size = nrow(cellCor_Pblue_orig))
  tmp_cellCor = cellCor_Pblue_all[rows,]
  
  colnames(tmp_cellCor) = paste(i, colnames(tmp_cellCor), sep = "-")
  
  tmp_cellCor
}
cellCor_Pblue_random = abs(cellCor_Pblue_random)

cellCor_Pblue = cbind(cellCor_Pblue_orig,cellCor_Pblue_random)

cellCor_Pblue_long = melt(cellCor_Pblue)
cellCor_Pblue_long$Var2 = as.character(cellCor_Pblue_long$Var2)
cellCor_Pblue_long$celltype = gsub("[0-9]*-", "",cellCor_Pblue_long$Var2)

for(celltype in unique(cellCor_Pblue_long$celltype)){
  tmp_long = cellCor_Pblue_long[cellCor_Pblue_long$celltype == celltype,]
  tmp_long$group = "Sig"
  tmp_long[grep("[0-9]",tmp_long$Var2),]$group = "Random"
  tmp_long$group = factor(tmp_long$group, levels = c("Sig", "Random"))
  
  tmp_long = tmp_long[complete.cases(tmp_long),]
  
  med_random = median(tmp_long[tmp_long$group == "Random",]$value)
  
  png(paste0("Celltype_cor/PITT/blue", celltype, ".png"), width = 1200, height = 800, units = "px")
  plot(ggplot(tmp_long, aes(x=Var2, y=value, fill = group)) +
      geom_boxplot() +
      geom_hline(yintercept=med_random, 
                 color = "green", size=2) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}
#### Red ####

cellCor_Pred_orig = cellCor_Pred_all[rownames(cellCor_Pred_all) %in% 
                                         rownames(sig_methyl_450k_redC),]
cellCor_Pred_orig = abs(cellCor_Pred_orig)
cellCor_Pred_random = foreach(i = seq(1:100), .combine = cbind) %dopar% {
  rows = sample(seq(1:nrow(cellCor_Pred_all)), size = nrow(cellCor_Pred_orig))
  tmp_cellCor = cellCor_Pred_all[rows,]
  
  colnames(tmp_cellCor) = paste(i, colnames(tmp_cellCor), sep = "-")
  
  tmp_cellCor
}
cellCor_Pred_random = abs(cellCor_Pred_random)

cellCor_Pred = cbind(cellCor_Pred_orig,cellCor_Pred_random)

cellCor_Pred_long = melt(cellCor_Pred)
cellCor_Pred_long$Var2 = as.character(cellCor_Pred_long$Var2)
cellCor_Pred_long$celltype = gsub("[0-9]*-", "",cellCor_Pred_long$Var2)

for(celltype in unique(cellCor_Pred_long$celltype)){
  tmp_long = cellCor_Pred_long[cellCor_Pred_long$celltype == celltype,]
  tmp_long$group = "Sig"
  tmp_long[grep("[0-9]",tmp_long$Var2),]$group = "Random"
  tmp_long$group = factor(tmp_long$group, levels = c("Sig", "Random"))
  
  tmp_long = tmp_long[complete.cases(tmp_long),]
  
  med_random = median(tmp_long[tmp_long$group == "Random",]$value)
  
  png(paste0("Celltype_cor/PITT/red", celltype, ".png"), width = 1200, height = 800, units = "px")
  plot(ggplot(tmp_long, aes(x=Var2, y=value, fill = group)) +
         geom_boxplot() +
         geom_hline(yintercept=med_random, 
                    color = "green", size=2) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}


###### UKBBN ######
#### Blue ####
cellCor_Ublue_orig = cellCor_Ublue_all[rownames(cellCor_Ublue_all) %in% 
                                         rownames(sig_methyl_450k_blueC),]
cellCor_Ublue_orig = abs(cellCor_Ublue_orig)
cellCor_Ublue_random = foreach(i = seq(1:100), .combine = cbind) %dopar% {
  rows = sample(seq(1:nrow(cellCor_Ublue_all)), size = nrow(cellCor_Ublue_orig))
  tmp_cellCor = cellCor_Ublue_all[rows,]
  
  colnames(tmp_cellCor) = paste(i, colnames(tmp_cellCor), sep = "-")
  
  tmp_cellCor
}
cellCor_Ublue_random = abs(cellCor_Ublue_random)

cellCor_Ublue = cbind(cellCor_Ublue_orig,cellCor_Ublue_random)

cellCor_Ublue_long = melt(cellCor_Ublue)
cellCor_Ublue_long$Var2 = as.character(cellCor_Ublue_long$Var2)
cellCor_Ublue_long$celltype = gsub("[0-9]*-", "",cellCor_Ublue_long$Var2)

for(celltype in unique(cellCor_Ublue_long$celltype)){
  tmp_long = cellCor_Ublue_long[cellCor_Ublue_long$celltype == celltype,]
  tmp_long$group = "Sig"
  tmp_long[grep("[0-9]",tmp_long$Var2),]$group = "Random"
  tmp_long$group = factor(tmp_long$group, levels = c("Sig", "Random"))
  
  tmp_long = tmp_long[complete.cases(tmp_long),]
  
  med_random = median(tmp_long[tmp_long$group == "Random",]$value)
  
  png(paste0("Celltype_cor/UKBBN/blue", celltype, ".png"), width = 1200, height = 800, units = "px")
  plot(ggplot(tmp_long, aes(x=Var2, y=value, fill = group)) +
         geom_boxplot() +
         geom_hline(yintercept=med_random, 
                    color = "green", size=2) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}

#### Red ####
cellCor_Ured_orig = cellCor_Ured_all[rownames(cellCor_Ured_all) %in% 
                                         rownames(sig_methyl_450k_redC),]
cellCor_Ured_orig = abs(cellCor_Ured_orig)
cellCor_Ured_random = foreach(i = seq(1:100), .combine = cbind) %dopar% {
  rows = sample(seq(1:nrow(cellCor_Ured_all)), size = nrow(cellCor_Ured_orig))
  tmp_cellCor = cellCor_Ured_all[rows,]
  
  colnames(tmp_cellCor) = paste(i, colnames(tmp_cellCor), sep = "-")
  
  tmp_cellCor
}
cellCor_Ured_random = abs(cellCor_Ured_random)

cellCor_Ured = cbind(cellCor_Ured_orig,cellCor_Ured_random)

cellCor_Ured_long = melt(cellCor_Ured)
cellCor_Ured_long$Var2 = as.character(cellCor_Ured_long$Var2)
cellCor_Ured_long$celltype = gsub("[0-9]*-", "",cellCor_Ured_long$Var2)

for(celltype in unique(cellCor_Ured_long$celltype)){
  tmp_long = cellCor_Ured_long[cellCor_Ured_long$celltype == celltype,]
  tmp_long$group = "Sig"
  tmp_long[grep("[0-9]",tmp_long$Var2),]$group = "Random"
  tmp_long$group = factor(tmp_long$group, levels = c("Sig", "Random"))
  
  tmp_long = tmp_long[complete.cases(tmp_long),]
  
  med_random = median(tmp_long[tmp_long$group == "Random",]$value)
  
  png(paste0("Celltype_cor/UKBBN/red", celltype, ".png"), width = 1200, height = 800, units = "px")
  plot(ggplot(tmp_long, aes(x=Var2, y=value, fill = group)) +
         geom_boxplot() +
         geom_hline(yintercept=med_random, 
                    color = "green", size=2) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}

###### ROSMAP ######
#### Blue ####
cellCor_Rblue_orig = cellCor_Rblue_all[rownames(cellCor_Rblue_all) %in% 
                                         rownames(sig_methyl_450k_blueC),]
cellCor_Rblue_orig = abs(cellCor_Rblue_orig)
cellCor_Rblue_random = foreach(i = seq(1:100), .combine = cbind) %dopar% {
  rows = sample(seq(1:nrow(cellCor_Rblue_all)), size = nrow(cellCor_Rblue_orig))
  tmp_cellCor = cellCor_Rblue_all[rows,]
  
  colnames(tmp_cellCor) = paste(i, colnames(tmp_cellCor), sep = "-")
  
  tmp_cellCor
}
cellCor_Rblue_random = abs(cellCor_Rblue_random)

cellCor_Rblue = cbind(cellCor_Rblue_orig,cellCor_Rblue_random)

cellCor_Rblue_long = melt(cellCor_Rblue)
cellCor_Rblue_long$Var2 = as.character(cellCor_Rblue_long$Var2)
cellCor_Rblue_long$celltype = gsub("[0-9]*-", "",cellCor_Rblue_long$Var2)

for(celltype in unique(cellCor_Rblue_long$celltype)){
  tmp_long = cellCor_Rblue_long[cellCor_Rblue_long$celltype == celltype,]
  tmp_long$group = "Sig"
  tmp_long[grep("[0-9]",tmp_long$Var2),]$group = "Random"
  tmp_long$group = factor(tmp_long$group, levels = c("Sig", "Random"))
  
  tmp_long = tmp_long[complete.cases(tmp_long),]
  
  med_random = median(tmp_long[tmp_long$group == "Random",]$value)
  
  png(paste0("Celltype_cor/ROSMAP/blue", celltype, ".png"), width = 1200, height = 800, units = "px")
  plot(ggplot(tmp_long, aes(x=Var2, y=value, fill = group)) +
         geom_boxplot() +
         geom_hline(yintercept=med_random, 
                    color = "green", size=2) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}
#### Red ####
cellCor_Rred_orig = cellCor_Rred_all[rownames(cellCor_Rred_all) %in% 
                                         rownames(sig_methyl_450k_blueC),]
cellCor_Rred_orig = abs(cellCor_Rred_orig)
cellCor_Rred_random = foreach(i = seq(1:100), .combine = cbind) %dopar% {
  rows = sample(seq(1:nrow(cellCor_Rred_all)), size = nrow(cellCor_Rred_orig))
  tmp_cellCor = cellCor_Rred_all[rows,]
  
  colnames(tmp_cellCor) = paste(i, colnames(tmp_cellCor), sep = "-")
  
  tmp_cellCor
}
cellCor_Rred_random = abs(cellCor_Rred_random)

cellCor_Rred = cbind(cellCor_Rred_orig,cellCor_Rred_random)

cellCor_Rred_long = melt(cellCor_Rred)
cellCor_Rred_long$Var2 = as.character(cellCor_Rred_long$Var2)
cellCor_Rred_long$celltype = gsub("[0-9]*-", "",cellCor_Rred_long$Var2)

for(celltype in unique(cellCor_Rred_long$celltype)){
  tmp_long = cellCor_Rred_long[cellCor_Rred_long$celltype == celltype,]
  tmp_long$group = "Sig"
  tmp_long[grep("[0-9]",tmp_long$Var2),]$group = "Random"
  tmp_long$group = factor(tmp_long$group, levels = c("Sig", "Random"))
  
  tmp_long = tmp_long[complete.cases(tmp_long),]
  
  med_random = median(tmp_long[tmp_long$group == "Random",]$value)
  
  png(paste0("Celltype_cor/ROSMAP/red", celltype, ".png"), width = 1200, height = 800, units = "px")
  plot(ggplot(tmp_long, aes(x=Var2, y=value, fill = group)) +
         geom_boxplot() +
         geom_hline(yintercept=med_random, 
                    color = "green", size=2) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  }

##### All #####

betas_cohort_blue = cbind(betas_cohort_PITT_blue, betas_cohort_UKBBN_blue, betas_cohort_ROSMAP_blue)
betas_cohort_red = cbind(betas_cohort_PITT_red, betas_cohort_UKBBN_red, betas_cohort_ROSMAP_red)
cellSet_blue = rbind(cellSet_PITT_blue,cellSet_UKBBN_blue,cellSet_ROSMAP_blue)
cellSet_red = rbind(cellSet_PITT_red,cellSet_UKBBN_red,cellSet_ROSMAP_red)

corNeunI = pbapply::pbapply(betas_cohort_blue, 1, function(x){cor(x, cellSet_blue$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_blue, 1, function(x){cor(x, cellSet_blue$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_blue, 1, function(x){cor(x, cellSet_blue$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_blue, 1, function(x){cor(x, cellSet_blue$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_blue, 1, function(x){cor(x, cellSet_blue$Oligodendrocyte)})

cellCor_blue_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_blue_all) <- colnames(cellSet_blue)[-1]
rownames(cellCor_blue_all) <- rownames(betas_cohort_blue)

corNeunI = pbapply::pbapply(betas_cohort_red, 1, function(x){cor(x, cellSet_red$Neuronal_Inhibitory)})
corNeunE = pbapply::pbapply(betas_cohort_red, 1, function(x){cor(x, cellSet_red$Neuronal_Excitatory)})
corMic = pbapply::pbapply(betas_cohort_red, 1, function(x){cor(x, cellSet_red$Microglia)})  
corAst = pbapply::pbapply(betas_cohort_red, 1, function(x){cor(x, cellSet_red$Astrocyte)})
corOlig = pbapply::pbapply(betas_cohort_red, 1, function(x){cor(x, cellSet_red$Oligodendrocyte)})

cellCor_red_all <- cbind(corNeunI,corNeunE,corMic,corAst,corOlig)
colnames(cellCor_red_all) <- colnames(cellSet_red)[-1]
rownames(cellCor_red_all) <- rownames(betas_cohort_red)


# save.image("RNAseq/current.Rdata")
##### FACS celltype enrichment #####

load("FACS/IRF8_allcpg.rdata")
AllLME_IRF8 = AllLME
load("FACS/NeuN_allcpg.rdata")
AllLME_NeuN = AllLME
load("FACS/Sox10_allcpg.rdata")
AllLME_Sox10 = AllLME
load("FACS/Trip neg_allcpg.rdata")
AllLME_TripNeg = AllLME
remove(AllLME)

AllLME_IRF8 = t(AllLME_IRF8)
AllLME_NeuN = t(AllLME_NeuN)
AllLME_Sox10 = t(AllLME_Sox10)
AllLME_TripNeg = t(AllLME_TripNeg)

head(AllLME_IRF8)
head(AllLME_NeuN)
head(AllLME_Sox10)
head(AllLME_TripNeg)

identical(rownames(AllLME_IRF8), rownames(AllLME_Sox10))

All_CT = cbind(AllLME_IRF8[,c(1,5)], AllLME_Sox10[,c(1,5)], AllLME_NeuN[,c(1,5)], AllLME_TripNeg[,c(1,5)])
head(All_CT)
colnames(All_CT) = c("Mic_est", "Mic_pval",
                     "Olig_est", "Olig_pval",
                     "NeuN_est", "NeuN_pval",
                     "Ast_est", "Ast_pval")
All_CT = All_CT[rownames(All_CT) %in% rownames(betas_ROSMAP),]
All_CT = as.data.frame(All_CT)

All_CT$Mic_pval = p.adjust(All_CT$Mic_pval, method = "bonf")
range(All_CT$Mic_est)

All_CT$Ast_pval = p.adjust(All_CT$Ast_pval, method = "bonf")
range(All_CT$Ast_est)

All_CT$NeuN_pval = p.adjust(All_CT$NeuN_pval, method = "bonf")
range(All_CT$NeuN_est)

All_CT$Olig_pval = p.adjust(All_CT$Olig_pval, method = "bonf")
range(All_CT$Olig_est)

All_CT_blue = All_CT[rownames(All_CT) %in% rownames(sig_methyl_450k_blueC),]

All_CT_red = All_CT[rownames(All_CT) %in% rownames(sig_methyl_450k_redC),]

All_CT_grey = All_CT[rownames(All_CT) %in% rownames(sig_methyl_450k_greyC),]

# save(All_CT, All_CT_red, All_CT_blue, All_CT_AD, file = "FACS/All_LME_FACS450k.rdata")
load("FACS/All_LME_FACS450k.rdata")

pvalue_threshold = 1e-10
effect_threshold = 0.5

AllLME_IRF8 = as.data.frame(AllLME_IRF8)
AllLME_IRF8 = AllLME_IRF8[AllLME_IRF8$"p-value" < pvalue_threshold & abs(AllLME_IRF8$Value) > effect_threshold,]
AllLME_IRF8 = AllLME_IRF8[rownames(AllLME_IRF8) %in% rownames(betas_ROSMAP),]
AllLME_IRF8 = AllLME_IRF8[order(AllLME_IRF8$`p-value`),]
AllLME_IRF8$p.bonf = p.adjust(AllLME_IRF8$`p-value`, method = "bonf")
write.csv(AllLME_IRF8, file = "SUBTYPING-PAPER-TABLEFIG/Sig_MIC_FANS.csv")

AllLME_NeuN = as.data.frame(AllLME_NeuN)
AllLME_NeuN = AllLME_NeuN[AllLME_NeuN$"p-value" < pvalue_threshold & abs(AllLME_NeuN$Value) > effect_threshold,]
AllLME_NeuN = AllLME_NeuN[rownames(AllLME_NeuN) %in% rownames(betas_ROSMAP),]
AllLME_NeuN = AllLME_NeuN[order(AllLME_NeuN$`p-value`),]
AllLME_NeuN$p.bonf = p.adjust(AllLME_NeuN$`p-value`, method = "bonf")
write.csv(AllLME_NeuN, file = "SUBTYPING-PAPER-TABLEFIG/Sig_NEUN_FANS.csv")

AllLME_Sox10 = as.data.frame(AllLME_Sox10)
AllLME_Sox10 = AllLME_Sox10[AllLME_Sox10$"p-value" < pvalue_threshold & abs(AllLME_Sox10$Value) > effect_threshold,]
AllLME_Sox10 = AllLME_Sox10[rownames(AllLME_Sox10) %in% rownames(betas_ROSMAP),]
AllLME_Sox10 = AllLME_Sox10[order(AllLME_Sox10$`p-value`),]
AllLME_Sox10$p.bonf = p.adjust(AllLME_Sox10$`p-value`, method = "bonf")
write.csv(AllLME_Sox10, file = "SUBTYPING-PAPER-TABLEFIG/Sig_OLIG_FANS.csv")

AllLME_TripNeg = as.data.frame(AllLME_TripNeg)
AllLME_TripNeg = AllLME_TripNeg[AllLME_TripNeg$"p-value" < pvalue_threshold & abs(AllLME_TripNeg$Value) > effect_threshold,]
AllLME_TripNeg = AllLME_TripNeg[rownames(AllLME_TripNeg) %in% rownames(betas_ROSMAP),]
AllLME_TripNeg = AllLME_TripNeg[order(AllLME_TripNeg$`p-value`),]
AllLME_TripNeg$p.bonf = p.adjust(AllLME_TripNeg$`p-value`, method = "bonf")
write.csv(AllLME_TripNeg, file = "SUBTYPING-PAPER-TABLEFIG/Sig_AST_FANS.csv")

table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)
table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)
table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)
table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)

table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)
table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)
table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)
table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)

table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)
table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)
table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)
table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)

table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)
table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)
table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)
table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)


table_res_fisher = data.frame(matrix(ncol = 6, nrow = 4))
colnames(table_res_fisher) = c("Blue_pval","Blue_odds", "Red_pval","Red_odds", "Grey_pval", "Grey_odds")
rownames(table_res_fisher) = c("NeuN","Ast","Mic","Oligo")

####  blue ####
table_blue_Olig = matrix(c(table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                           table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[1] - table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2] - table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1] + table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2]),
                         nrow = 2)
fisher_blue_Olig = fisher.test(table_blue_Olig, alternative = "greater")
table_res_fisher[4,1] = fisher_blue_Olig$p.value
table_res_fisher[4,2] = fisher_blue_Olig$estimate

table_blue_NeuN = matrix(c(table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[1] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1] + table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2]),
                         nrow = 2)
fisher_blue_NeuN = fisher.test(table_blue_NeuN, alternative = "greater")
table_res_fisher[1,1] = fisher_blue_NeuN$p.value
table_res_fisher[1,2] = fisher_blue_NeuN$estimate

table_blue_Ast = matrix(c(table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[2],
                           table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[1] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1] + table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2]),
                         nrow = 2)
fisher_blue_Ast = fisher.test(table_blue_Ast, alternative = "greater")
table_res_fisher[2,1] = fisher_blue_Ast$p.value
table_res_fisher[2,2] = fisher_blue_Ast$estimate

table_blue_Mic = matrix(c(table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                           table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[1] - table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                           table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2] - table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                           table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1]),
                         nrow = 2)
fisher_blue_Mic = fisher.test(table_blue_Mic, alternative = "greater")
table_res_fisher[3,1] = fisher_blue_Mic$p.value
table_res_fisher[3,2] = fisher_blue_Mic$estimate

#### red ####
table_red_Olig = matrix(c(table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                           table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[1] - table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2] - table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1]),
                         nrow = 2)
fisher_red_Olig = fisher.test(table_red_Olig, alternative = "greater")
table_res_fisher[4,3] = fisher_red_Olig$p.value
table_res_fisher[4,4] = fisher_red_Olig$estimate

table_red_NeuN = matrix(c(table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                           table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[1] - table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2] - table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1]),
                         nrow = 2)
fisher_red_NeuN = fisher.test(table_red_NeuN, alternative = "greater")
table_res_fisher[1,3] = fisher_red_NeuN$p.value
table_res_fisher[1,4] = fisher_red_NeuN$estimate

table_red_Ast = matrix(c(table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[2],
                          table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[1] - table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2] - table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1]),
                        nrow = 2)
# fisher_red_Ast = fisher.test(table_red_Ast, alternative = "greater")
table_res_fisher[2,3] = 0
table_res_fisher[2,4] = 0

table_red_Mic = matrix(c(table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                          table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[1] - table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2] - table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1]),
                        nrow = 2)
fisher_red_Mic = fisher.test(table_red_Mic, alternative = "greater")
table_res_fisher[3,3] = fisher_red_Mic$p.value
table_res_fisher[3,4] = fisher_red_Mic$estimate

#### grey ####
table_grey_Olig = matrix(c(table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2],
                           table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[1] - table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2] - table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1]),
                         nrow = 2)
fisher_grey_Olig = fisher.test(table_grey_Olig, alternative = "greater")
table_res_fisher[4,5] = fisher_grey_Olig$p.value
table_res_fisher[4,6] = fisher_grey_Olig$estimate

table_grey_NeuN = matrix(c(table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                           table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[1] - table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2] - table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1]),
                         nrow = 2)
fisher_grey_NeuN = fisher.test(table_grey_NeuN, alternative = "greater")
table_res_fisher[1,5] = fisher_grey_NeuN$p.value
table_res_fisher[1,6] = fisher_grey_NeuN$estimate

table_grey_Ast = matrix(c(table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[2],
                          table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[1] - table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2] - table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1]),
                        nrow = 2)
# fisher_grey_Ast = fisher.test(table_grey_Ast, alternative = "greater")
table_res_fisher[2,5] = 0
table_res_fisher[2,6] = 0

table_grey_Mic = matrix(c(table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                          table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[1] - table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2] - table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1]),
                        nrow = 2)
fisher_grey_Mic = fisher.test(table_grey_Mic)
table_res_fisher[3,5] = fisher_grey_Mic$p.value
table_res_fisher[3,6] = fisher_grey_Mic$estimate



# Bar

to_barplot = data.frame(matrix(ncol = 8, nrow = 12))
colnames(to_barplot) = c("subtype","cell_type","total_non_sig","total_sig",
                         "celltype_sig","celltype_non_sig","pvalue","odds")
to_barplot$subtype = c(rep("blue", 4), rep("red", 4), rep("grey", 4))
to_barplot$cell_type = rep(c("NeuN", "Ast", "Mic", "Olig"), 3)

to_barplot$total_non_sig = c(table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1],
                             table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1],
                             table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1],
                             table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1],
                             table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1],
                             table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1],
                             table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1],
                             table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1],
                             table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1],
                             table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1],
                             table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1],
                             table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1])

to_barplot$total_sig = c(table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2],
                             table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2],
                             table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2],
                             table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2],
                         table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2],
                         table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2],
                         table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2],
                         table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2],
                         table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2],
                         table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2],
                         table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2],
                         table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2])

to_barplot$celltype_sig = c(table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                            table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[2],
                            table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                            table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                            table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                            0,
                            table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                            table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                            table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                            0,
                            table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                            table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2])

to_barplot$celltype_non_sig = c(table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[1],
                            table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[1],
                            table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[1],
                            table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[1],
                            table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[1],
                            table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[1],
                            table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[1],
                            table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[1],
                            table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[1],
                            table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[1],
                            table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[1],
                            table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[1])

to_barplot$pvalue = c(table_res_fisher$Blue_pval, table_res_fisher$Red_pval, table_res_fisher$Grey_pval)
to_barplot$odds = c(table_res_fisher$Blue_odds, table_res_fisher$Red_odds, table_res_fisher$Grey_odds)
to_barplot$total_cpg = rowSums(to_barplot[,3:4])
to_barplot$total_subtype = rowSums(to_barplot[,5:6])

to_barplot$total_perc_sig = to_barplot$total_sig / to_barplot$total_cpg
to_barplot$total_perc_nonsig = to_barplot$total_non_sig / to_barplot$total_cpg
to_barplot$celltype_perc_sig = to_barplot$celltype_sig / to_barplot$total_subtype
to_barplot$celltype_perc_nonsig = to_barplot$celltype_non_sig / to_barplot$total_subtype

to_barplot_percs = to_barplot[,c(1,2,11,12,13,14)]

to_barplot_melt = melt(to_barplot_percs)
to_barplot_melt$cpg_level = str_split_fixed(as.character(to_barplot_melt$variable), "_",2)[,1]
to_barplot_melt$subtype = as.character(to_barplot_melt$subtype)
to_barplot_melt[to_barplot_melt$cpg_level == "celltype",]$cpg_level = to_barplot_melt[to_barplot_melt$cpg_level == "celltype",]$subtype

to_barplot_melt = to_barplot_melt[-(5:12),]

to_barplot_pvalue = to_barplot[,c(1,2,7)]
to_barplot_melt = merge(to_barplot_melt, to_barplot_pvalue)

to_barplot_melt[to_barplot_melt$variable %in% c("total_perc_sig","total_perc_nonsig","celltype_perc_nonsig"),]$pvalue = NA

to_barplot_melt$pvalue = as.numeric(format(to_barplot_melt$pvalue, digits = 3, scientific = TRUE))

to_barplot_melt$significant = factor(str_split_fixed(as.character(to_barplot_melt$variable), "_",3)[,3],
                                     levels = c("sig", "nonsig"))

to_barplot_melt$cell_type = factor(to_barplot_melt$cell_type, levels = c("NeuN","Olig","Mic","Ast"))

to_barplot_melt = to_barplot_melt[to_barplot_melt$significant == "sig",]
to_barplot_melt[to_barplot_melt$cpg_level == "total",]$cpg_level = "Complete set"
to_barplot_melt[to_barplot_melt$cpg_level == "red",]$cpg_level = "RED subtype"
to_barplot_melt[to_barplot_melt$cpg_level == "blue",]$cpg_level = "BLUE subtype"
to_barplot_melt[to_barplot_melt$cpg_level == "grey",]$cpg_level = "Unassigned AD"

to_barplot_melt$cpg_level = factor(to_barplot_melt$cpg_level, levels = c("Complete set", "BLUE subtype", "RED subtype", "Unassigned AD"))

to_barplot_melt[9,]$pvalue = NA
to_barplot_melt[13,]$pvalue = NA
to_barplot_melt[is.na(to_barplot_melt$pvalue),]$pvalue = 1

to_barplot_melt = to_barplot_melt[to_barplot_melt$cell_type != "Ast",]

to_barplot_melt$asterisks = NA

ggplot(to_barplot_melt, aes(x = cell_type, y = value, fill = cpg_level)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = scales::percent(value, 0.1)), position = position_dodge(width = 0.9), vjust = -0.5) +
  geom_text(aes(label = asterisks), position = position_dodge(width = 0.9), vjust = -0.5, size = 8) +
  # geom_text(aes(y = 0.15, label = scales::scientific(pvalue)), position = position_dodge(width = 0.9), vjust = -0.5) +
  # geom_text(aes(y = 0.16, label = if_else(cpg_level == "total", "Fisher's test P-value :", ""))) +
  geom_text(aes(y = -0.007, label =cpg_level ), position = position_dodge(width = 0.9), vjust = -0.5, angle = 45, size = 3) +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15)) + 
  scale_x_discrete(labels=c("NeuN" = "Neuronal", "Olig" = "Oligodendrocytes",
                            "Mic" = "Microglia","Ast" = "Astrocytes")) +
  xlab("Cell type") + ylab("Percentage of significance in the CpG set") +
  # geom_hline(yintercept=1.05) + geom_hline(yintercept=1.15) +
  scale_fill_discrete(type = c("black", "steelblue1", "brown2", "grey"), name = "CpG set") + 
  theme_bw() +
  theme(axis.text.x = element_text(size=10))
   
BG_Mic = All_CT[All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold,]
BG_Neun = All_CT[All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold,]
BG_Olig = All_CT[All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold,]
BG_Ast = All_CT[All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold,]
write.csv(BG_Mic, file ="SUBTYPING-PAPER-TABLEFIG/Sig_MIC_FANS.csv")
write.csv(BG_Neun, file ="SUBTYPING-PAPER-TABLEFIG/Sig_NEUN_FANS.csv")
write.csv(BG_Olig, file ="SUBTYPING-PAPER-TABLEFIG/Sig_OLIG_FANS.csv")
write.csv(BG_Ast, file ="SUBTYPING-PAPER-TABLEFIG/Sig_AST_FANS.csv")
## MIC
All_CT_red_mic = All_CT_red[All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold,]
All_CT_blue_mic = All_CT_blue[All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold,]
intersect(rownames(All_CT_blue_mic), rownames(All_CT_red_mic))
FANS_Gometh_Red = gometh(sig.cpg = rownames(All_CT_red_mic),
       all.cpg = rownames(BG_Mic$X),
       collection = c("GO", "KEGG"),
       array.type = c("450K", "EPIC"),
       sig.genes = FALSE)
FANS_Gometh_Red = FANS_Gometh_Red[FANS_Gometh_Red$P.DE < 0.05,]
FANS_Gometh_Red = FANS_Gometh_Red[order(FANS_Gometh_Red$P.DE),]
write.csv(FANS_Gometh_Red, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Red_MIC.csv")

FANS_Gometh_Blue = gometh(sig.cpg = rownames(All_CT_blue_mic),
       all.cpg = rownames(BG_Mic$X),
       collection = c("GO", "KEGG"),
       array.type = c("450K", "EPIC"),
       sig.genes = FALSE)
FANS_Gometh_Blue = FANS_Gometh_Blue[FANS_Gometh_Blue$P.DE < 0.05,]
FANS_Gometh_Blue = FANS_Gometh_Blue[order(FANS_Gometh_Blue$P.DE),]
write.csv(FANS_Gometh_Blue, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Blue_MIC.csv")

## NeuN RED
All_CT_red_neun = All_CT_red[All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold,]
FANS_Gometh_Red = gometh(sig.cpg = rownames(All_CT_red_neun),
                         all.cpg = rownames(BG_Neun$X),
                         collection = c("GO", "KEGG"),
                         array.type = c("450K", "EPIC"),
                         sig.genes = FALSE)
FANS_Gometh_Red = FANS_Gometh_Red[FANS_Gometh_Red$P.DE < 0.05,]
FANS_Gometh_Red = FANS_Gometh_Red[order(FANS_Gometh_Red$P.DE),]
write.csv(FANS_Gometh_Red, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Red_NEUN.csv")



## Oligo RED
All_CT_red_olig = All_CT_red[All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold,]
FANS_Gometh_Red = gometh(sig.cpg = rownames(All_CT_red_olig),
                         all.cpg = rownames(BG_Olig$X),
                         collection = c("GO", "KEGG"),
                         array.type = c("450K", "EPIC"),
                         sig.genes = FALSE)
FANS_Gometh_Red = FANS_Gometh_Red[FANS_Gometh_Red$P.DE < 0.05,]
FANS_Gometh_Red = FANS_Gometh_Red[order(FANS_Gometh_Red$P.DE),]
write.csv(FANS_Gometh_Red, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Red_OLIG.csv")




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
