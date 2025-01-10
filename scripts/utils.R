# Load all relevant libraries to the project
library(SNFtool)
library(Spectrum)
library(cluster)
library(ConsensusClusterPlus)
library(CancerSubtypes)
library(iClusterPlus)
library(ANF)
library(ggplot2)
library(factoextra)
library(grid)
library(gridExtra)
library(gtable) 
library(gplots)
library(r.jive)
library(mixOmics)
library(iClusterPlus)
library(reshape2)
library(WGCNA)
library(corrplot)
library(sva)
library(splitstackshape)
library(parallel)
library(snow)
library("pbapply")
library(aricode)
library(prospectr)
library(caret)
library(dplyr)
library(pheatmap)
library(biomaRt)
library(lsr)
library(foreach)
library(doParallel)
library(GeneOverlap)
library(edgeR)
library(ggrepel)
library(ggpubr)
library(coloc)
library(moments)
library(gap)
library(bacon)
library(limma)
library(ggvenn)
library(rstatix)
library(WGCNA)
library(sva)
library(stringr)
library("QuantPsyc")
library(grDevices)
library(pracma)
library(metafor)
library(ggfortify)  
library(missMethyl)
library(pbapply)
library(RRtest)
library(caret)
library(cvms)
library(pROC)
library(VennDiagram)

# Generate standard names for clustering pipeline
make_dataset <- function(data_matrix, metadata, dirname, db, cpgs_intersect){
    if(db == "UKBBN"){
    metadata_AD = metadata[metadata$AD == 1,]
    metadata_AD = metadata_AD[order(metadata_AD$Row.names),]
    
    data_matrix_AD = data_matrix[,colnames(data_matrix) %in% metadata_AD$Row.names]
    data_matrix_AD = data_matrix_AD[rownames(data_matrix_AD) %in% cpgs_intersect,]
    
  } else if(db == "PITT"){
    metadata_AD = metadata[metadata$Phenotype %in% c("AD-P", "AD+P"),]
    metadata_AD = metadata_AD[order(metadata_AD$Basename),]
    
    data_matrix_AD = data_matrix[,colnames(data_matrix) %in% metadata_AD$Basename]
    data_matrix_AD = data_matrix_AD[rownames(data_matrix_AD) %in% cpgs_intersect,]
    
  } else if(db == "ROSMAP"){
    metadata_AD = metadata[metadata$diag == "AD",]
    metadata_AD = metadata_AD[order(metadata_AD$Full_Sentrix),]
    
    data_matrix_AD = data_matrix[,colnames(data_matrix) %in% metadata_AD$Full_Sentrix]
    data_matrix_AD = data_matrix_AD[rownames(data_matrix_AD) %in% cpgs_intersect,]
    
  }
  
  
  data_matrix_AD = data_matrix_AD[,order(colnames(data_matrix_AD))]
  
  dir.create("results/")
  dir.create(paste0("results/", dirname))
  dir.create(paste0("results/", dirname, "/HClust"))
  dir.create(paste0("results/", dirname, "/Kmeans"))
  
  output = list(data_matrix_AD, metadata_AD)
  return(output)
}

# Make clustering models
wrapper_models <- function(data_matrix_AD, metadata_AD, cluster_vector = c(), dirname){
  
  list_datasets = list(data_matrix_AD) #Input of models is a list (Because of multiomics)
  
  if(length(cluster_vector) == 0){ # Base number of clusters is 3
    cluster_vector = c(3,3)
  }
  df_assignment = data.frame(matrix(nrow = ncol(data_matrix_AD), ncol = 2))
  colnames(df_assignment) = c("HC", "KM")
 #Build Hclust clustering (CC generates best number of cluster figures auto)
  multi_HC = CC_modelization(list_datasets[[1]],
                             title_model = paste0("results/", dirname, "/HClust/"), 
                             cluster_alg = "hc", distance = "euclidean")
  plot_model(multi_HC, "CC", colnames(data_matrix_AD), paste0("results/", dirname, "/HClust/"),
             numCC_clusters = cluster_vector[1])
  df_assignment$HC = multi_HC[[cluster_vector[1]]]$consensusClass
  
  #Build Kmeans clustering (CC generates best number of cluster figures auto)
  multi_KM = CC_modelization(list_datasets[[1]],
                             title_model = paste0("results/", dirname, "/Kmeans/"), 
                             cluster_alg = "km", distance = "euclidean")
  plot_model(multi_KM, "CC", colnames(data_matrix_AD), paste0("results/", dirname, "/Kmeans/"),
             numCC_clusters = cluster_vector[2])
  df_assignment$KM = multi_KM[[cluster_vector[2]]]$consensusClass
  
  colnames(df_assignment) = c("HClust", "Kmeans")
  df_NMI_models = calcNMI_all(df_assignment)
  models_tables(df_assignment, df_NMI_models, paste0("results/", dirname, "/"))
  
  return(df_assignment)
  
}

# Generate heatmap of the clustering models results (unused)
model_heatmap = function (model_matrix, model_clusters, label = NULL){
  sorted = sort(table(model_clusters))
  o.stol = as.numeric(names(sorted))
  o = NULL
  for (i in o.stol) {
    o = c(o, which(model_clusters == i))
  }

  label = label[o]
  model_matrix = model_matrix[o,o]
  rownames(model_matrix) = rev(label)
  colnames(model_matrix) = label
  
  factor(label, levels = label)
  melt_matrix = melt(model_matrix)
  melt_matrix$value <- (melt_matrix$value - min(melt_matrix$value)) / (max(melt_matrix$value) - min(melt_matrix$value))  
  plot_clusters = ggplot(data = melt_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(plot_clusters)
}

# Add clustering results (HC/KM) labels to the complet pheno
add_subtype_pheno <- function(pheno, df_assignment){
  pheno$HC_clusters = NA
  pheno$KM_clusters = NA
  
  pheno[match(rownames(df_assignment), pheno$Basename),]$HC_clusters = df_assignment$HClust
  pheno[match(rownames(df_assignment), pheno$Basename),]$KM_clusters = df_assignment$Kmeans

  
  pheno[is.na(pheno$HC_clusters),]$HC_clusters = "Control"
  pheno[is.na(pheno$KM_clusters),]$KM_clusters = "Control"
  return(pheno)
}

# Generate EWASs based on SVs as covariates and returning the results (per CPG)
EWAS <- function(x, group,tmppheno){
  sv1 = tmppheno$sv1
  sv2 = tmppheno$sv2
  sv3 = tmppheno$sv3
  # sv4 = tmppheno$sv4
  # sv5 = tmppheno$sv5
  # sv6 = tmppheno$sv6
  # sv7 = tmppheno$sv7
  # sv8 = tmppheno$sv8
  # sv9 = tmppheno$sv9
  # sv10 = tmppheno$sv10
  xx = as.numeric(x)
  fit<-try(lm(group ~ xx + sv1 + sv2 + sv3,na.action=na.omit))
  if(inherits(fit,'try-error')) return(rep(NA,4))
  return(c(coef(summary(fit))[2,],as.numeric(QuantPsyc::lm.beta(fit)[1]), coef(summary(fit))[2,2] /sd(group)))
}

# Generate EWASs for SNP results (unused)
EWAS_SNP <- function(x, group,tmppheno){
  PC1 = tmppheno$PC1
  PC2 = tmppheno$PC2
  PC3 = tmppheno$PC3
  PC4 = tmppheno$PC4
  PC5 = tmppheno$PC5
  PC6 = tmppheno$PC6
  PC7 = tmppheno$PC7
  PC8 = tmppheno$PC8
  PC9 = tmppheno$PC9
  PC10 = tmppheno$PC10
  
  xx = as.numeric(x)
  fit<-try(lm(group ~ xx + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,na.action=na.omit))
  if(inherits(fit,'try-error')) return(rep(NA,4))
  return(c(coef(summary(fit))[2,],as.numeric(QuantPsyc::lm.beta(fit)[1]), coef(summary(fit))[2,2] /sd(group)))
}

# Calculate surrogate variables for the methylation level
calc_sva_methyl <- function(tmpcounts, tmppheno, cohort){
  
  if(cohort == "UKBBN"){
    tmppheno$Brain.Bank = factor(tmppheno$Brain.Bank, levels = unique(tmppheno$Brain.Bank))
    mod1 = model.matrix(~Sex + Age + Brain.Bank +
                        Astrocyte + Microglia + Oligodendrocyte + Neuronal_Inhibitory + Neuronal_Excitatory, data=tmppheno)  
  } else if (cohort == "BDR"){
    mod1 = model.matrix(~Gender + Age + Plate +
                        Astrocyte + Microglia + Oligodendrocyte + Neuronal_Inhibitory + Neuronal_Excitatory, data=tmppheno)  
  } else if (cohort == "PITT"){
    mod1 = model.matrix(~Sex + Age + Plate +
                        Astrocyte + Microglia + Oligodendrocyte + Neuronal_Inhibitory + Neuronal_Excitatory, data=tmppheno)  
  } else if (cohort == "ROSMAP"){
    mod1 = model.matrix(~msex + age_death + Sample_Plate +
                        Astrocyte + Microglia + Oligodendrocyte + Neuronal_Inhibitory + Neuronal_Excitatory, data=tmppheno)  
  }
  mod0 = model.matrix(~1,data=tmppheno)
  set.seed(1234)
  
  tmpcounts = na.omit(tmpcounts)
  tmpcounts = as.matrix(tmpcounts)
  
  svseq = sva::sva(tmpcounts,mod1,mod0, n.sv = 10)$sv
  svseq = data.frame(svseq)
  colnames(svseq) = paste0("sv",seq(1,ncol(svseq)))
  tmppheno = cbind(tmppheno, svseq)
  
  return(tmppheno)
}

#Calculate surrogate variables for the FANS level
calc_sva_FANS <- function(tmpcounts, tmppheno){
  
  tmppheno$Brain.Bank = factor(tmppheno$Brain.Bank, levels = unique(tmppheno$Brain.Bank))
  mod1 = model.matrix(~Sex + Age + Brain.Bank, data=tmppheno)  
  mod0 = model.matrix(~1,data=tmppheno)
  set.seed(1234)
  
  tmpcounts = na.omit(tmpcounts)
  tmpcounts = as.matrix(tmpcounts)
  
  svseq = sva::sva(tmpcounts,mod1,mod0, n.sv = 10)$sv
  svseq = data.frame(svseq)
  colnames(svseq) = paste0("sv",seq(1,ncol(svseq)))
  tmppheno = cbind(tmppheno, svseq)
  
  return(tmppheno)
}

#Calculate surrogate variables for the FANS level
calc_sva <- function(tmpcounts, tmppheno, cohort){
  
  if(cohort == "UKBBN"){
    tmppheno$Brain.Bank = factor(tmppheno$Brain.Bank, levels = unique(tmppheno$Brain.Bank))
    mod1 = model.matrix(~Gender + Age + Brain.Bank + RIN + plate + well, data=tmppheno)  
  } else if (cohort == "PITT"){
    mod1 = model.matrix(~Sex + Age + pool + plate + RIN, data=tmppheno)  
  } else if (cohort == "ROSMAP"){
    mod1 = model.matrix(~msex + age_death + seqbatch + RIN, data=tmppheno)  
  }
  mod0 = model.matrix(~1,data=tmppheno)
  set.seed(1234)
  
  tmpcounts = na.omit(tmpcounts)
  tmpcounts = as.matrix(tmpcounts)
  
  svseq = sva::sva(tmpcounts,mod1,mod0, n.sv = 10)$sv
  svseq = data.frame(svseq)
  colnames(svseq) = paste0("sv",seq(1,ncol(svseq)))
  tmppheno = cbind(tmppheno, svseq)
  
  return(tmppheno)
}

# Associate biomart gene IDs (unused)
biomart_gene_assoc <- function(res_table, biomart_names, Pval_threshold, FC_threshold){
  
  res_table = res_table[-log10(res_table$PValue) > Pval_threshold & abs(res_table$logFC) > FC_threshold,]
  
  setdiff(res_table$ENSG,biomart_names$Gene_stable_ID_version) # Release 104 vs 109
  
  biomart_names = biomart_names[order(biomart_names$Gene_stable_ID_version),]
  res_table = res_table[order(res_table$ENSG),]
  
  unique(biomart_names$Gene_name)
  biomart_names[sort(match(res_table$ENSG, biomart_names$Gene_stable_ID_version)),]
  res_table$genenames = NA
  res_table$genenames = biomart_names[match(res_table$ENSG, biomart_names$Gene_stable_ID_version),]$Gene_name
  
  res_table[!res_table$ENSG %in% biomart_names$Gene_stable_ID_version,]
  
  return(res_table)
}

# Plotting the convex hulls when establishing the subtypes
plot_convex_hull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

# Calculating effect sizes for the methylation data
cohenD <- function(cpg, pheno, groups, AD = FALSE){
  if(AD){
    pheno_1 = pheno[pheno$diag %in% groups[1],]
    pheno_2 = pheno[pheno$diag %in% groups[2],]
  } else {
    pheno_1 = pheno[pheno$subtype %in% groups[1],]
    pheno_2 = pheno[pheno$subtype %in% groups[2],]
    
  }
  cpg_group1 = as.numeric(cpg[names(cpg) %in% pheno_1$Basename])
  cpg_group2 = as.numeric(cpg[names(cpg) %in% pheno_2$Basename])
  
  mean_group1 = mean(cpg_group1)
  mean_group2 = mean(cpg_group2)
  
  s1 <- sd(cpg_group1)
  s2 <- sd(cpg_group2)
  
  n1 <- length(cpg_group1)
  n2 <- length(cpg_group2)
  
  poolsd_cpg = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  cohend = (mean_group1 - mean_group2) / poolsd_cpg
  return(cohend)
}