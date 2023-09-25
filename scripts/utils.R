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

make_dataset <- function(data_matrix, metadata, dirname, db, cpgs_intersect){
  if(db == "BDR"){
    metadata_AD = metadata[metadata$diag == "AD",]
    metadata_AD = metadata_AD[order(metadata_AD$Basename),]
    
    data_matrix_AD = data_matrix[,colnames(data_matrix) %in% metadata_AD$Basename]
    data_matrix_AD = data_matrix_AD[rownames(data_matrix_AD) %in% cpgs_intersect,]
    
  } else if(db == "UKBBN"){
    metadata_AD = metadata[metadata$AD == 1,]
    metadata_AD = metadata_AD[order(metadata_AD$Row.names),]
    
    data_matrix_AD = data_matrix[,colnames(data_matrix) %in% metadata_AD$Row.names]
    data_matrix_AD = data_matrix_AD[rownames(data_matrix_AD) %in% cpgs_intersect,]
    
  } else if(db == "PITT"){
    metadata_AD = metadata[metadata$AD == 1,]
    metadata_AD = metadata_AD[order(metadata_AD$Basename),]
    
    data_matrix_AD = data_matrix[,colnames(data_matrix) %in% metadata_AD$Basename]
    data_matrix_AD = data_matrix_AD[rownames(data_matrix_AD) %in% cpgs_intersect,]
    
  }# Add PITT here why is it missing ? ? ?
  
  
  data_matrix_AD = data_matrix_AD[,order(colnames(data_matrix_AD))]
  
  dir.create("results/")
  dir.create(paste0("results/", dirname))
  dir.create(paste0("results/", dirname, "/Spectrum"))
  dir.create(paste0("results/", dirname, "/HClust"))
  dir.create(paste0("results/", dirname, "/Kmeans"))
  dir.create(paste0("results/", dirname, "/PAM"))
  
  output = list(data_matrix_AD, metadata_AD)
  return(output)
}

wrapper_models <- function(data_matrix_AD, metadata_AD, cluster_vector = c(), dirname){
  
  list_datasets = list(data_matrix_AD) #Input of models is a list (Because of multiomics)
  
  if(length(cluster_vector) == 0){ # Base number of clusters is 3
    cluster_vector = c(3,3,3,3)
  }
  #Build spectrum clustering and eigenvalue check for best number of clusters
  multi_Spectrum = Spectrum_modelization(list_datasets, spectrum_method = 2, 
                                         cluster_alg = "GMM", num_clusters = 2)
  print(multi_Spectrum)
  png(paste0("results/", dirname, "/Spectrum/eigenvector_analysis.png"), width = 800, height = 600)
  plot(multi_Spectrum$eigenvector_analysis$K, multi_Spectrum$eigensystem$values[1:11])
  dev.off()
  #Build optimal spectrum
  multi_Spectrum = Spectrum_modelization(list_datasets, spectrum_method = 3, 
                                         cluster_alg = "GMM", num_clusters = cluster_vector[1])
  print(multi_Spectrum)
  #Plot results
  plot_model(multi_Spectrum, "Spectrum", colnames(data_matrix_AD), paste0("results/", dirname, "/Spectrum/"))
  df_assignment = data.frame("Spectrum" = multi_Spectrum$assignments,
                             row.names = names(multi_Spectrum$assignments))  
  #Build Hclust clustering (CC generates best number of cluster figures auto)
  multi_HC = CC_modelization(list_datasets[[1]],
                             title_model = paste0("results/", dirname, "/HClust/"), 
                             cluster_alg = "hc", distance = "euclidean")
  plot_model(multi_HC, "CC", colnames(data_matrix_AD), paste0("results/", dirname, "/HClust/"),
             numCC_clusters = cluster_vector[2])
  df_assignment$HC = multi_HC[[cluster_vector[2]]]$consensusClass
  
  #Build Kmeans clustering (CC generates best number of cluster figures auto)
  multi_KM = CC_modelization(list_datasets[[1]],
                             title_model = paste0("results/", dirname, "/Kmeans/"), 
                             cluster_alg = "km", distance = "euclidean")
  plot_model(multi_KM, "CC", colnames(data_matrix_AD), paste0("results/", dirname, "/Kmeans/"),
             numCC_clusters = cluster_vector[3])
  df_assignment$KM = multi_KM[[cluster_vector[3]]]$consensusClass
  
  #Build PAM clustering (CC generates best number of cluster figures auto)
  multi_PAM = CC_modelization(list_datasets[[1]],
                              title_model = paste0("results/", dirname, "/PAM/"), 
                              cluster_alg = "pam", distance = "euclidean")
  plot_model(multi_PAM, "CC", colnames(data_matrix_AD), paste0("results/", dirname, "/PAM/"),
             numCC_clusters = cluster_vector[4])
  df_assignment$PAM = multi_PAM[[cluster_vector[4]]]$consensusClass
  
  colnames(df_assignment) = c("Spectrum", "HClust", "Kmeans", "PAM")
  df_NMI_models = calcNMI_all(df_assignment)
  models_tables(df_assignment, df_NMI_models, paste0("results/", dirname, "/"))
  
  return(df_assignment)
  
}

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

resid = function(row, age, sex, prop){
  fit = try(
    lm( row ~ age + sex + prop),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals

  dat.reg = {
    sex	= 	as.factor(pheno$Gender)
    age	= 	pheno$Ageatdeath
    prop =  as.numeric(pheno$prop)
    t(apply(mydat, 1, resid, age, sex, prop))
  }
}

CDF <- function(CDF_model, breaks = 100){
  plot(c(0), xlim = c(0, 1), ylim = c(0, 1), col = "white", 
       bg = "white", xlab = "consensus index", ylab = "CDF", 
       main = "consensus CDF", las = 2)
  k = length(CDF_model)
  this_colors = rainbow(k - 1)
  areaK = c()
  for (i in 2:length(CDF_model)) {
    v = triangle(CDF_model[[i]], mode = 1)
    h = hist(v, plot = FALSE, breaks = seq(0, 2, by = 1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)
    thisArea = 0
    for (bi in 1:(length(h$breaks) - 1)) {
      thisArea = thisArea + h$counts[bi] * (h$breaks[bi + 
                                                       1] - h$breaks[bi])
      bi = bi + 1
    }
    areaK = c(areaK, thisArea)
    lines(h$mids, h$counts, col = this_colors[i - 1], lwd = 2, 
          type = "l")
  }
  legend(0.8, 0.5, legend = paste(rep("", k - 1), seq(2, k, 
                                                      by = 1), sep = ""), fill = this_colors)
  deltaK = areaK[1]
  for (i in 2:(length(areaK))) {
    deltaK = c(deltaK, (areaK[i] - areaK[i - 1])/areaK[i - 
                                                         1])
  }
  plot(1 + (1:length(deltaK)), y = deltaK, xlab = "k", ylab = "relative change in area under CDF curve", 
       main = "Delta area", type = "b")
}


triangle <- function(CDF_model, mode = 1){
  n = dim(CDF_model)[1]
  nm = matrix(0, ncol = n, nrow = n)
  fm = CDF_model
  nm[upper.tri(nm)] = CDF_model[upper.tri(CDF_model)]
  fm = t(nm) + nm
  diag(fm) = diag(CDF_model)
  nm = fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = CDF_model[lower.tri(nm)]
  if (mode == 1) {
    return(vm)
  }
  else if (mode == 3) {
    return(fm)
  }
  else if (mode == 2) {
    return(nm)
  }
}

calc_sva <- function(tmpcounts, tmppheno, cohort){
  
  if(cohort == "UKBBN"){
    mod1 = model.matrix(~Gender + Age + plate + Brain.Bank, data=tmppheno)  
  } else if (cohort == "PITT"){
    mod1 = model.matrix(~Sex.y + Age + Plate, data=tmppheno)  
  }
  mod0 = model.matrix(~1, data=tmppheno)
  set.seed(1234)
  
  
  tmpcounts = as.matrix(tmpcounts)
  tmpcounts = na.omit(tmpcounts)
  
  count_mean = apply(tmpcounts,1,mean)
  filter_var0 = names(subset(count_mean ,count_mean==0))
  tmpcounts = tmpcounts[-which(rownames(tmpcounts) %in% filter_var0),]
  
  svseq = sva::svaseq(tmpcounts,mod1,mod0)$sv
  svseq = data.frame(svseq)
  colnames(svseq) = paste0("sv",seq(1,ncol(svseq)))
  tmppheno = cbind(tmppheno, svseq)
  
  return(tmppheno)
}


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
