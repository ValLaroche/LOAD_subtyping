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

add_subtype_pheno <- function(pheno, df_assignment){
  pheno$HC_clusters = NA
  pheno$KM_clusters = NA
  
  pheno[match(rownames(df_assignment), pheno$Basename),]$HC_clusters = df_assignment$HClust
  pheno[match(rownames(df_assignment), pheno$Basename),]$KM_clusters = df_assignment$Kmeans

  
  pheno[is.na(pheno$HC_clusters),]$HC_clusters = "Control"
  pheno[is.na(pheno$KM_clusters),]$KM_clusters = "Control"
  return(pheno)
}

calc_sva <- function(tmpcounts, tmppheno, cohort){
  
  if(cohort == "UKBBN"){
    tmppheno$Brain.Bank = factor(tmppheno$Brain.Bank, levels = unique(tmppheno$Brain.Bank))
    mod1 = model.matrix(~Gender + Age + Brain.Bank + RINscale + plate, data=tmppheno)  
  } else if (cohort == "PITT"){
    mod1 = model.matrix(~Sex + Age + RINscale + Plate, data=tmppheno)  
  } else if (cohort == "ROSMAP"){
    mod1 = model.matrix(~msex + age_death + RINscale + Sample_Plate, data=tmppheno)  
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

plot_convex_hull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

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
  
  #find sample size of each sample
  n1 <- length(cpg_group1)
  n2 <- length(cpg_group2)
  
  # sd_cpg = sd(c(cpg_group1, cpg_group2))
  
  poolsd_cpg = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n1-2))
  
  cohend = (mean_group1 - mean_group2) / poolsd_cpg
  return(cohend)
}

make_RF <- function(training){
  set.seed(4)
  
  fitControl <- trainControl(method = "repeatedcv", 
                             number = 10, 
                             repeats = 10, 
                             classProb=TRUE,  
                             savePredictions = TRUE, 
                             summaryFunction = twoClassSummary)
  
  rf_OMICS = train(subtype ~ ., data = training,
                   method = "rf",
                   trControl = fitControl,
                   verbose = FALSE,
                   tuneLength = 20)
  
  return(rf_OMICS)
}

make_tableres <- function(model, training, testing){
  tableres_test = data.frame("RID" = rownames(testing),
                             "Target_test" = testing$subtype,
                             "Predicted_test" = predict(model, newdata = testing))
  tableres_train = data.frame("RID" = rownames(training),
                              "Target_train" = training$subtype,
                              "Predicted_train" = predict(model, newdata = training))
  
  
  tableres_test = cbind(tableres_test, as.data.frame(predict(model, newdata = testing, type = "prob"))[,,1])
  tableres_train = cbind(tableres_train, as.data.frame(predict(model, newdata = training, type = "prob"))[,,1])
  
  colnames(tableres_test)[4] = unique(testing$subtype)[1]
  colnames(tableres_test)[5] = unique(testing$subtype)[2]
  colnames(tableres_train)[4] = unique(training$subtype)[1]
  colnames(tableres_train)[5] = unique(training$subtype)[2]
  
  return(list(tableres_test, tableres_train))
}

make_ROC <- function(tableres, roc_data, MCC_data, outpath){
  pdf(file = outpath, width = 7, height = 7)
  plot.roc(roc_data)
  legend("topleft", legend = paste("MCC =", MCC_data[16]))
  dev.off()
}

make_CM <- function(tableres, outpath){
  conf_mat = confusion_matrix(targets = tableres$Target,
                              predictions = tableres$Predicted)
  
  plot_cm = plot_confusion_matrix(
    conf_mat$`Confusion Matrix`[[1]],
    add_sums = TRUE,
    sums_settings = sum_tile_settings(
      palette = "Oranges",
      label = "Total",
      tc_tile_border_color = "black",
    ),
    font_row_percentages = font(size = 5),
    font_col_percentages = font(size = 5),
    font_counts = font(size = 7))
  plot_cm = plot_cm + theme(axis.title = element_text(size = 18),
                            axis.text = element_text(size = 14))   
  
  pdf(file = outpath, width = 7, height = 7)
  plot(plot_cm)
  dev.off()
  
}

make_impt <- function(model, num_features, outpath, model_type = "nonplsda"){
  if(model_type == "nonplsda"){
    
    top_features = min(50, num_features)
    plot_impt = plot(varImp(model, scale = F), top = top_features)
    
    pdf(file = outpath, width = 7, height = 7)
    plot(plot_impt)
    dev.off()
  } else {
    var_imp = data.frame(feature = rownames(model$betahat),
                         value = model$betahat[,1])
    var_imp = var_imp[order(abs(var_imp$value), decreasing = T),]
    var_imp$feature = factor(var_imp$feature, levels = var_imp$feature)
    
    top_features = min(50, num_features)
    var_imp = var_imp[seq(1:top_features),]
    
    plot_impt = ggplot(var_imp, aes(x = rev(feature), y = value)) +
      ggtitle("Overall variable importance of SPLSDA model") +
      geom_bar(stat = "identity", width = 0.1) +
      geom_point(colour = "#529EFF") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      xlab("Features") +
      coord_flip()
    
    pdf(file = outpath, width = 7, height = 7)
    plot(plot_impt)
    dev.off()
  }
}

make_opti <- function(df_AUC_MCC, outpath){
  df_AUC_MCC_melt = melt(df_AUC_MCC, id.vars = c("model", "nfeatures", "train_test"), measure.vars = c("AUC", "MCC"))
  df_AUC_MCC_melt$value = as.numeric(df_AUC_MCC_melt$value)
  df_AUC_MCC_melt$nfeatures = as.numeric(df_AUC_MCC_melt$nfeatures)
  AUCMCC_plot = ggplot(df_AUC_MCC_melt, aes(x = nfeatures, y = value, group = variable, color = variable)) +
    geom_point() +
    ylim(-1, 1) +
    facet_wrap(~train_test, nrow = 2) +
    theme_bw()
  
  
  pdf(file = outpath, width = 8, height = 7)
  plot(AUCMCC_plot)
  dev.off()
}