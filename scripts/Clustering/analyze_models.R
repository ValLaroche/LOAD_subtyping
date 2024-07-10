# Plot output graphs
plot_model <- function(model, model_type, label_samples, outdir, numCC_clusters = 0){
  
  model_clusters = model[[numCC_clusters]]$consensusClass
  model_matrix = model[[numCC_clusters]]$consensusMatrix
  
  sil_model = silhouette(model_clusters, model_matrix)
  pdf(paste0(outdir, "silhouette_plot.pdf"), width = 7, height = 7)
  plot(sil_model,col = sil_model[,1][order(sil_model[,1])])
  dev.off()
  
  plot_clusters = model_heatmap(model_matrix = model_matrix,
                                model_clusters = model_clusters,
                                label = label_samples)
  
  pdf(paste0(outdir, "clusters_heatmap.pdf"), width = 7, height = 7)
  plot(plot_clusters)
  dev.off()
  
}

# Plot NMI table and 2x2 tables for cluster comparison
models_tables <- function(df_assignment, df_NMI_models, outdir){
  tables_res = list()
  
  for(i in c(1:length(df_assignment))){
    for(j in c(2:length(df_assignment))){
      if(i != j && i < j){
        
        table_res = tableGrob(table(df_assignment[[i]], df_assignment[[j]]), 
                              theme = ttheme_minimal(base_size = 7, padding = unit(c(4,2),"mm"),
                                rowhead=list(fg_params=list(col="black", fontface=2L,
                                                            hjust = 1.5))))
        
        table_res = gtable_add_grob(table_res,
                             grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                             t = 2, b = nrow(table_res), l = 1, r = ncol(table_res))
        table_res = gtable_add_grob(table_res,
                             grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                             t = 1, b = nrow(table_res), l = 2, r = ncol(table_res))
        
        h = grobHeight(table_res)
        w = grobWidth(table_res)
        
        
        table1_label = textGrob(names(df_assignment)[i], x=unit(0.5, "npc") - w*1, 
                                y=unit(0.5,"npc"), vjust=1, gp=gpar(fontsize=8, fontface = "bold"), rot = 90)
        table2_label = textGrob(names(df_assignment)[j], x=unit(0.5, "npc") + w*0.1, 
                                y=unit(0.5,"npc") + 1.05*h, vjust=0, gp=gpar(fontsize=8, fontface = "bold"))
        
        final_res = gTree(children = gList(table_res, table1_label, table2_label))
        tables_res[[length(tables_res)+1]] = final_res
      }
    }
  }
  
  plot_clusters = ggplot(data = df_NMI_models, aes(x = model_1, y = model_2, fill = nmi)) +
    ggtitle("NMI between models") +
    geom_tile() +
    geom_text(aes(label = round(nmi, 2), size = 13)) +
    scale_y_discrete(limits = rev(levels(factor(colnames(df_assignment), levels = colnames(df_assignment))))) +
    scale_x_discrete(limits = levels(factor(colnames(df_assignment), levels = colnames(df_assignment)))) +
    scale_fill_gradient2(low="#ffffb2", mid="#fd8d3c", high="firebrick", #colors in the scale
                         midpoint = 0.5, limits=c(0, 1)) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          panel.border = element_blank(),
          legend.position = "none")
  
  tables_res[[length(tables_res)+1]] = plot_clusters
  
  if(length(tables_res) == 7) {
    final_plot = arrangeGrob(
      grobs = tables_res,
      layout_matrix = rbind(c(7,7,7,1,4),
                            c(7,7,7,2,5),
                            c(7,7,7,3,6)),
      top = textGrob("Comparative results of clustering methods", gp = gpar(fontsize = 14))
    )
    ggsave(paste0(outdir, "Comparative_clustering.pdf"), final_plot, units = "px", width  = 3000, height = 1500)
  } else if(length(tables_res) == 11){
    final_plot = arrangeGrob(
      grobs = tables_res,
      layout_matrix = rbind(c(11,11,11,11,1,2),
                            c(11,11,11,11,3,4),
                            c(11,11,11,11,5,6),
                            c(11,11,11,11,7,8),
                            c(11,11,11,11,9,10)),
      
      top = textGrob("Comparative results of clustering methods", gp = gpar(fontsize = 16))
    )
    ggsave(paste0(outdir, "Comparative_clustering.pdf"), final_plot, units = "px", width  = 3000, height = 1500)
  }
}

# Compute NMI for a given clustering result
calcNMI_all <- function(df_assignment){
  df_NMI_models = data.frame(matrix(ncol = 3))
  colnames(df_NMI_models) = c("model_1","model_2","nmi")
  
  if(ncol(df_assignment >= 2)){
    row_index = 1
    for(i in c(1:ncol(df_assignment))){
      for(j in c(i:ncol(df_assignment))){  
        df_NMI_models[row_index,] = c(colnames(df_assignment)[i],
                               colnames(df_assignment)[j],
                               calNMI(df_assignment[,i], df_assignment[,j]))
        row_index = row_index + 1
      }
    }
  }
  
  df_NMI_models$nmi = as.numeric(df_NMI_models$nmi)
  return(df_NMI_models)
}

# Generate all plots relevant for DIABLO results
plot_mixOmics_all <- function(model, outdir){
  
  pdf(paste0(outdir, "diablo_indiv_plotbackground12.pdf"), width = 7, height = 7)
  plotIndiv(model,ind.names = FALSE, comp = c(1,2), 
            legend=TRUE, ellipse = TRUE,
            title = 'Methyl') 
  dev.off()
  
  pdf(paste0(outdir, "diablo_indiv_plotbackground13.pdf"), width = 7, height = 7)
  plotIndiv(model,ind.names = FALSE, comp = c(1,3), 
            legend=TRUE, ellipse = TRUE,
            title = 'Methyl')
  dev.off()
  
  pdf(paste0(outdir, "diablo_indiv_plotbackground14.pdf"), width = 7, height = 7)
  plotIndiv(model,ind.names = FALSE, comp = c(1,4), 
            legend=TRUE, ellipse = TRUE,
            title = 'Methyl') 
  dev.off()
  
  pdf(paste0(outdir, "diablo_var_plot.pdf"), width = 7, height = 7)
  plotVar(model)
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp1.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 1, contrib = "max")
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp2.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 2, contrib = "max")
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp3.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 3, contrib = "max")
  dev.off()
  
  pdf(paste0(outdir, "diablo_loading_comp4.pdf"), width = 7, height = 7)
  plotLoadings(model, comp = 4, contrib = "max")
  dev.off()

  rocroc = auroc(model, roc.comp = 6)
  pdf(paste0(outdir, "diablo_roc_curves.pdf"), width = 7, height = 7)
  plot(rocroc$graph.Comp6 + theme_minimal())
  dev.off()
  
  
  colorscim = c("#388ECC","#F68B33","#C2C2C2","#009E73","#CC79A7","#F0E442")
  length(colorscim) = length(levels(model$Y))
  # set the styling of the legend to be homogeneous with previous plots
  legend=list(legend = levels(model$Y), # set of classes
              col = colorscim, # set of colours
              title = "Predicted subtypes", # legend title
              cex = 0.7) # legend size
  pdf(paste0(outdir, "diablo_heatmap.pdf"), width = 20, height = 17)
  # generate the CIM, using the legend and colouring rows by each sample's class
  
  cim(model, row.sideColors = color.mixo(model$Y), 
             legend = legend, mapping = "X")
  
  dev.off()
}

# Loop for bootstrapping random permutations
random_permutations <- function(data_matrix_AD, metadata_AD, diablo_Methyls){
  #Get original 1200 features stats
  HC_roc = auroc(diablo_Methyls[[1]], plot = FALSE, print = FALSE)
  KM_roc = auroc(diablo_Methyls[[2]], plot = FALSE, print = FALSE)
  PAM_roc = auroc(diablo_Methyls[[3]], plot = FALSE, print = FALSE)
  Spec_roc = auroc(diablo_Methyls[[4]], plot = FALSE, print = FALSE)
  final_aucs = c(mean(HC_roc$Comp6[,1]),
                 mean(KM_roc$Comp6[,1]),
                 mean(PAM_roc$Comp6[,1]),
                 mean(Spec_roc$Comp6[,1]))
  
  #Get top features
  HC_var = plotVar(diablo_Methyls[[1]], comp.select = c(1:6), plot = FALSE)
  HC_features = HC_var$names
  
  KM_var = plotVar(diablo_Methyls[[2]], comp.select = c(1:6), plot = FALSE)
  KM_features = KM_var$names
  
  PAM_var = plotVar(diablo_Methyls[[3]], comp.select = c(1:6), plot = FALSE)
  PAM_features = PAM_var$names
  
  Spec_var = plotVar(diablo_Methyls[[4]], comp.select = c(1:6), plot = FALSE)
  Spec_features = Spec_var$names
  
  df_shuffle = data.frame(matrix(nrow = 1, ncol = 12))
  colnames(df_shuffle) = c("HC_auc","HC_ratio", "HC_trueprop",
                           "KM_auc","KM_ratio", "KM_trueprop",
                           "PAM_auc","PAM_ratio", "PAM_trueprop",
                           "Spec_auc","Spec_ratio", "Spec_trueprop")
  #Generate random permutations of labeling from 100% random to 1% random
for(i in c(rep(seq(from = 1, to = 0.1, by = -0.1), each = 10), 
                             rep(seq(from = 0.09, to = 0.01, by = -0.01),each = 10))){
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$HC_clusters)
      HC_shuffle = metadata_AD$HC_clusters
      HC_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$HC_clusters, HC_shuffle)
      HC_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$KM_clusters)
      KM_shuffle = metadata_AD$KM_clusters
      KM_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$KM_clusters, KM_shuffle)
      KM_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$PAM_clusters)
      PAM_shuffle = metadata_AD$PAM_clusters
      PAM_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$PAM_clusters, PAM_shuffle)
      PAM_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Randomise outcomes
      to_random = sample(seq(1:nrow(metadata_AD)), size = i*nrow(metadata_AD))
      shuffle_clust = sample(metadata_AD[to_random,]$Spec_clusters)
      Spec_shuffle = metadata_AD$Spec_clusters
      Spec_shuffle[to_random] = shuffle_clust
      df_random = data.frame(metadata_AD$Spec_clusters, Spec_shuffle)
      Spec_prop = 1 - length(which(df_random[,1] == df_random[,2])) / nrow(metadata_AD)
      #Make new models with random outcome
      diablo_shuffles = list()
      diablo_shuffles[[1]] = mixOmics::splsda(t(data_matrix_AD), HC_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      diablo_shuffles[[2]] = mixOmics::splsda(t(data_matrix_AD), KM_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      diablo_shuffles[[3]] = mixOmics::splsda(t(data_matrix_AD), PAM_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      diablo_shuffles[[4]] = mixOmics::splsda(t(data_matrix_AD), Spec_shuffle, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
      #Calculate stats
      HC_roc = auroc(diablo_shuffles[[1]], plot = FALSE, print = FALSE)
      KM_roc = auroc(diablo_shuffles[[2]], plot = FALSE, print = FALSE)
      PAM_roc = auroc(diablo_shuffles[[3]], plot = FALSE, print = FALSE)
      Spec_roc = auroc(diablo_shuffles[[4]], plot = FALSE, print = FALSE)
      
      #Extract top features of random models
      HC_shuffle_var = plotVar(diablo_shuffles[[1]], comp.select = c(1:6), plot = FALSE)
      HC_shuffle_features = HC_shuffle_var$names
      
      KM_shuffle_var = plotVar(diablo_shuffles[[2]], comp.select = c(1:6), plot = FALSE)
      KM_shuffle_features = KM_shuffle_var$names
      
      PAM_shuffle_var = plotVar(diablo_shuffles[[3]], comp.select = c(1:6), plot = FALSE)
      PAM_shuffle_features = PAM_shuffle_var$names
      
      Spec_shuffle_var = plotVar(diablo_shuffles[[4]], comp.select = c(1:6), plot = FALSE)
      Spec_shuffle_features = Spec_shuffle_var$names
      
      #Intersect of features between original model and randomized model
      HC_inter = intersect(HC_features, HC_shuffle_features)
      KM_inter = intersect(KM_features, KM_shuffle_features)
      PAM_inter = intersect(PAM_features, PAM_shuffle_features)
      Spec_inter = intersect(Spec_features, Spec_shuffle_features)
      
      shuffle_res = c(mean(HC_roc$Comp4[,1]), length(HC_inter)/length(HC_features), HC_prop,
                      mean(KM_roc$Comp4[,1]), length(KM_inter)/length(KM_features), KM_prop,
                      mean(PAM_roc$Comp4[,1]), length(PAM_inter)/length(PAM_features), PAM_prop,
                      mean(Spec_roc$Comp4[,1]), length(Spec_inter)/length(Spec_features), Spec_prop)
      df_shuffle = rbind(df_shuffle,shuffle_res)
  }
  return(df_shuffle) 
}

# Output the graphs for the bootstrap method
output_permutation <- function(df_shuffle, dirname){
  
  # Output models of accuracy vs intersect of features, useless?
  ggplot(data = df_shuffle, aes(x = HC_accuracy, y = HC_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/HClust/HC_acc.jpeg"))
  ggplot(data = df_shuffle, aes(x = KM_accuracy, y = KM_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/Kmeans/KM_acc.jpeg"))
  ggplot(data = df_shuffle, aes(x = PAM_accuracy, y = PAM_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/PAM/PAM_acc.jpeg"))
  ggplot(data = df_shuffle, aes(x = Spec_accuracy, y = Spec_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/",dirname,"/Spectrum/Spec_acc.jpeg"))
  
  # Output models of proportion of randomized outcome vs prop of shared features
  df_shuffle$prop = rep(c(seq(from = 1, to = 0.1, by = -0.1), seq(from = 0.09, to = 0.01, by = -0.01)), each = 10)
  
  ggplot(data = df_shuffle, aes(x = prop, y = HC_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/HClust/HC_prop.jpeg"))
  ggplot(data = df_shuffle, aes(x = prop, y = KM_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/Kmeans/KM_prop.jpeg"))
  ggplot(data = df_shuffle, aes(x = prop, y = PAM_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/PAM/PAM_prop.jpeg"))
  ggplot(data = df_shuffle, aes(x = prop, y = Spec_intersect)) +
    geom_point() +
    xlim(1,0) + ylim(0,1) +
    theme_minimal()
  ggsave(paste0("./results/", dirname, "/Spectrum/Spec_prop.jpeg"))
}