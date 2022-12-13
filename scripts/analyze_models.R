plot_model <- function(model, model_type, label_samples, outdir, numCC_clusters = 0){
  
  if(model_type == "Spectrum"){
    model_clusters = model$assignments
    model_matrix = model$similarity_matrix
  } else if(model_type == "CC"){
    model_clusters = model[[numCC_clusters]]$consensusClass
    model_matrix = model[[numCC_clusters]]$consensusMatrix
  } else if(model_type == "SNF"){
    model_clusters = model[[2]]
    model_matrix = model[[1]]
  } else if(model_type == "ANF"){
    model_clusters = model[[2]]
    model_matrix = model[[1]]
  } else if(model_type == "SNF.CC"){
    model_clusters = model$group
    model_matrix = model$distanceMatrix
  }
  
  if(model_type == "iCB"){
    plot_clusters = iCB_sim_matrix(model, label_samples)
    
    png(paste0(outdir,"clusters_heatmap.png"), width = 1000, height = 800)
    plot(plot_clusters)
    dev.off()
    
  } else if(model_type == "JIVE") {
    n_indiv = rep(1, length(jive_model$data))
    
    png(paste0(outdir, "VarExplained.png"),height=1000,width=1300)  
    showVarExplained(jive_model)  
    dev.off()
    
    png(paste0(outdir, "HeatmapsBRCA.png"),height=1000,width=1500)
    showHeatmaps(jive_model)
    dev.off()
    
    png(paste0(outdir, "JointPCA.png"),height=1000,width=1000)
    showPCA(jive_model, n_joint = 2)
    dev.off()
    
    png(paste0(outdir, "MorePCA.png"),height=1600,width=1600)
    showPCA(jive_model, n_joint=1, n_indiv= n_indiv)
    dev.off()
    
  } else {
    
    sil_model = silhouette(model_clusters, model_matrix)
    png(paste0(outdir, "silhouette_plot.png"), width = 1000, height = 800)
    plot(sil_model,col = sil_model[,1][order(sil_model[,1])])
    dev.off()
    
    plot_clusters = model_heatmap(model_matrix = model_matrix,
                                  model_clusters = model_clusters,
                                  label = label_samples)
    
    png(paste0(outdir, "clusters_heatmap.png"), width = 1000, height = 800)
    plot(plot_clusters)
    dev.off()
  }
  
}

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
    geom_text(aes(label = round(nmi, 2))) +
    scale_y_discrete(limits = rev(levels(factor(colnames(df_assignment), levels = colnames(df_assignment))))) +
    scale_x_discrete(limits = levels(factor(colnames(df_assignment), levels = colnames(df_assignment)))) +
    scale_fill_gradient2(low="#ffffb2", mid="#fd8d3c", high="#bd0026", #colors in the scale
                         midpoint = 0.5, limits=c(0, 1)) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 12),
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
    ggsave(paste0(outdir, "Comparative_clustering.png"), final_plot, units = "px", width  = 3000, height = 1500)
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
    ggsave(paste0(outdir, "Comparative_clustering.png"), final_plot, units = "px", width  = 3000, height = 1500)
  }
}


clinical_compare <- function(df_assignment, BDR_metadata, df_NMI_models, outdir){
  tables_res = list()
  
  melt_NMI_models = melt(df_NMI_models, id.vars = c("Model","Clinical_var"))
  
  dcast_models = dcast(data = melt_NMI_models,formula = Model~variable+Clinical_var,fun.aggregate = sum,value.var = "value")
  
  grid.table(melt_NMI_models)
  
  # example data & header row
  tab <- tableGrob(df_NMI_models)
  header <- tableGrob(df_NMI_models[1,], rows=NULL, cols = unique(df_NMI_models[,2])) 
  
  grid.draw(tab)
  
  tables_res[[length(tables_res)+1]] = plot_clusters
  
  final_plot = plot_clusters
  ggsave(paste0(outdir, "Clinical_clustering.png"), final_plot, units = "px", width  = 3000, height = 1500)
}


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

calcNMI_clinical <- function(df_assignment, BDR_metadata){
  df_NMI_models = data.frame(matrix(ncol = 8))
  colnames(df_NMI_models) = c("Model","Clinical_var","nmi","AMI","Chi2","RI","ARI", "spearman")
  
  df_BDR = data.frame("Gender" = BDR_metadata_train$Gender,
                      "Age" = BDR_metadata_train$Age,
                      "Age_conc" = BDR_metadata_train$age_conc,
                      "Plate" = BDR_metadata_train$Plate,
                      "Institute" = BDR_metadata_train$Institute,
                      "Diag" = BDR_metadata_train$diag,
                      "Dementia" = BDR_metadata_train$dementia,
                      "CoD" = BDR_metadata_train$Cause_of_Death,
                      "BraakTangle" = BDR_metadata_train$Braak_tangle,
                      "BraakLB" = BDR_metadata_train$Braak_LB,
                      "Cerad" = BDR_metadata_train$Cerad,
                      "APOE" = BDR_metadata_train$APOE,
                      "CDR" = BDR_metadata_train$CDR,
                      "MMSE" = BDR_metadata_train$MMSE)
  
  if(ncol(df_assignment_train >= 2)){
    row_index = 1
    for(i in c(1:ncol(df_assignment_train))){
      for(j in c(1:ncol(df_BDR))){  
        if(colnames(df_BDR)[j] %in% c("Age", "MMSE")){
          df_NMI_models[row_index,1] = colnames(df_assignment_train)[i]
          df_NMI_models[row_index,2] = colnames(df_BDR)[j]
          df_NMI_models[row_index,8] = cor(df_assignment_train[!is.na(df_BDR[,j]),i], df_BDR[!is.na(df_BDR[,j]),j], method = "spearman")
          row_index = row_index + 1
        } else {
          df_NMI_models[row_index,1] = colnames(df_assignment_train)[i]
          df_NMI_models[row_index,2] = colnames(df_BDR)[j]
          df_NMI_models[row_index,3] = NMI(df_assignment_train[!is.na(df_BDR[,j]),i], df_BDR[!is.na(df_BDR[,j]),j])
          df_NMI_models[row_index,4] = AMI(df_assignment_train[!is.na(df_BDR[,j]),i], df_BDR[!is.na(df_BDR[,j]),j])
          df_NMI_models[row_index,5] = Chi2(df_assignment_train[!is.na(df_BDR[,j]),i], df_BDR[!is.na(df_BDR[,j]),j])
          df_NMI_models[row_index,6] = RI(df_assignment_train[!is.na(df_BDR[,j]),i], df_BDR[!is.na(df_BDR[,j]),j])
          df_NMI_models[row_index,7] = ARI(df_assignment_train[!is.na(df_BDR[,j]),i], df_BDR[!is.na(df_BDR[,j]),j])
          row_index = row_index + 1
        }
      }
    }
  }
  
  return(df_NMI_models)
}

plot_mixOmics_all <- function(model, outdir){
  
  background = background.predict(model, comp.predicted=2, dist = "max.dist")
  
  png(paste0(outdir, "diablo_indiv_plotbackground12.png"), width = 1000, height = 800)
  plotIndiv(model,ind.names = FALSE, comp = c(1,2), 
            legend=TRUE,
            title = 'Methyl',
            background = background) 
  dev.off()
  
  png(paste0(outdir, "diablo_indiv_plotbackground13.png"), width = 1000, height = 800)
  plotIndiv(model,ind.names = FALSE, comp = c(1,3), 
            legend=TRUE,
            title = 'Methyl',
            background = background) 
  dev.off()
  
  png(paste0(outdir, "diablo_indiv_plotbackground14.png"), width = 1000, height = 800)
  plotIndiv(model,ind.names = FALSE, comp = c(1,4), 
            legend=TRUE,
            title = 'Methyl',
            background = background) 
  dev.off()
  
  png(paste0(outdir, "diablo_indiv_plotellipse.png"), width = 1000, height = 800)
  plotIndiv(model, # plot samples from final model
            ind.names = FALSE, # colour by class label
            ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
            title = 'Methyl')
  dev.off()
  
  png(paste0(outdir, "diablo_var_plot.png"), width = 1000, height = 800)
  plotVar(model)
  dev.off()
  
  png(paste0(outdir, "diablo_loading_comp1.png"), width = 1000, height = 800)
  plotLoadings(model, comp = 1, contrib = "max")
  dev.off()
  
  png(paste0(outdir, "diablo_loading_comp2.png"), width = 1000, height = 800)
  plotLoadings(model, comp = 2, contrib = "max")
  dev.off()
  
  png(paste0(outdir, "diablo_loading_comp3.png"), width = 1000, height = 800)
  plotLoadings(model, comp = 3, contrib = "max")
  dev.off()
  
  png(paste0(outdir, "diablo_loading_comp4.png"), width = 1000, height = 800)
  plotLoadings(model, comp = 4, contrib = "max")
  dev.off()
  
  png(paste0(outdir, "diablo_roc_curves.png"), width = 1000, height = 800)
  auroc(model, roc.comp = 4)
  dev.off()
  
  
  # set the styling of the legend to be homogeneous with previous plots
  legend=list(legend = levels(model$Y), # set of classes
              col = unique(color.mixo(model$Y)), # set of colours
              title = "Predicted subtypes", # legend title
              cex = 0.7) # legend size
  png(paste0(outdir, "diablo_heatmap.png"), width = 2000, height = 1500)
  # generate the CIM, using the legend and colouring rows by each sample's class
  
  cim(model, row.sideColors = color.mixo(model$Y), 
             legend = legend, mapping = "X", )
  
  dev.off()
}

