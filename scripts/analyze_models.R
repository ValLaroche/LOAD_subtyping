plot_model <- function(model, model_type, label_samples, outdir){
  
  if(model_type == "Spectrum"){
    model_clusters = model$assignments
    model_matrix = model$similarity_matrix
    
    png(paste0(outdir, "eigenvector_analysis.png"), width = 2000, height = 1600)
    plot(model$eigenvector_analysis$K, model$eigenvector_analysis$evals)
    dev.off()
    
  } else if(model_type == "CC"){
    model_clusters = model[[4]]$consensusClass
    model_matrix = model[[4]]$consensusMatrix
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
    
    png(paste0(outdir,"clusters_heatmap.png"), width = 2000, height = 1600)
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
    png(paste0(outdir, "silhouette_plot.png"), width = 2000, height = 1600)
    plot(silhouette(model_clusters, model_matrix))
    dev.off()
    
    plot_clusters = model_heatmap(model_matrix = model_matrix,
                                  model_clusters = model_clusters,
                                  label = label_samples)
    
    png(paste0(outdir, "clusters_heatmap.png"), width = 2000, height = 1600)
    plot(plot_clusters)
    dev.off()
  }
  
}

models_tables <- function(df_assignment, df_NMI, outdir){
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
  
  plot_clusters = ggplot(data = df_NMI, aes(x = model_1, y = model_2, fill = nmi)) +
    ggtitle("NMI between models") +
    geom_tile() +
    geom_text(aes(label = round(nmi, 2))) +
    scale_y_discrete(limits = rev(levels(factor(colnames(df_assignment), levels = colnames(df_assignment))))) +
    scale_x_discrete(limits = levels(factor(colnames(df_assignment), levels = colnames(df_assignment)))) +
    scale_fill_gradient2(low="#ffffb2", mid="#fd8d3c", high="#bd0026", #colors in the scale
                         midpoint = 0.5, limits=c(0, 1)) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 14),
          panel.border = element_blank(),
          legend.position = "none")
  
  tables_res[[length(tables_res)+1]] = plot_clusters
  
  if(length(tables_res) == 4) {
    final_plot = arrangeGrob(
      grobs = tables_res,
      layout_matrix = rbind(c(4,4,4,4,1),
                            c(4,4,4,4,2),
                            c(4,4,4,4,3)),
      top = textGrob("Comparative results of clustering methods", gp = gpar(fontsize = 16))
    )
    ggsave(paste0(outdir, "Comparative_clustering.png"), final_plot, units = "px", width  = 2000, height = 1000)
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


clinical_compare <- function(df_assignment, BDR_metadata, df_NMI, outdir){
  tables_res = list()
  
  for(i in c(1:length(df_assignment))){
    
    table_res = tableGrob(table(BDR_metadata$diag, df_assignment[[i]]), 
                          theme = ttheme_minimal(base_size = 7, padding = unit(c(4,2),"mm"),
                                                 rowhead=list(fg_params=list(col="black", fontface=2L))))
    table_res = gtable_add_grob(table_res,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 2, b = nrow(table_res), l = 1, r = ncol(table_res))
    table_res = gtable_add_grob(table_res,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 1, b = nrow(table_res), l = 2, r = ncol(table_res))
    
    h = grobHeight(table_res)
    w = grobWidth(table_res)
    
    table1_label = textGrob(names(df_assignment)[i], x=unit(0.5, "npc") + w*0.1, 
                            y=unit(0.5,"npc") + 1.05*h, vjust=0, gp=gpar(fontsize=10, fontface = "bold"))
    
    final_res = gTree(children = gList(table_res, table1_label))
    tables_res[[length(tables_res)+1]] = final_res
  }
  
  
  plot_clusters = ggplot(data = df_NMI, aes(x = model_1, y = model_2, fill = nmi)) +
    ggtitle("NMI between models") +
    geom_tile() +
    geom_text(aes(label = round(nmi, 2))) +
    scale_fill_gradient2(low="#ffffb2", mid="#fd8d3c", high="#bd0026", #colors in the scale
                         midpoint = 0.5, limits=c(0, 1)) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5, size = 14),
          panel.border = element_blank(),
          legend.position = "none")

  tables_res[[length(tables_res)+1]] = plot_clusters
  
  final_plot = plot_clusters
  # 
  # if(length(tables_res) == 4) {
  #   final_plot = arrangeGrob(
  #     grobs = tables_res,
  #     layout_matrix = rbind(c(4,4,4,4,1),
  #                           c(4,4,4,4,2),
  #                           c(4,4,4,4,3)),
  #     top = textGrob("Models versus phenotype", gp = gpar(fontsize = 16))
  #   )
  #   ggsave(paste0(outdir, "Clinical_clustering.png"), final_plot, units = "px", width  = 2000, height = 1000)
  # } else if(length(tables_res) == 6){
  #   final_plot = arrangeGrob(
  #     grobs = tables_res,
  #     layout_matrix = rbind(c(6,6,6,6,1),
  #                           c(6,6,6,6,2),
  #                           c(6,6,6,6,3),
  #                           c(6,6,6,6,4),
  #                           c(6,6,6,6,5)),
  #     
  #     top = textGrob("Models versus phenotype", gp = gpar(fontsize = 16))
  #   )
    ggsave(paste0(outdir, "Clinical_clustering.png"), final_plot, units = "px", width  = 3000, height = 1500)
  # }
}


calcNMI_all <- function(df_assignment){
  df_NMI = data.frame(matrix(ncol = 3))
  colnames(df_NMI) = c("model_1","model_2","nmi")
  
  if(ncol(df_assignment >= 2)){
    row_index = 1
    for(i in c(1:ncol(df_assignment))){
      for(j in c(i:ncol(df_assignment))){  
        df_NMI[row_index,] = c(colnames(df_assignment)[i],
                               colnames(df_assignment)[j],
                               calNMI(df_assignment[,i], df_assignment[,j]))
        row_index = row_index + 1
      }
    }
  }
  
  df_NMI$nmi = as.numeric(df_NMI$nmi)
  return(df_NMI)
}

calcNMI_clinical <- function(df_assignment, BDR_metadata){
  df_NMI = data.frame(matrix(ncol = 3))
  colnames(df_NMI) = c("model_1","model_2","nmi")
  
  df_BDR = data.frame("Gender" = BDR_metadata$Gender,
                      "BraakTangle" = BDR_metadata$Braak_Tangle,
                      "Cerad" = BDR_metadata$CERAD,
                      "Thal" = BDR_metadata$Thal_phase,
                      "BraakLB" = BDR_metadata$Braak_LB,
                      "Diag" = BDR_metadata$diag,
                      "Dementia" = BDR_metadata$dementia)
  
  if(ncol(df_assignment >= 2)){
    row_index = 1
    for(i in c(1:ncol(df_assignment))){
      for(j in c(1:ncol(df_BDR))){  
        df_NMI[row_index,] = c(colnames(df_assignment)[i],
                               colnames(df_BDR)[j],
                               calNMI(df_assignment[,i], df_BDR[,j]))
        row_index = row_index + 1
      }
    }
  }
  
  df_NMI$nmi = as.numeric(df_NMI$nmi)
  return(df_NMI)
}
