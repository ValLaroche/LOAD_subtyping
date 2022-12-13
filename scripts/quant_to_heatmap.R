library(ggdendro)

quant_to_heatmap=function(res_clustering, discrimination, label_disc, label_disc2, 
                          transformation="scale", rainbow_color_palette, 
                          all_genes, output, sublist_biomart = NULL){
  
  tmp_counts = res_clustering[,2:ncol(res_clustering)]
  tmp_counts = data.frame(lapply(tmp_counts, function(x) {as.numeric(x)}))
  res_clustering = data.frame(cbind(res_clustering[,1], tmp_counts))
  colnames(res_clustering)[1] = "Transcript"
  
  scaled_res_clustering=res_clustering
  
  
  if(transformation=="scale"){
    
    tmp_matrix = scaled_res_clustering[,2:ncol(scaled_res_clustering)]
    
    tmp_matrix = t(tmp_matrix)
    tmp_matrix = scale(tmp_matrix, center = TRUE, scale = TRUE)
    tmp_matrix = t(tmp_matrix)
    
    tmp_matrix[is.na(tmp_matrix)]=0
    tmp_matrix[tmp_matrix > 2]=2
    tmp_matrix[tmp_matrix < -2]=-2
    
    scaled_res_clustering = data.frame(scaled_res_clustering[,1], tmp_matrix)
    colnames(scaled_res_clustering)[1] = "Transcript"
  }
  if(transformation=="log"){
    scaled_res_clustering[, c(2:ncol(res_clustering))]=
      log2(res_clustering[, 2:ncol(res_clustering)]+1)
    scaled_res_clustering[is.na(scaled_res_clustering)]=0
    scaled_res_clustering=scaled_res_clustering %>% 
      mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  }
  
  
  
  new_res_clustering=as.matrix(scaled_res_clustering[-1])
  rownames(new_res_clustering)=res_clustering$Transcript
  quant_h_clust=hclust(d=dist(x=new_res_clustering),method="ward.D2")
  quant_dendro=as.dendrogram(quant_h_clust)
  
  dendro_plot=ggdendrogram(data=quant_dendro, rotate=TRUE)
  
  data_dendro=dendro_data(quant_dendro)     
  
  labels = data_dendro$labels
  labels$label = as.character(labels$label)
  
  melted_res_clustering=melt(scaled_res_clustering, id.vars = "Transcript")
  
  dendro_plot=dendro_plot + 
    theme_void()+
    scale_x_continuous(expand = c(0.5/length(unique(melted_res_clustering$Transcript)), 
                                  0.5/length(unique(melted_res_clustering$Transcript)))) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none')
  
  transposed_res_clustering=t(as.matrix(new_res_clustering))
  label_h_clust=hclust(d=dist(transposed_res_clustering),method="ward.D2")
  label_dendro=as.dendrogram(label_h_clust)
  
  clusters_evaluate = cutree(label_dendro, k=2)
  conds_evaluate = gsub("\\_.*","",names(clusters_evaluate))
  
  clust_eval_k2 = round(cluster.evaluation(conds_evaluate, clusters_evaluate), digits = 3)
  
  clusters_evaluate = cutree(label_dendro, k=3)
  clust_eval_k3 = round(cluster.evaluation(conds_evaluate, clusters_evaluate), digits = 3)
  clusters_evaluate = cutree(label_dendro, k=4)
  clust_eval_k4 = round(cluster.evaluation(conds_evaluate, clusters_evaluate), digits = 3)
  
  clust_eval_hmap = data.frame("Groups" = c("k = 2", "k = 3", "k = 4"),
                               "Similarity_index" = c(clust_eval_k2, clust_eval_k3, clust_eval_k4))
  
  clust_eval_hmap = tableGrob(clust_eval_hmap, theme = ttheme_default(base_size = 14), rows = NULL)
  
  label_dendro_plot=ggdendrogram(data=label_dendro)
  
  label_data_dendro=dendro_data(label_dendro)     
  
  row_labels = label_data_dendro$labels
  row_labels$label = as.character(row_labels$label)
  
  label_dendro_plot=label_dendro_plot +
    theme_void()+
    scale_x_continuous(expand = c(0,0))+
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none')
  
  colnames(melted_res_clustering)[1] = "Transcript"
  melted_res_clustering$Transcript=factor(x=melted_res_clustering$Transcript, 
                                        levels=unique(scaled_res_clustering$Transcript[order.dendrogram(quant_dendro)]))
  
  melted_res_clustering$variable=factor(x=melted_res_clustering$variable, 
                                      levels=rownames(transposed_res_clustering[order.dendrogram(label_dendro),]))
  
  res=ggplot(melted_res_clustering, aes(x=variable, y=Transcript, fill=value))+
    geom_tile(size=0.2)+
    {if(transformation == "log")scale_fill_gradientn(colours  = rainbow_color_palette, na.value = "black")} +
    {if(transformation == "scale")scale_fill_gradient2(low = rainbow_color_palette[[1]], 
                                                       mid = rainbow_color_palette[[2]], 
                                                       high = rainbow_color_palette[[3]])} +
    theme_void()+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none')
  
  panel=as.factor(as.vector(sapply(melted_res_clustering$variable,function(x) 
    paste(strsplit(as.character(x),"_")[[1]][discrimination],collapse="_"))))
  cond1=lapply(panel, function(x){
    strsplit(as.character(x), '_')[[1]][1]
  })
  cond2=lapply(panel, function(x){
    strsplit(as.character(x), '_')[[1]][2]
  })
  cond3=lapply(panel, function(x){
    strsplit(as.character(x), '_')[[1]][3]
  })
  
  panel_df=data.frame(cond1=as.character(cond1), 
                      cond2=as.character(cond2), 
                      cond3=as.character(cond3), 
                      sample=melted_res_clustering$variable)
  
  cond1_nb_cols=length(unique(cond1))
  
  cond1_colors=colorRampPalette(c("coral","palegreen2","mediumpurple4","lightblue"))(cond1_nb_cols)
  
  cond1_label=ggplot(panel_df)+
    geom_bar(aes(x=sample, y=1, fill=cond1), 
             stat='identity', width=1)+
    theme_void()+
    labs(fill=label_disc)+
    theme(panel.spacing.x = unit(1, "mm"))+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values=cond1_colors)+
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none')
  
  legend_cond1_label=ggplot(panel_df)+geom_bar(aes(x=sample, y=1, fill=cond1), 
                                               stat='identity', width=1)+
    theme_void()+
    labs(fill=label_disc)+
    theme(panel.spacing.x = unit(1, "mm"), 
          legend.text = element_text(size = 16), 
          legend.title = element_text(size = 17))+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))+
    scale_fill_manual(values=cond1_colors)
  
  cond2_nb_cols=length(unique(cond2))
  
  cond2_colors=colorRampPalette(c("darkorchid3","darkorange3"))(cond2_nb_cols)
  
  cond2_label=ggplot(panel_df)+
    geom_bar(aes(x=sample, y=1, fill=cond2), 
             stat='identity', width=1)+
    theme_void()+
    labs(fill=label_disc2)+
    theme(panel.spacing.x = unit(1, "mm"))+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values=cond2_colors)+
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none')
  
  legend_cond2_label=ggplot(panel_df)+geom_bar(aes(x=sample, y=1, fill=cond2), 
                                               stat='identity', width=1)+
    theme_void()+
    labs(fill=label_disc2)+
    theme(panel.spacing.x = unit(1, "mm"), 
          legend.text = element_text(size = 16), 
          legend.title = element_text(size = 17))+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))+
    scale_fill_manual(values=cond2_colors)
  
  path_label=ggplot(panel_df) +
    geom_bar(aes(x=sample, y=1, fill=cond1), 
             stat='identity', width=1) +
    theme_void() +
    labs(fill=label_disc) +
    theme(panel.spacing.x = unit(1, "mm")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values=cond1_colors) +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none')
  
  if(!all_genes){
    data_dendro$labels$label = sublist_biomart[match(data_dendro$labels$label, sublist_biomart$Transcript_stable_ID),]$Gene_name
    
    
    label_row=ggplot()+geom_text(data=data_dendro$labels, aes(x=x, label=label), 
                                 size=4, y=1,hjust = "inward",
                                 vjust = 0, angle = 0)+
      theme_void()+
      scale_x_continuous(expand = c(0.5/length(unique(melted_res_clustering$Transcript)), 
                                    0.5/length(unique(melted_res_clustering$Transcript)))) +
      coord_flip()
  } else {
    label_row = ggplot() + theme_minimal()
  }
  
  label_data_dendro$labels$label = lapply(label_data_dendro$labels$label, function(x){
    strsplit(as.character(x), '_')[[1]][3]
  })
  label_col=ggplot()+geom_text(data=label_data_dendro$labels, 
                               aes(x=x, label=label), 
                               size=4, y=1, hjust = 1,vjust = 0.5,angle = 90)+
    theme_void()+
    scale_x_continuous(expand = c(0.5/length(unique(melted_res_clustering$variable)), 
                                  0.5/length(unique(melted_res_clustering$variable)))) 
  
  plot_legend=ggplot(melted_res_clustering, aes(x=variable, y=Transcript)) +
    labs(fill="Relative expression") +
    geom_tile(aes(fill=value),color="black",size=0.1) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    {if(transformation == "log")scale_fill_gradientn(colours = rainbow_color_palette, na.value = "black")} +
    {if(transformation == "scale")scale_fill_gradient2(low = rainbow_color_palette[[1]], 
                                                       mid = rainbow_color_palette[[2]], 
                                                       high = rainbow_color_palette[[3]])}
  
  if(length(discrimination)==1){
    grob = list(res, label_dendro_plot, dendro_plot, label_col, label_row, 
                cond1_label, g_legend(plot_legend), 
                g_legend(legend_cond1_label), clust_eval_hmap)
    legend_hm=8
    legend_panel1=9
    legend_panel2=6
    legend_panel3=6
    clust_coeff=10
  }
  if(length(discrimination)==2){
    grob = list(res, label_dendro_plot, dendro_plot, label_col, label_row, 
                cond1_label, cond2_label, g_legend(plot_legend), 
                g_legend(legend_cond1_label), g_legend(legend_cond2_label), 
                clust_eval_hmap)
    
    panel_2=8
    legend_hm=9
    legend_panel1=10
    legend_panel2=11
    legend_panel3=6
    clust_coeff=12
  }
  # if(length(discrimination)==3){
  #   grob = list(res, label_dendro_plot, dendro_plot, label_col, label_row, 
  #               cond1_label, cond2_label, cond3_label, g_legend(plot_legend), 
  #               g_legend(legend_cond1_label), g_legend(legend_cond2_label), 
  #               g_legend(legend_cond3_label), clust_eval_hmap)
  #   
  #   panel_2=8
  #   panel_3=9
  #   legend_hm=10
  #   legend_panel1=11
  #   legend_panel2=12
  #   legend_panel3=13
  #   clust_coeff=14
  #   
  # }
  
  png(output, width = 1920, height = 1080, units = "px")
  
  final_plot = grid.arrange(
    grobs = grob,
    layout_matrix = rbind(c(NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,legend_hm,legend_hm),
                          c(NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,legend_hm,legend_hm),
                          c(NA,clust_coeff,NA,2,2,2,2,2,2,2,2,2,2,2,2,legend_hm,legend_hm),
                          c(NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,legend_hm,legend_hm),
                          c(NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,legend_hm,legend_hm),
                          c(NA,NA,NA,7,7,7,7,7,7,7,7,7,7,7,7,NA,NA),
                          if(length(discrimination)>=2){ c(NA,NA,NA,panel_2,panel_2,panel_2,panel_2,panel_2,panel_2,panel_2,panel_2,panel_2,panel_2,panel_2,panel_2,NA,NA)},
                          if(length(discrimination)==3){ 
                            c(NA,NA,NA,panel_3,panel_3,panel_3,panel_3,panel_3,panel_3,panel_3,panel_3,panel_3,panel_3,panel_3,panel_3,NA,NA)
                          },
                          c(legend_panel1,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(legend_panel1,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(legend_panel2,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(legend_panel2,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(legend_panel3,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(legend_panel3,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3),
                          c(NA,NA,NA,5,5,5,5,5,5,5,5,5,5,5,5,NA,NA),
                          c(NA,NA,NA,5,5,5,5,5,5,5,5,5,5,5,5,NA,NA),
                          c(NA,NA,NA,5,5,5,5,5,5,5,5,5,5,5,5,NA,NA),
                          c(NA,NA,NA,5,5,5,5,5,5,5,5,5,5,5,5,NA,NA)))
  
  dev.off()
  
}