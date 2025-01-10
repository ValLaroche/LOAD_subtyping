gen_eigengenes <- function(data_1, data_2, data_3 = NULL, net_colors){
  
  UKBBN_eigen = moduleEigengenes(data_1, net_colors)$eigengenes
  PITT_eigen = moduleEigengenes(data_2, net_colors)$eigengenes
  ret = list(UKBBN_eigen, PITT_eigen)
  if(!is.null(data_3)){
    ROSMAP_eigen = moduleEigengenes(data_3, net_colors)$eigengenes  
    ret = list(UKBBN_eigen, PITT_eigen,ROSMAP_eigen)
  }
  
  return(ret)
}

gen_phenos <- function(pheno_1, pheno_2, pheno_3, clustering){

  pheno_1 = pheno_1[pheno_1[[clustering]] != "Control",]
  pheno_2 = pheno_2[pheno_2[[clustering]] != "Control",]
  pheno_3 = pheno_3[pheno_3[[clustering]] != "Control",]
  return(list(pheno_1,pheno_2,pheno_3))
}

gen_subtype_matrix <- function(pheno, clustering, nsubtypes){
  
  if(nsubtypes == 3){
    pheno_trait = as.data.frame(cbind(pheno[[clustering]],pheno[[clustering]],pheno[[clustering]]))
    rownames(pheno_trait) = pheno$Basename
    colnames(pheno_trait) = c("sub1","sub2","sub3")
    
    pheno_trait[pheno_trait$sub1 == 1,]$sub1 = "1"
    pheno_trait[pheno_trait$sub1 != 1,]$sub1 = "0"
    pheno_trait[pheno_trait$sub2 != 2,]$sub2 = "0"
    pheno_trait[pheno_trait$sub2 == 2,]$sub2 = "1"
    pheno_trait[pheno_trait$sub3 != 3,]$sub3 = "0"
    pheno_trait[pheno_trait$sub3 == 3,]$sub3 = "1"
  } else if(nsubtypes == 2){
    pheno_trait = as.data.frame(cbind(pheno[[clustering]],pheno[[clustering]]))
    rownames(pheno_trait) = pheno$Basename
    colnames(pheno_trait) = c("sub1","sub2")
    
    pheno_trait[pheno_trait$sub1 != 1,]$sub1 = "0"
    pheno_trait[pheno_trait$sub1 == 1,]$sub1 = "1"
    pheno_trait[pheno_trait$sub2 != 2,]$sub2 = "0"
    pheno_trait[pheno_trait$sub2 == 2,]$sub2 = "1"
    
  }
  return(pheno_trait)
}

heatmap_module_to_subtype <- function(eigengenes, pheno_trait, filename){
  ME_eigen = orderMEs(eigengenes)
  ME_eigen = ME_eigen[rownames(ME_eigen) %in% rownames(pheno_trait),]
  moduleTraitCor <- cor(ME_eigen, pheno_trait, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(pheno_trait))
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  pdf(file = filename, width = 10, height = 10)
  ## display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(moduleTraitCor),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = rownames(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.7,
                 main = paste("Module-trait relationships"))
  dev.off()
}

gen_anova_tukey_tests <- function(pheno, clustering, eigengenes, filename){
  eigengenes = eigengenes[rownames(eigengenes) %in% pheno$Basename,]
  
  to_anova = cbind(pheno[[clustering]], eigengenes)
  colnames(to_anova)[1] = "clusters"

  tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
  for(i in c(2:(ncol(to_anova)))){
    
    tmptukey_table <- aov(to_anova[,i] ~ clusters, data = to_anova) %>%
      tukey_hsd() 
    
    vec_posy = seq(from = max(to_anova[,-1]), to = max(to_anova[,-1]*2), length.out = nrow(tmptukey_table))
    
    tmptukey_table = tmptukey_table %>% mutate(y.position = vec_posy)
    
    tukey_table = rbind(tukey_table, tmptukey_table)
  }
  tukey_table$variable = rep(colnames(to_anova)[c(2:ncol(to_anova))], each = nrow(tmptukey_table))
  
  to_anova$Basename = rownames(to_anova)
  to_anovamelt = melt(to_anova, id.vars = c("Basename","clusters"))
  
  pdf(file = filename, width = 18, height = 10)
  print(ggplot(to_anovamelt, aes(x = clusters, y = value, color = clusters)) +
    geom_boxplot() +
    geom_jitter(size = .5) +
    facet_wrap(~ variable, nrow = 3, scales = "free") +
    stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                       method = "anova", label.y = max(to_anova[, i]*1.2, na.rm = TRUE)) +
    stat_pvalue_manual(tukey_table, size = 3, label = "p.adj", tip.length = 0.03, bracket.shorten = 0.05) +
    theme_bw() +
    theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside"))
  dev.off()
}

double_anovas <- function(pheno_1, eigen_1, pref_1, pheno_2, eigen_2, pref_2, clustering, filename){
  eigen_1 = eigen_1[rownames(eigen_1) %in% pheno_1$Basename,]
  to_anova_1 = cbind(pheno_1[[clustering]], eigen_1)
  
  colnames(to_anova_1)[1] = "clusters"
  to_anova_1$clusters = paste0(pref_1, to_anova_1$clusters)
  
  eigen_2 = eigen_2[rownames(eigen_2) %in% pheno_2$Basename,]
  to_anova_2 = cbind(pheno_2[[clustering]], eigen_2)
  colnames(to_anova_2)[1] = "clusters"
  to_anova_2$clusters = paste0(pref_2, to_anova_2$clusters)
  
  full_anova = rbind(to_anova_1, to_anova_2)
  
  
  tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
  for(i in c(2:(ncol(full_anova)))){
    
    tmptukey_table <- aov(full_anova[,i] ~ clusters, data = full_anova) %>%
      tukey_hsd() 
  
    vec_posy = seq(from = max(full_anova[,-1]), to = max(full_anova[,-1]*2), length.out = nrow(tmptukey_table))
    
    tmptukey_table = tmptukey_table %>% mutate(y.position = vec_posy)
    tukey_table = rbind(tukey_table, tmptukey_table)
  }
  tukey_table$variable = rep(colnames(full_anova)[c(2:ncol(full_anova))], each = nrow(tmptukey_table))
  tukey_table = tukey_table[c(intersect(grep("EPIC1",tukey_table$group1),grep("450k1", tukey_table$group2)),
                              intersect(grep("450k1",tukey_table$group1),grep("EPIC1", tukey_table$group2)),
                              intersect(grep("EPIC2",tukey_table$group1),grep("450k2", tukey_table$group2)),
                              intersect(grep("450k2",tukey_table$group1),grep("EPIC2", tukey_table$group2)),
                              intersect(grep("EPIC3",tukey_table$group1),grep("450k3", tukey_table$group2)),
                              intersect(grep("450k3",tukey_table$group1),grep("EPIC3", tukey_table$group2))),]
  tukey_table$y.position = min(tukey_table$y.position)
  full_anovamelt = melt(full_anova, id.vars = "clusters")
  
  full_anovamelt_EPIC = full_anovamelt[grep("EPIC", full_anovamelt$clusters),]
  full_anovamelt_450k = full_anovamelt[grep("450k", full_anovamelt$clusters),]
  
  stat_EPIC = full_anovamelt_EPIC %>% group_by(variable) %>% anova_test(value ~ clusters)
  stat_EPIC$anovares = paste0("ANOVA EPIC, p-value :", stat_EPIC$p)
  stat_450k = full_anovamelt_450k %>% group_by(variable) %>% anova_test(value ~ clusters)
  stat_450k$anovares = paste0("ANOVA 450k, p-value :", stat_450k$p)
  
  full_anovamelt$clusters = factor(full_anovamelt$clusters, levels = c(paste0(pref_1, "1"),paste0(pref_2, "1"),
                                                                       paste0(pref_1, "2"),paste0(pref_2, "2"),
                                                                       paste0(pref_1, "3"),paste0(pref_2, "3")))
  pdf(file = filename, width = 9.5, height = 15.27)
  print(ggplot(full_anovamelt, aes(x = clusters, y = value, color = clusters)) +
    geom_boxplot() +
    geom_jitter(size = .5) +
    geom_text(aes(x = 3.6, y = 0.55, label = anovares), size = 2.7, data = stat_EPIC, inherit.aes = FALSE, position = position_dodge(.9)) +
    geom_text(aes(x = 3.6, y = 0.6, label = anovares), size = 2.7, data = stat_450k, inherit.aes = FALSE, position = position_dodge(.9)) +
    facet_wrap(~ variable, nrow = 4, scales = "free") +
    stat_pvalue_manual(tukey_table, size = 4, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
    theme_bw() +
    ylim(-0.4, 0.65) +
    xlab("WGCNA network modules") + ylab("Eigen values") + 
    theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
          axis.text.x=element_text(size=7,angle = 45, vjust = 1, hjust=1),
          strip.text = element_text(size = 11)))
  dev.off()
  
}

triple_anovas <- function(pheno_1, eigen_1, pref_1, 
                          pheno_2, eigen_2, pref_2, 
                          pheno_3, eigen_3, pref_3, 
                          pheno_4, eigen_4, pref_4, 
                          clustering, filename){
  eigen_1 = eigen_1[rownames(eigen_1) %in% pheno_1$Basename,]
  to_anova_1 = cbind(pheno_1[[clustering]], eigen_1)
  colnames(to_anova_1)[1] = "clusters"
  to_anova_1$clusters = paste0(pref_1, to_anova_1$clusters)
  
  eigen_2 = eigen_2[rownames(eigen_2) %in% pheno_2$Basename,]
  to_anova_2 = cbind(pheno_2[[clustering]], eigen_2)
  colnames(to_anova_2)[1] = "clusters"
  to_anova_2$clusters = paste0(pref_2, to_anova_2$clusters)
  
  eigen_3 = eigen_3[rownames(eigen_3) %in% pheno_3$Basename,]
  to_anova_3 = cbind(pheno_3[[clustering]], eigen_3)
  colnames(to_anova_3)[1] = "clusters"
  to_anova_3$clusters = paste0(pref_3, to_anova_3$clusters)
  
  eigen_4 = eigen_4[rownames(eigen_4) %in% pheno_4$Basename,]
  to_anova_4 = cbind(pheno_4[[clustering]], eigen_4)
  colnames(to_anova_4)[1] = "clusters"
  to_anova_4$clusters = paste0(pref_4, to_anova_4$clusters)
  
  full_anova = rbind(to_anova_1, to_anova_2, to_anova_4)
  
  
  tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
  for(i in c(2:(ncol(full_anova)))){
    
    tmptukey_table <- aov(full_anova[,i] ~ clusters, data = full_anova) %>%
      tukey_hsd() 
    
    vec_posy = seq(from = max(full_anova[,-1]), to = max(full_anova[,-1]*2), length.out = nrow(tmptukey_table))
    
    tmptukey_table = tmptukey_table %>% mutate(y.position = vec_posy)
    tukey_table = rbind(tukey_table, tmptukey_table)
  }
  tukey_table$variable = rep(colnames(full_anova)[c(2:ncol(full_anova))], each = nrow(tmptukey_table))
  
  full_anovamelt = melt(full_anova, id.vars = "clusters")
  
  pdf(file = filename, width = 18, height = 10)
  print(ggplot(full_anovamelt, aes(x = clusters, y = value, color = clusters)) +
    geom_boxplot() +
    geom_jitter(size = .5) +
    facet_wrap(~ variable, nrow = 3, scales = "free") +
    stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                       method = "anova", label.y = max(full_anovamelt$value*1, na.rm = TRUE)) +
    stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5, hide.ns = TRUE) +
    theme_bw() +
    theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
          axis.text.x=element_text(size=7)))
  dev.off()
  
}


gen_comploading_plots <- function(pheno_1, eigen_1, cohort_1, 
                                  pheno_2, eigen_2, cohort_2, 
                                  pheno_3, eigen_3, cohort_3, 
                                  pheno_4, eigen_4, cohort_4, 
                                  clustering, title, filename){
  
  eigen_1 = eigen_1[rownames(eigen_1) %in% pheno_1$Basename,]
  splsda_1 = mixOmics::splsda(eigen_1, pheno_1[[clustering]], ncomp = 6)

  loadings_1 = data.frame(plotLoadings(splsda_1, comp = 1, contrib = "max", legend = FALSE))
  loadings_1$module = as.factor(rownames(loadings_1))
  
  loadings_1$cohort = cohort_1
  
  eigen_2 = eigen_2[rownames(eigen_2) %in% pheno_2$Basename,]
  splsda_2 = mixOmics::splsda(eigen_2, pheno_2[[clustering]], ncomp = 6)
 
  loadings_2 = data.frame(plotLoadings(splsda_2, comp = 1, contrib = "max", size.legend = 0.7, title = "HCUKBBN"))
  loadings_2$module = as.factor(rownames(loadings_2))
  
  loadings_2$cohort = cohort_2
  
  eigen_3 = eigen_3[rownames(eigen_3) %in% pheno_3$Basename,]
  splsda_3 = mixOmics::splsda(eigen_3, pheno_3[[clustering]], ncomp = 6)
  
  loadings_3 = data.frame(plotLoadings(splsda_3, comp = 1, contrib = "max", size.legend = 0.7, title = "HCUKBBN"))
  loadings_3$module = as.factor(rownames(loadings_3))
  
  loadings_3$cohort = cohort_3
  
  eigen_4 = eigen_4[rownames(eigen_4) %in% pheno_4$Basename,]
  splsda_4 = mixOmics::splsda(eigen_4, pheno_4[[clustering]], ncomp = 6)
  
  loadings_4 = data.frame(plotLoadings(splsda_4, comp = 1, contrib = "max", size.legend = 0.7, title = "HCUKBBN"))
  loadings_4$module = as.factor(rownames(loadings_4))
  
  loadings_4$cohort = cohort_4
  
  merge_loadings = rbind(loadings_1, loadings_2, loadings_3, loadings_4)
  
  merge_loadings$cohort = factor(merge_loadings$cohort, levels = c(cohort_1, cohort_2, cohort_3, cohort_4))
  
  pdf(file = filename, width = 18, height = 10)
  print(ggplot(merge_loadings, aes(x = X.importance, y = module, fill = X.GroupContrib)) +
    ggtitle(title) +
    geom_bar(stat = "identity") +
    xlab("Loadings") + ylab("Module") +
    facet_wrap(cohort~.) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

gen_random_anovas <- function(pheno_1, eigen_1, pref_1, 
                              pheno_2, eigen_2, pref_2, 
                              pheno_3, eigen_3, pref_3, 
                              pheno_4, eigen_4, pref_4, 
                              clustering, filename){
  pdf(file = filename, width = 18, height = 10)
  
  for(i in seq(1:5)){
    
    eigen_1 = eigen_1[rownames(eigen_1) %in% pheno_1$Basename,]
    random_1 = sample(pheno_1[[clustering]])
    to_anova_1 = cbind(random_1, eigen_1)
    colnames(to_anova_1)[1] = "clusters"
    to_anova_1$clusters = paste0(pref_1, to_anova_1$clusters)
    
    
    eigen_2 = eigen_2[rownames(eigen_2) %in% pheno_2$Basename,]
    random_2 = sample(pheno_2[[clustering]])
    to_anova_2 = cbind(random_2, eigen_2)
    colnames(to_anova_2)[1] = "clusters"
    to_anova_2$clusters = paste0(pref_2, to_anova_2$clusters)
    
    eigen_3 = eigen_3[rownames(eigen_3) %in% pheno_3$Basename,]
    random_3 = sample(pheno_3[[clustering]])
    to_anova_3 = cbind(random_3, eigen_3)
    colnames(to_anova_3)[1] = "clusters"
    to_anova_3$clusters = paste0(pref_3, to_anova_3$clusters)
    
    eigen_4 = eigen_4[rownames(eigen_4) %in% pheno_4$Basename,]
    random_4 = sample(pheno_4[[clustering]])
    to_anova_4 = cbind(random_4, eigen_4)
    colnames(to_anova_4)[1] = "clusters"
    to_anova_4$clusters = paste0(pref_4, to_anova_4$clusters)
    
    full_anova = rbind(to_anova_1, to_anova_2, to_anova_4)
    
    
    tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
    for(i in c(2:(ncol(full_anova)))){
      
      tmptukey_table <- aov(full_anova[,i] ~ clusters, data = full_anova) %>%
        tukey_hsd() 
      
      vec_posy = seq(from = max(full_anova[,-1]), to = max(full_anova[,-1]*2), length.out = nrow(tmptukey_table))
      
      tmptukey_table = tmptukey_table %>% mutate(y.position = vec_posy)
      tukey_table = rbind(tukey_table, tmptukey_table)
    }
    tukey_table$variable = rep(colnames(full_anova)[c(2:ncol(full_anova))], each = nrow(tmptukey_table))
    
    full_anovamelt = melt(full_anova, id.vars = "clusters")
    
    print(ggplot(full_anovamelt, aes(x = clusters, y = value, color = clusters)) +
      geom_boxplot() +
      geom_jitter(size = .5) +
      facet_wrap(~ variable, nrow = 3, scales = "free") +
      stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                         method = "anova", label.y = max(full_anovamelt$value*1.6, na.rm = TRUE)) +
      stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
      theme_bw() +
      theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
            axis.text.x=element_text(size=7)))
    
  }
  dev.off()
  
}

pred_mapping_1pred <- function(data_disc, pheno_disc, data_pred, pheno_pred, cohort, disc_cohort, clustering, filename){
  
  pheno_disc = pheno_disc[pheno_disc[[clustering]] != "Excluded",]
  data_disc = data_disc[rownames(data_disc) %in% pheno_disc$Basename,]
  
  pheno_pred = pheno_pred[pheno_pred[[clustering]] != "Excluded",]
  data_pred = data_pred[rownames(data_pred) %in% pheno_pred$Basename,]
  
  data_disc = data_disc[order(rownames(data_disc)),]
  pheno_disc = pheno_disc[order(pheno_disc$Basename),]
  pheno_pred = pheno_pred[order(pheno_pred$Basename),]
  data_pred = data_pred[order(rownames(data_pred)),]
  
  splsda_discovery = mixOmics::splsda(data_disc, pheno_disc[[clustering]], ncomp = 6)
  
  pred = predict(splsda_discovery, newdata = data_pred) 
  
  plot_disc = plotIndiv(splsda_discovery, comp = 1:2, rep.space = "X-variate",
                        style="graphics",ind.names=FALSE, legend = FALSE,
                        title = paste0(disc_cohort, " mapping in the ", cohort, " space"))
  
  
  pdf(filename, width = 8.27, height = 7.5)
  plot_disc = plotIndiv(splsda_discovery, comp = 1:2, rep.space = "X-variate",
                        style="graphics",ind.names=FALSE, legend = FALSE,
                        title = paste0(disc_cohort, " mapping in the ", cohort, " space"),
                        xlim = c(min(plot_disc$df$x)-20,max(plot_disc$df$x)+20),
                        ylim = c(min(plot_disc$df$y)-20,max(plot_disc$df$y)+20), pch = c(15, 16, 17))
  
  df_s1 = plot_disc$df[plot_disc$df$group == "1",]
  plot_convex_hull(df_s1$x, df_s1$y, lcolor = unique(df_s1$col))
  
  df_s2 = plot_disc$df[plot_disc$df$group == "2",]
  plot_convex_hull(df_s2$x, df_s2$y, lcolor = unique(df_s2$col))
  
  df_s3 = plot_disc$df[plot_disc$df$group == "3",]
  plot_convex_hull(df_s3$x, df_s3$y, lcolor = unique(df_s3$col))
  
  points(pred$variates[, 1], pred$variates[, 2], 
         pch = as.numeric(pheno_pred[[clustering]])+3, cex = 1, col = pheno_pred[[clustering]])
  par(xpd = TRUE)
  legend("topleft",legend = c(1,2,3),
         pch = c(4,5,6), cex = 1, col = c(1,2,3), title = "ROSMAP", inset=c(0,0.2))
  legend("topleft",legend = c(1,2,3),
         pch = c(15,16,17), cex = 1, col = c("#388ECC","#F68B33","#C2C2C2"), title = cohort)
  
  dev.off()
}

test_chulls <- function(data_disc, pheno_disc, data_pred, pheno_pred, clustering, disc){
  
  pheno_disc = pheno_disc[pheno_disc[[clustering]] != "Excluded",]
  data_disc = data_disc[rownames(data_disc) %in% pheno_disc$Basename,]
  pheno_pred = pheno_pred[pheno_pred[[clustering]] != "Excluded",]
  data_pred = data_pred[rownames(data_pred) %in% pheno_pred$Basename,]
  
  data_disc = data_disc[order(rownames(data_disc)),]
  pheno_disc = pheno_disc[order(pheno_disc$Basename),]
  pheno_pred = pheno_pred[order(pheno_pred$Basename),]
  data_pred = data_pred[order(rownames(data_pred)),]
  
  splsda_discovery = mixOmics::splsda(data_disc, pheno_disc[[clustering]], ncomp = 6)
  
  pred = predict(splsda_discovery, newdata = data_pred) 
  
  plot_disc = plotIndiv(splsda_discovery, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
  
  df_subtype = plot_disc$df[plot_disc$df$group == subtype_hull,]
  xp_s1 = df_subtype[chull(df_subtype$x, df_subtype$y),]$x
  yp_s1 = df_subtype[chull(df_subtype$x, df_subtype$y),]$y
  
  pred_s1 = pred$variates[rownames(pred$variates) %in% pheno_pred[pheno_pred[[clustering]] == 1,]$Basename,]
  inpoly_s1 = rownames(pred_s1[inpolygon(x = pred_s1[,1], y = pred_s1[,2],
                  xp = xp_s1, yp = yp_s1),])
  pred_s2 = pred$variates[rownames(pred$variates) %in% pheno_pred[pheno_pred[[clustering]] == 2,]$Basename,]
  inpoly_s2 = rownames(pred_s2[inpolygon(x = pred_s2[,1], y = pred_s2[,2],
                                         xp = xp_s1, yp = yp_s1),])
  pred_s3 = pred$variates[rownames(pred$variates) %in% pheno_pred[pheno_pred[[clustering]] == 3,]$Basename,]
  inpoly_s3 = rownames(pred_s3[inpolygon(x = pred_s3[,1], y = pred_s3[,2],
                                         xp = xp_s1, yp = yp_s1),])

  sum_inpoly = c(inpoly_s1,inpoly_s2,inpoly_s3)
  
  tmp_GO = GeneOverlap::newGeneOverlap(listA = sum_inpoly,
                                       listB = pheno_pred[pheno_pred[[clustering]] == 1,]$Basename,
                                       genome.size = nrow(pheno_pred))
  fisher_s1 = testGeneOverlap(tmp_GO)
  
  tmp_GO = GeneOverlap::newGeneOverlap(listA = sum_inpoly,
                                       listB = pheno_pred[pheno_pred[[clustering]] == 2,]$Basename,
                                       genome.size = nrow(pheno_pred))
  fisher_s2 = testGeneOverlap(tmp_GO)
  
  tmp_GO = GeneOverlap::newGeneOverlap(listA = sum_inpoly,
                                       listB = pheno_pred[pheno_pred[[clustering]] == 3,]$Basename,
                                       genome.size = nrow(pheno_pred))
  fisher_s3 = testGeneOverlap(tmp_GO)

  
  res_s1 = c(nrow(pheno_pred),
             nrow(pheno_pred[pheno_pred[[clustering]] == 1,]),
             length(sum_inpoly),
             fisher_s1@pval,
             fisher_s1@odds.ratio)
  res_s2 = c(nrow(pheno_pred),
             nrow(pheno_pred[pheno_pred[[clustering]] == 2,]),
             length(sum_inpoly),
             fisher_s2@pval,
             fisher_s2@odds.ratio)
  res_s3 = c(nrow(pheno_pred),
             nrow(pheno_pred[pheno_pred[[clustering]] == 3,]),
             length(sum_inpoly),
             fisher_s3@pval,
             fisher_s3@odds.ratio)
  
  df_res_fishers = data.frame(matrix(ncol = 5, nrow = 0))
  df_res_fishers = rbind(df_res_fishers, res_s1, res_s2, res_s3)
  
  colnames(df_res_fishers) = c("Total_cohort","Total_in_cluster","Total_in_hull","P_value_Fisher","Odds_ratio_Fisher")
  rownames(df_res_fishers) = c(paste0("Hull_S",subtype_hull,"_",disc,"_Subtype_1"),
                               paste0("Hull_S",subtype_hull,"_",disc,"_Subtype_2"),
                               paste0("Hull_S",subtype_hull,"_",disc,"_Subtype_3"))

  return(df_res_fishers)
}


pred_mapping_multpred <- function(data_disc, pheno_disc, 
                               data_pred_1, pheno_pred_1, 
                               data_pred_2, pheno_pred_2, 
                               data_pred_3 = NULL, pheno_pred_3 = NULL, 
                               clustering, filename){
  
  pheno_disc = pheno_disc[pheno_disc[[clustering]] != "Excluded",]
  data_disc = data_disc[rownames(data_disc) %in% pheno_disc$Basename,]
  
  pheno_pred_1 = pheno_pred_1[pheno_pred_1[[clustering]] != "Excluded",]
  data_pred_1 = data_pred_1[rownames(data_pred_1) %in% pheno_pred_1$Basename,]
  
  pheno_pred_2 = pheno_pred_2[pheno_pred_2[[clustering]] != "Excluded",]
  data_pred_2 = data_pred_2[rownames(data_pred_2) %in% pheno_pred_2$Basename,]
  
  
  splsda_discovery = mixOmics::splsda(data_disc, pheno_disc[[clustering]], ncomp = 6)
  
  pred_1 = predict(splsda_discovery, newdata = data_pred_1) 
  
  pred_2 = predict(splsda_discovery, newdata = data_pred_2) 
  
  if(!is.null(data_pred_3)){
    pheno_pred_3 = pheno_pred_3[pheno_pred_3[[clustering]] != "Excluded",]
    data_pred_3 = data_pred_3[rownames(data_pred_3) %in% pheno_pred_3$Basename,]
    
    
    pred_3 = predict(splsda_discovery, newdata = data_pred_3) 
  }
  pdf(filename, width = 15, height = 15)
  plotIndiv(splsda_discovery, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
  points(pred_1$variates[, 1], pred_1$variates[, 2], 
         pch = 19, cex = 1.2, col = pheno_pred_1[[clustering]])
  points(pred_2$variates[, 1], pred_2$variates[, 2], 
         pch = 3, cex = 1.2, col = pheno_pred_2[[clustering]])
  if(!is.null(data_pred_3)){
    points(pred_3$variates[, 1], pred_3$variates[, 2], 
           pch = 7, cex = 1.2, col = pheno_pred_3[[clustering]])
  }
  
  plotIndiv(splsda_discovery, comp = 3:4, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
  points(pred_1$variates[, 3], pred_1$variates[, 4], 
         pch = 19, cex = 1.2, col = pheno_pred_1[[clustering]])
  points(pred_2$variates[, 3], pred_2$variates[, 4], 
         pch = 3, cex = 1.2, col = pheno_pred_2[[clustering]])
  
  if(!is.null(data_pred_3)){
  points(pred_3$variates[, 3], pred_3$variates[, 4], 
         pch = 7, cex = 1.2, col = pheno_pred_3[[clustering]])
  }
  dev.off()
}

median_samples <- function(pheno_1, data_1, pref_1, pheno_2, data_2, pref_2,
                           pheno_3, data_3, pref_3, clustering){
  cor_median = data.frame(matrix(NA, nrow = 0, ncol = ncol(data_1)))
  allrownames = c()
  for(subtype in unique(pheno_1[[clustering]])){
    if(subtype %in% unique(pheno_1[[clustering]])){
      subtype_samples = pheno_1[pheno_1[[clustering]] == subtype,]$Basename
      to_median = data_1[,colnames(data_1) %in% subtype_samples]
      to_median_1 = rowMedians(as.matrix(to_median))
      cor_median = rbind(cor_median, to_median_1) 
      allrownames = c(allrownames, paste0(pref_1, "_S", subtype))
                      
    }
    if(subtype %in% unique(pheno_2[[clustering]])){
      subtype_samples = pheno_2[pheno_2[[clustering]] == subtype,]$Basename
      to_median = data_2[,colnames(data_2) %in% subtype_samples]
      to_median_2 = rowMedians(as.matrix(to_median))
      cor_median = rbind(cor_median, to_median_2)
      allrownames = c(allrownames, paste0(pref_2, "_S", subtype))
      
    }
    if(subtype %in% unique(pheno_3[[clustering]])){
      subtype_samples = pheno_3[pheno_3[[clustering]] == subtype,]$Basename
      to_median = data_3[,colnames(data_3) %in% subtype_samples]
      to_median_3 = rowMedians(as.matrix(to_median))
      cor_median = rbind(cor_median, to_median_3)
      allrownames = c(allrownames, paste0(pref_3, "_S", subtype))
      
    }
  }
  
  rownames(cor_median) = allrownames
  colnames(cor_median) = rownames(data_1)
  cor_median = cor_median[complete.cases(cor_median),]
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(cor(t(cor_median)), method = "color", order = "hclust", hclust.method = "ward.D2",
           col = col(200), tl.col = "black", addCoef.col = "black", insig = "blank")
  
  return(cor_median)
}
  
red_vs_blue <- function(){
  plot_disc = plotIndiv(splsda_discovery, comp = 3:4, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
  
  df_s1 = plot_disc$df[plot_disc$df$group == "1",]
  plot_convex_hull(df_s1$x, df_s1$y, lcolor = unique(df_s1$col))
  
  df_s2 = plot_disc$df[plot_disc$df$group == "2",]
  plot_convex_hull(df_s2$x, df_s2$y, lcolor = unique(df_s2$col))
  
  df_s3 = plot_disc$df[plot_disc$df$group == "3",]
  plot_convex_hull(df_s3$x, df_s3$y, lcolor = unique(df_s3$col))
  
  
  points(pred$variates[, 3], pred$variates[, 4], 
         pch = 19, cex = 1.2, col = pheno_pred[[clustering]])
  dev.off()
}
    

