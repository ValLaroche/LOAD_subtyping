run_edgeR <- function(counts, pheno, subtype, cohort){
  
  # Generate surrogate variable on metadata
  # Optionnally include celltypes (Did not yeald great results ?)
  tmppheno = pheno[pheno$subtype %in% c(subtype, "Control", "Ctrl"),]
  if(cohort == "UKBBN"){
    tmpcounts = counts[,colnames(counts) %in% tmppheno$BBNId]
    group = factor(tmppheno$subtype, levels = c("Ctrl", subtype))
    tmppheno = calc_sva(tmpcounts, tmppheno, cohort)
    
    design = model.matrix(~group + Gender + Age + plate + Brain.Bank + 
                            # Astrocytes + Endothelial + Microglia + Neurons + OPC + Oligodendrocytes +
                            sv1 + sv2 + sv3 + sv4 + sv5, data = tmppheno)
    
  } else if(cohort == "PITT"){
    tmpcounts = counts[,colnames(counts) %in% tmppheno$Individual_ID]
    group = factor(tmppheno$subtype, levels = c("Control", subtype))
    colnames(tmppheno)[5] = "Age"

    tmppheno = calc_sva(tmpcounts, tmppheno, cohort)
    
    design = model.matrix(~group + Sex.y + Age + Plate + 
                            # Astrocytes + Endothelial + Microglia + Neurons + OPC + Oligodendrocytes +
                            sv1 + sv2 + sv3 + sv4 + sv5, data = tmppheno)
  }
  
  y = DGEList(counts=tmpcounts,group=group)
  # Filter for minimal expression and proportion of samples expressing
  keep = filterByExpr(y, min.count = 10, min.prop = 0.75)
  
  y = y[keep,,keep.lib.sizes=FALSE]
  # TMM normalization
  y = normLibSizes(y, method = "TMM")
  
  y = estimateDisp(y,design)

  fit = glmQLFit(y,design)
  qlf = glmQLFTest(fit,coef=2)
  
  final_table = topTags(qlf, n = nrow(tmpcounts))
  final_table = final_table$table
  
  return(final_table)
}

vplot_edgeR <- function(DEA_table, stat_value, stat_threshold, FC_threshold){
  tmpDEA_table = DEA_table
  tmpDEA_table$color = "black"
  tmpDEA_table[-log10(tmpDEA_table[stat_value]) > stat_threshold & tmpDEA_table$logFC < -FC_threshold,]$color = "green"
  tmpDEA_table[-log10(tmpDEA_table[stat_value]) > stat_threshold & tmpDEA_table$logFC > FC_threshold,]$color = "red"
  
  if(stat_value == "PValue"){
    ggplot(data = tmpDEA_table, aes(x = logFC, y = -log10(PValue), color = color)) +
      geom_point() +
      scale_color_manual(values=c("black","green", "red")) +
      geom_hline(yintercept = stat_threshold, color = "red") +
      geom_vline(xintercept = -FC_threshold, color = "red") +
      geom_vline(xintercept = FC_threshold, color = "red") +
      geom_text_repel(data = tmpDEA_table[tmpDEA_table$color %in% c("green", "red"),],
                      aes(label = ENSG),
                      size = 3.5) +
      theme_bw()
    
  } else if (stat_value == "FDR"){
    ggplot(data = tmpDEA_table, aes(x = logFC, y = -log10(FDR), color = color)) +
      geom_point() +
      scale_color_manual(values=c("black","green", "red")) +
      geom_hline(yintercept = stat_threshold, color = "red") +
      geom_vline(xintercept = -FC_threshold, color = "red") +
      geom_vline(xintercept = FC_threshold, color = "red") +
      geom_text_repel(data = tmpDEA_table[tmpDEA_table$color %in% c("green", "red"),],
                      aes(label = ENSG),
                      size = 2) +
      theme_bw()
    
  }

}

vplot_edgeR_intersect <- function(DEA_table, stat_value, stat_threshold, FC_threshold){
  tmpDEA_table = DEA_table
  
  tmpDEA_table$color = "black"
  tmpDEA_table[-log10(tmpDEA_table[stat_value]) > stat_threshold & 
                 tmpDEA_table$logFC < -FC_threshold,]$color = "green"
  tmpDEA_table[-log10(tmpDEA_table[stat_value]) > stat_threshold & tmpDEA_table$logFC > FC_threshold,]$color = "red"
  
  if(stat_value == "PValue"){
    ggplot(data = tmpDEA_table, aes(x = logFC, y = -log10(PValue), color = color)) +
      geom_point() +
      scale_color_manual(values=c("black","green", "red")) +
      geom_hline(yintercept = 2, color = "red") +
      geom_vline(xintercept = -0.5, color = "red") +
      geom_vline(xintercept = 0.5, color = "red") +
      geom_text_repel(data = tmpDEA_table[tmpDEA_table$color %in% c("green", "red"),],
                      aes(label = ENSG),
                      size = 3.5) +
      theme_bw()
    
  } else if (stat_value == "FDR"){
    ggplot(data = tmpDEA_table, aes(x = logFC, y = -log10(FDR), color = color)) +
      geom_point() +
      scale_color_manual(values=c("black","green", "red")) +
      geom_hline(yintercept = 2, color = "red") +
      geom_vline(xintercept = -0.5, color = "red") +
      geom_vline(xintercept = 0.5, color = "red") +
      geom_text_repel(data = tmpDEA_table[tmpDEA_table$color %in% c("green", "red"),],
                      aes(label = ENSG),
                      size = 3.5) +
      theme_bw()
    
  }
}