iCB_sim_matrix <- function (fit, label = NULL){
  cl = fit$clusters
  sorted = sort(table(cl))
  o.stol = as.numeric(names(sorted))
  o = NULL
  for (i in o.stol) {
    o = c(o, which(cl == i))
  }
  s.matrix = fit$meanZ %*% t(fit$meanZ)
  diag.elements = diag(s.matrix)
  n = length(diag.elements)
  denom = matrix(rep(diag.elements, n), nrow = n, byrow = T)
  a = s.matrix/sqrt(denom)/sqrt(t(denom))
  a = replace(a, a < 0, 0)
  a = a[o, o]
  n = dim(a)[1]
  f.a = t(as.matrix(rev(as.data.frame(t(a)))))
  label = label[o]
  rownames(f.a) = rev(label)
  colnames(f.a) = label

  factor(label, levels = label)
  meltfa = melt(f.a)
  plot_clusters = ggplot(data = meltfa, aes(x = Var1, y = Var2, fill = value)) +
                         geom_tile() +
                         scale_y_discrete(limits = rev(levels(factor(label, levels = label)))) +
                         scale_x_discrete(limits = levels(factor(label, levels = label))) +
                         scale_fill_distiller(palette = "YlOrRd", direction = 1) +
                         theme_bw() + xlab("") + ylab("") +
                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(plot_clusters)
}

model_heatmap <- function (model_matrix, model_clusters, label = NULL){
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
  plot_clusters = ggplot(data = melt_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_y_discrete(limits = rev(levels(factor(label, levels = label)))) +
    scale_x_discrete(limits = levels(factor(label, levels = label))) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(plot_clusters)
}
