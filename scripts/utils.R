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

iCB_sim_matrix = function (fit, label = NULL){
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
  plot_clusters = ggplot(data = melt_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_y_discrete(limits = rev(levels(factor(label, levels = label)))) +
    scale_x_discrete(limits = levels(factor(label, levels = label))) +
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
