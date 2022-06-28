##### Sourcing scripts and libraries #####
setwd("C:/Users/P70077107/Desktop/PhD")

sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}

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

##### Loading data #####

RNA_gene = load_data("data/new_gene_counts.txt")
RNA_transc = load_data("data/new_transcript_counts.txt")
miRNA = load_data("data/new_miRNA_counts.txt")
circRNA = load_data("data/new_circ_counts.txt")
proteomics = load_data("data/new_prot_counts.txt")
lncRNA = load_data("data/new_lncrna_counts.txt")

BDR_metadata = load_data("phenotyping/new_complete_metadata.txt")


##### Normalization and PCA #####

RNA_gene_norm = normalize_data(expr_matrix = RNA_gene, norm_method = "log2", is_RNA1 = FALSE, 
                               threshold_read_counts = 10, threshold_perc_samples = 0.8)
pca_analysis(as.matrix(RNA_gene_norm))

RNA_transc_norm = normalize_data(expr_matrix = RNA_transc, norm_method = "log2", is_RNA1 = FALSE, 
                                 threshold_read_counts = 10, threshold_perc_samples = 0.8)
pca_analysis(as.matrix(RNA_transc_norm))

miRNA_norm = normalize_data(expr_matrix = miRNA, norm_method = "log2", is_RNA1 = FALSE, 
                            threshold_read_counts = 10, threshold_perc_samples = 0.8)
pca_analysis(as.matrix(miRNA_norm))

circRNA_norm = normalize_data(expr_matrix = circRNA, norm_method = "log2", is_RNA1 = FALSE, 
                              threshold_read_counts = 3, threshold_perc_samples = 0.8)
pca_analysis(as.matrix(circRNA_norm))

proteomics[is.na(proteomics)] <- 0
proteomics_norm = normalize_data(expr_matrix = proteomics, norm_method = "log2", is_RNA1 = FALSE, 
                              threshold_read_counts = 10, threshold_perc_samples = 0.8)       
pca_analysis(as.matrix(proteomics_norm))

lncRNA_norm = normalize_data(expr_matrix = lncRNA, norm_method = "log2", is_RNA1 = FALSE, 
                             threshold_read_counts = 10, threshold_perc_samples = 0.8)
pca_analysis(as.matrix(lncRNA_norm))


G_aff_matrix = aff_matrix_calc(RNA_gene_norm)
T_aff_matrix = aff_matrix_calc(RNA_transc_norm)
M_aff_matrix = aff_matrix_calc(miRNA_norm)
C_aff_matrix = aff_matrix_calc(circRNA_norm)
P_aff_matrix = aff_matrix_calc(proteomics_norm)
L_aff_matrix = aff_matrix_calc(lncRNA_norm)

##### Building models #####

combinations = list(c(T,F,F,F,F), c(F,T,F,F,F), c(F,F,T,F,F), c(F,F,F,T,F), c(F,F,F,F,T),
                    c(T,T,F,F,F),c(T,F,T,F,F),c(T,F,F,T,F),c(T,F,F,F,T),c(F,T,T,F,F),
                    c(F,T,F,T,F),c(F,T,F,F,T),c(F,F,T,T,F),c(F,F,T,F,T),c(F,F,F,T,T),
                    c(T,T,T,F,F),c(T,T,F,T,F),c(T,T,F,F,T),c(T,F,T,F,T),c(T,F,F,T,T),
                    c(F,T,T,T,F),c(F,T,T,F,T),c(F,T,F,T,T),c(F,F,T,T,T),c(T,T,T,T,F),c(T,T,T,F,T),
                    c(T,T,F,T,T),c(T,F,T,T,T),c(F,T,T,T,T),c(T,T,T,T,T))

dir.create("results")

fulldf_allclust = data.frame(matrix(nrow = 92, ncol = 0))

for(combination in combinations){
  print(combination)
  list_datasets = list()
  list_aff_matrix = list()
  list_keep_vals = list()
  
  cols_df = ""
  dirname = ""
  if (combination[1]) { 
    list_datasets[[length(list_datasets)+1]] = RNA_transc_norm
    names(list_datasets)[length(list_datasets)] = "RNA_transc"
    
    list_aff_matrix[[length(list_aff_matrix)+1]] = T_aff_matrix
    names(list_aff_matrix)[length(list_aff_matrix)] = "RNA_transc"
    
    list_keep_vals[[length(list_keep_vals)+1]] = c(50,50)
    names(list_keep_vals)[length(list_keep_vals)] = "RNA_transc"
    
    dirname = paste(dirname, "RNAtransc", sep = "_")
    cols_df = paste0(cols_df, "T")
  }
  if (combination[2]) { 
    list_datasets[[length(list_datasets)+1]] = miRNA_norm
    names(list_datasets)[length(list_datasets)] = "miRNA"
    
    list_aff_matrix[[length(list_aff_matrix)+1]] = M_aff_matrix
    names(list_aff_matrix)[length(list_aff_matrix)] = "miRNA"
    
    list_keep_vals[[length(list_keep_vals)+1]] = c(50,50)
    names(list_keep_vals)[length(list_keep_vals)] = "miRNA"
    
    dirname = paste(dirname, "miRNA", sep = "_")
    cols_df = paste0(cols_df, "M")
  }
  if (combination[3]) { 
    list_datasets[[length(list_datasets)+1]] = circRNA_norm
    names(list_datasets)[length(list_datasets)] = "circRNA"
    
    list_aff_matrix[[length(list_aff_matrix)+1]] = C_aff_matrix
    names(list_aff_matrix)[length(list_aff_matrix)] = "circRNA"
    
    list_keep_vals[[length(list_keep_vals)+1]] = c(50,50)
    names(list_keep_vals)[length(list_keep_vals)] = "circRNA"
    
    dirname = paste(dirname, "circRNA", sep = "_")
    cols_df = paste0(cols_df, "C")
  }
  if (combination[4]) { 
    list_datasets[[length(list_datasets)+1]] = proteomics_norm
    names(list_datasets)[length(list_datasets)] = "proteomics"
    
    list_aff_matrix[[length(list_aff_matrix)+1]] = P_aff_matrix
    names(list_aff_matrix)[length(list_aff_matrix)] = "proteomics"
    
    list_keep_vals[[length(list_keep_vals)+1]] = c(50,50)
    names(list_keep_vals)[length(list_keep_vals)] = "proteomics"
    
    dirname = paste(dirname, "proteomics", sep = "_")
    cols_df = paste0(cols_df, "P")
  }
  if (combination[5]) { 
    list_datasets[[length(list_datasets)+1]] = lncRNA_norm
    names(list_datasets)[length(list_datasets)] = "lncRNA"
    
    list_aff_matrix[[length(list_aff_matrix)+1]] = L_aff_matrix
    names(list_aff_matrix)[length(list_aff_matrix)] = "lncRNA"
    
    list_keep_vals[[length(list_keep_vals)+1]] = c(50,50)
    names(list_keep_vals)[length(list_keep_vals)] = "lncRNA"
    
    dirname = paste(dirname, "lncRNA", sep = "_")
    cols_df = paste0(cols_df, "L")
  }
  
  dirname = gsub("^_", "", dirname)
  dir.create(paste0("results/", dirname))
  
  dir.create(paste0("results/", dirname, "/Spectrum"))
  multi_Spectrum = Spectrum_modelization(list_datasets, 
                                     spectrum_method = 1, cluster_alg = "GMM")
  plot_model(multi_Spectrum, "Spectrum", colnames(RNA_gene_norm), paste0("results/", dirname, "/Spectrum/"))
  
  df_assignment = data.frame("Spectrum" = multi_Spectrum$assignments,
                             row.names = names(multi_Spectrum$assignments))  

  dir.create(paste0("results/", dirname, "/ConsensusClustering"))
  multi_CC = CC_modelization(multi_Spectrum$similarity_matrix,
                             title_model = paste0("results/", dirname, "/ConsensusClustering/"), 
                             cluster_alg = "hc", distance = "pearson")
  plot_model(multi_CC, "CC", colnames(RNA_gene_norm), paste0("results/", dirname, "/ConsensusClustering/"))
  
  df_assignment$CC = multi_CC[[4]]$consensusClass
  
  dir.create(paste0("results/", dirname, "/iClusterBayes"))
  multi_iCB = ICB_modelization(list_datasets,
                               data_distrib = rep("gaussian", length(list_datasets)), K = 3)
  png(paste0("results/", dirname, "/iClusterBayes/feature_importance.png"), width = 2000, height = 1600)
  # plotHMBayes(multi_iCB, list_datasets, type = rep("gaussian", length(list_datasets)),
              # sample.order = multi_iCB$clusters, width = 20)
  dev.off()
  
  plot_model(multi_iCB, "iCB", colnames(RNA_gene_norm), paste0("results/", dirname, "/iClusterBayes/"))
  
  df_assignment$iCB = multi_iCB$clusters
  
  if(length(list_aff_matrix) >= 2){
    dir.create(paste0("results/", dirname, "/SNF"))
    multi_SNF = SNF_modelization(list_aff_matrix, K = 20, t = 20, num_clusters = 4)
    plot_model(multi_SNF, "SNF", colnames(RNA_gene_norm), paste0("results/", dirname, "/SNF/"))
    
    df_assignment$SNF = multi_SNF[[2]]
    
    dir.create(paste0("results/", dirname, "/ANF"))
    multi_ANF = ANF_modelization(list_aff_matrix, K = 20, num_clusters = 4)
    plot_model(multi_ANF, "ANF", colnames(RNA_gene_norm), paste0("results/", dirname, "/ANF/"))
    
    df_assignment$ANF = multi_ANF[[2]]

    # multi_JIVE = JIVE_modelization(datasets = list_datasets)
    
    dir.create(paste0("results/", dirname, "/SNF.CC"))
    multi_SNF.CC = SNF.CC_modelization(datasets = list_datasets, num_clusters = 4)
    plot_model(multi_SNF.CC, "SNF.CC", colnames(RNA_gene_norm), paste0("results/", dirname, "/SNF.CC/"))
    
  }
  
  df_NMI_models = calcNMI_all(df_assignment)
  
  models_tables(df_assignment, df_NMI_models, paste0("results/", dirname, "/"))
  
  df_NMI_clinical = calcNMI_clinical(df_assignment, BDR_metadata)
  df_NMI_clinical$model_2 = factor(df_NMI_clinical$model_2, levels = unique(df_NMI_clinical$model_2))
  clinical_compare(df_assignment, BDR_metadata, df_NMI_clinical, paste0("results/", dirname, "/"))
  
  colnames(df_assignment) = paste(cols_df, colnames(df_assignment), sep = "_")
    
  fulldf_allclust = cbind(fulldf_allclust, df_assignment)
}

################################## Testing ##################################

multi_JIVE = JIVE_modelization(datasets = list("RNA" = RNA_gene_norm,
                                               "miRNA" = miRNA_norm, 
                                               "circRNA" = circRNA_norm, 
                                               "proteomics" = proteomics_norm))

JIVE_analysis(jive_model = multi_JIVE)


X_diablo = list("RNA" = t(RNA_gene_norm),
                "miRNA" = t(miRNA_norm), 
                "circRNA" = t(circRNA_norm), 
                "proteomics" = t(proteomics_norm))
Y_diablo = BDR_metadata$Case.Control
keepX_diablo = list("RNA" = c(50, 50),
                    "miRNA" = c(50, 50),
                    "circRNA" = c(50, 50),
                    "proteomics" = c(50, 50))
multi_DIABLO = block.splsda(X = X_diablo,
                            Y = Y_diablo,
                            keepX = keepX_diablo)

plotIndiv(multi_DIABLO) ## sample plot
plotVar(multi_DIABLO, var.names = F, legend = TRUE, pch = c(20,20,20,20)) ## variable plot

plotIndiv(multi_DIABLO, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2,3), pch = c(15,18,20),
          title = 'BRCA with DIABLO')
 
Y = multi_DIABLO$Y
plotDiablo(multi_DIABLO, ncomp = 1)

circosPlot(multi_DIABLO, cutoff=0.8)

png("test.png", width = 3000, height = 1600)
lol = cimDiablo(multi_DIABLO, color.blocks = c('darkorchid', 'brown1', 'lightgreen', "dodgerblue2"), margins = c(5, 22), comp = 1, legend.position = "right")
dev.off()

plotLoadings(multi_DIABLO, comp = 1, contrib = "max")

network(multi_DIABLO, blocks = c(1,2,3,4),
        color.node = c('darkorchid', 'brown1', 'lightgreen', "dodgerblue2"), 
        cutoff = 0.6, save = 'png', name.save = 'DIABLOnetwork')


Y_unsup_DIABLO = as.data.frame(matrix(0, ncol = 3, nrow = ncol(RNA_gene)))
colnames(Y_unsup_DIABLO) = c("case","intermediate","control")

for (i in seq(1:nrow(BDR_metadata))){
  if (BDR_metadata$Case.Control[i] == "case"){
    Y_unsup_DIABLO[i,1] = 1
  } else if (BDR_metadata$Case.Control[i] == "intermediate"){
    Y_unsup_DIABLO[i,2] = 1
  } else if (BDR_metadata$Case.Control[i] == "control"){
    Y_unsup_DIABLO[i,3] = 1
  }
}
rownames(Y_unsup_DIABLO) = BDR_metadata$No.
rownames(Y_unsup_DIABLO) = gsub("10230_BBN_", "BDR_", rownames(Y_unsup_DIABLO))

Y_unsup_DIABLO = as.matrix(Y_unsup_DIABLO)

multi_unsup_DIABLO = block.spls(X = X_diablo,
                                Y = Y_unsup_DIABLO,
                                ncomp = 2,
                                keepX = keepX_diablo,
                                scale = TRUE,
                                tol = 1e-06,
                                max.iter = 100,
                                near.zero.var = FALSE,
                                all.outputs = TRUE)

plotIndiv(multi_unsup_DIABLO) ## sample plot
plotVar(multi_unsup_DIABLO, var.names = F, legend = TRUE, pch = c(20,20,20,20)) ## variable plot

plotIndiv(multi_unsup_DIABLO, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1), pch = c(15,18,20,4,5),
          title = 'BRCA with DIABLO')

Y = multi_unsup_DIABLO$Y
plotDiablo(multi_unsup_DIABLO, ncomp = 1)

plotLoadings(multi_unsup_DIABLO, comp = 1, contrib = "max")

network(multi_unsup_DIABLO, blocks = c(1,2,3,4),
        color.node = c('darkorchid', 'brown1', 'lightgreen', "dodgerblue2"), 
        cutoff = 0.94, save = 'png', name.save = 'DIABLOnetwork')



multi_iCB = iClusterBayes(dt1 = t(RNA_gene_norm),
              dt2 = t(miRNA_norm),
              dt3 = t(circRNA_norm),
              dt4 = t(proteomics_norm),
              type = c("gaussian","gaussian","gaussian","gaussian"),
              K = 2)

png("testicb.png")
plotHeatmap(multi_iCB,
            datasets = X_diablo,
            type = c("gaussian","gaussian","gaussian","gaussian"))
dev.off()

iCB_sim_matrix(multi_iCB, label = colnames(RNA_gene_norm))



result=ExecuteSNF.CC(datasets=datasets, clusterNum=4, K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500, 
                   title = "GBM", plot = "png", finalLinkage ="average")



group1=result$group
distanceMatrix1=result$distanceMatrix
p_value=survAnalysis(mainTitle="GBM1",time,status,group1,
                     distanceMatrix1,similarity=TRUE)
