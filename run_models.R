##### Sourcing scripts, libraries and data #####
setwd("D:/valentin/main/")

sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}


load("betas/betas_UKBBN_BBcorr.Rdata")
load("betas/PITT_betas_platecorr.Rdata")
load("betas/ROSMAP_betas_platecorr.Rdata")
load("phenos_3cohorts.Rdata")
betas_UKBBN = betas_UKBBN_BBcorr
remove(betas_UKBBN_BBcorr)
##### Preparing data #####

cpgs_ROSMAP = rownames(betas_ROSMAP)
cpgs_UKBBN = rownames(betas_UKBBN)
cpgs_PITT = rownames(betas_PITT)
cpgs_intersect = Reduce(intersect, list(cpgs_PITT, cpgs_UKBBN))
cpgs_intersect_ROSMAP = intersect(cpgs_intersect, rownames(betas_ROSMAP))
##### PCA check #####
ROSMAP_PCA = prcomp(t(betas_ROSMAP))
autoplot(ROSMAP_PCA, data = pheno_ROSMAP, colour = "subtype")
autoplot(ROSMAP_PCA, data = pheno_ROSMAP, colour = "Sample_Plate")
autoplot(ROSMAP_PCA, data = pheno_ROSMAP, colour = "batch")
ROSMAP_PCAres = ROSMAP_PCA$x[order(rownames(ROSMAP_PCA$x)),]

identical(pheno_ROSMAP$Full_Sentrix, rownames(ROSMAP_PCAres))
pheno_ROSMAP[pheno_ROSMAP$age_first_ad_dx == "90+",]$age_first_ad_dx = 90
pheno_ROSMAP$age_first_ad_dx = as.numeric(pheno_ROSMAP$age_first_ad_dx)

pheno_ROSMAP = cbind(pheno_ROSMAP, ROSMAP_PCAres[,1:5])
ROSMAP_PCAmatrix = data.frame(matrix(nrow = 5, ncol = 0))
ROSMAP_txt = data.frame(matrix(nrow = 5, ncol = 0))
for(pheno in c("Gender", "batch","Age",
               "Neuronal_Inhibitory", "Neuronal_Excitatory", "Oligodendrocyte", "Microglia", "Astrocyte")){
  pheno_var = pheno_ROSMAP[[pheno]]
  if(typeof(pheno_var) == "character"){
    pheno_var = as.numeric(factor(pheno_var))
  }
  tmp_pheno = data.frame(cbind(ROSMAP_PCAres[,1:5],pheno_var))
  tmp_cors = c()
  tmp_res = c()
  for(pc in c(1:5)){

    cor_res = cor.test(tmp_pheno[,pc],tmp_pheno$pheno_var)
    tmp_cors = c(tmp_cors, cor_res$estimate)

    tmp_res = c(tmp_res,paste0(signif(cor_res$estimate, digits = 2), " (", signif(cor_res$p.value, digits = 2),")"))
  }
  ROSMAP_PCAmatrix = cbind(ROSMAP_PCAmatrix, tmp_cors)
  ROSMAP_txt = cbind(ROSMAP_txt, tmp_res)
}
colnames(ROSMAP_PCAmatrix) = c("Gender", "batch", "Age",
                               "Neuronal_Inhibitory", "Neuronal_Excitatory", "Oligodendrocyte", "Microglia", "Astrocyte")
rownames(ROSMAP_PCAmatrix) = paste0("PC", rownames(ROSMAP_PCAmatrix))
pheatmap::pheatmap(ROSMAP_PCAmatrix, treeheight_row = 0, treeheight_col = 0, 
                   color = redWhiteGreen(20), breaks = seq(-1, 1, by = 0.1),
                   cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = ROSMAP_txt, number_color = "black")

UKBBN_PCA = prcomp(t(betas_UKBBN_BBcorr))
autoplot(UKBBN_PCA, data = pheno_UKBBN, colour = "subtype")
autoplot(UKBBN_PCA, data = pheno_UKBBN, colour = "Brain.Bank")
UKBBN_PCAres = UKBBN_PCA$x[order(rownames(UKBBN_PCA$x)),]

UKBBN_PCAmatrix = data.frame(matrix(nrow = 5, ncol = 0))
UKBBN_txt = data.frame(matrix(nrow = 5, ncol = 0))
for(pheno in c("Gender", "Age", "Brain.Bank",
               "Neuronal_Inhibitory", "Neuronal_Excitatory", "Oligodendrocyte", "Microglia", "Astrocyte")){
  pheno_var = pheno_UKBBN[[pheno]]
  if(typeof(pheno_var) == "character"){
    pheno_var = as.numeric(factor(pheno_var))
  }
  tmp_pheno = data.frame(cbind(UKBBN_PCAres[,1:5],pheno_var))
  tmp_cors = c()
  tmp_res = c()
  for(pc in c(1:5)){
    
    cor_res = cor.test(tmp_pheno[,pc],tmp_pheno$pheno_var)
    tmp_cors = c(tmp_cors, cor_res$estimate)
    
    tmp_res = c(tmp_res,paste0(signif(cor_res$estimate, digits = 2), " (", signif(cor_res$p.value, digits = 2),")"))
  }
  UKBBN_PCAmatrix = cbind(UKBBN_PCAmatrix, tmp_cors)
  UKBBN_txt = cbind(UKBBN_txt, tmp_res)
}
colnames(UKBBN_PCAmatrix) = c("Gender", "Age", "Brain.Bank",
                               "Neuronal_Inhibitory", "Neuronal_Excitatory", "Oligodendrocyte", "Microglia", "Astrocyte")
rownames(UKBBN_PCAmatrix) = paste0("PC", rownames(UKBBN_PCAmatrix))
pheatmap::pheatmap(UKBBN_PCAmatrix, treeheight_row = 0, treeheight_col = 0, 
                   color = redWhiteGreen(20), breaks = seq(-1, 1, by = 0.1),
                   cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = UKBBN_txt, number_color = "black")

PITT_PCA = prcomp(t(betas_PITT))
autoplot(PITT_PCA, data = pheno_PITT, colour = "subtype")
PITT_PCAres = PITT_PCA$x[order(rownames(PITT_PCA$x)),]

PITT_PCAmatrix = data.frame(matrix(nrow = 5, ncol = 0))
PITT_txt = data.frame(matrix(nrow = 5, ncol = 0))
for(pheno in c("Sex", "Age", "Plate",
               "Neuronal_Inhibitory", "Neuronal_Excitatory", "Oligodendrocyte", "Microglia", "Astrocyte")){
  pheno_var = pheno_PITT[[pheno]]
  if(typeof(pheno_var) == "character"){
    pheno_var = as.numeric(factor(pheno_var))
  }
  tmp_pheno = data.frame(cbind(PITT_PCAres[,1:5],pheno_var))
  tmp_cors = c()
  tmp_res = c()
  for(pc in c(1:5)){
    
    cor_res = cor.test(tmp_pheno[,pc],tmp_pheno$pheno_var)
    tmp_cors = c(tmp_cors, cor_res$estimate)
    
    tmp_res = c(tmp_res,paste0(signif(cor_res$estimate, digits = 2), " (", signif(cor_res$p.value, digits = 2),")"))
  }
  PITT_PCAmatrix = cbind(PITT_PCAmatrix, tmp_cors)
  PITT_txt = cbind(PITT_txt, tmp_res)
}
colnames(PITT_PCAmatrix) = c("Sex", "Age", "Plate",
                              "Neuronal_Inhibitory", "Neuronal_Excitatory", "Oligodendrocyte", "Microglia", "Astrocyte")
rownames(PITT_PCAmatrix) = paste0("PC", rownames(PITT_PCAmatrix))
pheatmap::pheatmap(PITT_PCAmatrix, treeheight_row = 0, treeheight_col = 0, 
                   color = redWhiteGreen(20), breaks = seq(-1, 1, by = 0.1),
                   cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = PITT_txt, number_color = "black")

#####

dirname = "PITT"
output = make_dataset(data_matrix = betas_PITT, metadata = pheno_PITT, dirname = dirname, 
                      db = "PITT", cpgs_intersect = cpgs_intersect)  
data_matrix_AD = as.matrix(output[[1]])
metadata_AD = output[[2]]
remove(output)

##### Make models #####

#Preview to select number of clusters
# df_assignment = wrapper_models(data_matrix_AD, metadata_AD,cluster_vector = c(), dirname = dirname)

#Optimal number of clusters
df_assignment = wrapper_models(data_matrix_AD, metadata_AD,cluster_vector = c(3,3), dirname = dirname)

  
##### SPLSDA output of the identified subtypes ##### 
  rownames(df_assignment) = metadata_AD$Basename
  pheno_PITT = add_subtype_pheno(pheno_PITT, df_assignment)

  metadata_AD = pheno_PITT[pheno_PITT$HC_clusters != "Control",]
  
  diablo_Methyls = list()
  diablo_Methyls[[1]] = mixOmics::splsda(t(data_matrix_AD), metadata_AD$HC_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  diablo_Methyls[[2]] = mixOmics::splsda(t(data_matrix_AD), metadata_AD$KM_clusters, keepX=c(2000,2000,2000,2000,2000,2000), ncomp = 6)
  
  for(dMethyl in seq(1:length(diablo_Methyls))){
    if(dMethyl == 1){
      outdir = paste0("results/", dirname, "/HClust/")
    } else if (dMethyl == 2){
      outdir = paste0("results/", dirname, "/Kmeans/")
    }
    
    plot_mixOmics_all(model = diablo_Methyls[[dMethyl]], outdir = outdir)
  }
  # save.image("3cohort_clustering.Rdata")
  save(betas_PITT, pheno_PITT, betas_ROSMAP, pheno_ROSMAP, betas_UKBBN, pheno_UKBBN,
       file = "new_corr_allcohorts.Rdata")
  
  ##### Random permutations for specificity #####

  df_shuffle = random_permutations(data_matrix_AD, metadata_AD, diablo_Methyls)

  colnames(df_shuffle) = c("HC_accuracy","HC_intersect","HC_prop",
                           "KM_accuracy","KM_intersect","KM_prop",
                           "PAM_accuracy","PAM_intersect","PAM_prop",
                           "Spec_accuracy","Spec_intersect","Spec_prop")
  
  df_shuffle = df_shuffle[-1,]

  output_permutation(df_shuffle, dirname)
    
  save(data_matrix_AD, metadata_AD, diablo_Methyls, df_shuffle, file = paste0("results/", dirname, "/res_clustering.Rdata"))