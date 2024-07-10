##### Sourcing scripts, libraries and data #####
setwd("D:/valentin/main/")

sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}


load("betas/ROSMAP_final_set_27-11-23.Rdata")
load("betas/UKBBN_final_set_27-9-2023.Rdata")
load("betas/PITT_final_set_15-11-23.Rdata")
##### Preparing data #####

cpgs_ROSMAP = rownames(betas_ROSMAP)
cpgs_UKBBN = rownames(betas_UKBBN)
cpgs_PITT = rownames(betas_PITT)
cpgs_intersect = Reduce(intersect, list(cpgs_PITT, cpgs_UKBBN))
cpgs_intersect_ROSMAP = intersect(cpgs_intersect, rownames(betas_ROSMAP))

##### Uniformize data and generate output directory

dirname = "PITT"
output = make_dataset(data_matrix = betas_PITT, metadata = pheno_PITT, dirname = dirname, 
                      db = "PITT", cpgs_intersect = cpgs_intersect)  
data_matrix_AD = as.matrix(output[[1]])
metadata_AD = output[[2]]
remove(output)

##### Make models #####

#1st run to select number of clusters
df_assignment = wrapper_models(data_matrix_AD, metadata_AD,cluster_vector = c(), dirname = dirname)

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