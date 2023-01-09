##### Sourcing scripts and libraries #####
setwd("C:/Users/P70077107/Desktop/PhD/main")

sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}

##### Loading data #####

  BDR_all = read.csv("raw_data/BDR_FULL_PHENOTYPING.csv", sep = ";")
  pheno = read.csv("./BDR_pheno.csv", sep = "\t")
  # load("./processed.Rdata")
  load("./Met_BDR_preproc.Rdata")
  
##### Sample filtering #####
pheno_PFC = pheno[pheno$Basename %in% colnames(mSet_betas),]
# pheno_PFC = pheno_PFC[pheno_PFC$prop >= 0.3,]
table(pheno_PFC$BraakTangle_numeric)
remove(pheno)
gc()
pheno_PFC = pheno_PFC[order(pheno_PFC$Basename),]
mSet_betas = mSet_betas[,order(colnames(mSet_betas))]
##### Batch effect removal #####

mSet_betas = ComBat(dat=mSet_betas, batch=pheno_PFC$Institute, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
mSet_betas = ComBat(dat=mSet_betas, batch=pheno_PFC$Plate, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

##### Celltype correction #####
resid <- function(row, age, sex, prop){
  fit <- try(
    lm( row ~ age + sex + prop),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
}
##### Adjusting for array sensitivity #####

rowbetas = rowMedians(mSet_betas)

mSet_betas = mSet_betas[rowbetas > 0.1 & rowbetas < 0.9,]
colbetas = colnames(mSet_betas)
mSet_betas<- {
  sex    	<- 	as.factor(pheno_PFC$Gender)
  age    	<- 	pheno_PFC$Age
  prop		  <-  as.numeric(pheno_PFC$prop)
  t(apply(mSet_betas, 1, resid, age, sex, prop))
}

colnames(mSet_betas) = colbetas
##### MAD/Quantile filtering #####

mSet_betas = mSet_betas[!is.na(rownames(mSet_betas)),]

df_values = data.frame(matrix(ncol = 2, nrow = nrow(mSet_betas)))
### Use 80-90% quantile absolute deviation ?
for(i in 1:nrow(mSet_betas)){
  tmp_cpg = mSet_betas[i,]
  value = c(mad(tmp_cpg), mad(tmp_cpg, center = quantile(tmp_cpg, probs = 0.9)))
  df_values[i,] = value
}

rownames(df_values) = rownames(mSet_betas)
colnames(df_values) = c("mad","90pc")

# hist(df_values$mad)
# hist(df_values$"90pc")

df_values$"90pc" = df_values$"90pc"[order(df_values$"90pc", decreasing = T)]
nfeatures = round(0.2*length(df_values$"90pc"))
df_values = df_values[1:nfeatures,]
# hist(df_values$"90pc")

mSet_betas = mSet_betas[rownames(mSet_betas) %in% rownames(df_values),]
gc()

# ##### Correlation filtering #####
# 
# cpg_list = list()
# betas_to_cor = betas
# select_col = rownames(betas_to_cor)[1]
# tmp_cpg_list = select_col
# to_filter = c()
# 
# while (nrow(betas_to_cor) > 1){
# 
#   adjacency = cor(t(betas_to_cor), betas_to_cor[rownames(betas_to_cor) == select_col])
#   
#   colnames(adjacency) = select_col
#   head(adjacency)
#   
#   most_cor = which(adjacency == max(adjacency[adjacency < 0.9999,1]))
#   adjacency[most_cor,]
#   
#   if(adjacency[most_cor,] > 0.8){
#     
#     cor_names = names(adjacency[adjacency > 0.8,])
#     cor_betas = betas_to_cor[rownames(betas_to_cor) %in% cor_names,]
#     
#     adjacency = cor(t(cor_betas))
#     
#     index = which(adjacency == max(adjacency[adjacency <  0.9999]), arr.ind = T) 
# 
#     best_match = df_values[rownames(df_values) %in% rownames(index),]
#     most_var = rownames(best_match[best_match$`90pc` == max(best_match$`90pc`),])
#     
#     to_filter = c(to_filter, cor_names[!cor_names == most_var])
#     
#     betas_to_cor = betas_to_cor[!rownames(betas_to_cor) %in% cor_names,]
#     select_col = rownames(betas_to_cor)[1]
#     tmp_cpg_list = cor_names
#   } else {
#     cpg_list[[length(cpg_list) + 1]] = tmp_cpg_list
#     
#     betas_to_cor = betas_to_cor[rownames(betas_to_cor) != select_col,]
#     select_col = rownames(betas_to_cor)[1]
#     tmp_cpg_list = select_col
#   }
# }


##### Detailing metadata #####

  BDR_filter = BDR_all[BDR_all$BBNId %in% pheno_PFC$BBNId,]
  
  BDR_metadata_filter = pheno_PFC[pheno_PFC$BBNId %in% BDR_filter$BBNId,]
  
  BDR_filter = BDR_filter[order(BDR_filter$BBNId),]
  BDR_metadata_filter = BDR_metadata_filter[order(BDR_metadata_filter$BBNId),]
  
  BDR_metadata_filter$Braak_tangle = BDR_metadata_filter$BraakTangle_numeric
  BDR_metadata_filter$Braak_LB = BDR_filter$Braak.LB.stage
  BDR_metadata_filter$Cerad = BDR_filter$Cerad
  BDR_metadata_filter$CDR = BDR_filter$Final.CDR
  BDR_metadata_filter$Gender = BDR_filter$Gender
  BDR_metadata_filter$Age = BDR_filter$Age
  BDR_metadata_filter$Cause_of_Death = BDR_filter$Cause.Of.Death
  BDR_metadata_filter$MMSE = BDR_filter$Final.MMSE
  BDR_metadata_filter$APOE = BDR_filter$Polymorphisms
  BDR_metadata_filter$clinical = BDR_filter$Clinical.Diagnoses
  BDR_metadata_filter$clinical_plus = BDR_filter$Clin.Diagnosis.Details
  BDR_metadata_filter$pathology = BDR_filter$Pathology
  
  BDR_metadata_filter$clinical
  
  
  BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD no plaques "] = 0
  BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD sparse neuritic plaques "] = 1
  BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD moderate density of neuritic plaques "] = 2
  BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD high density of neuritic plaques "] = 3
  BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == ""] = NA
  
  BDR_metadata_filter$age_conc = NA ### Making age into simpler intervals
  BDR_metadata_filter[BDR_metadata_filter$Age < summary(BDR_metadata_filter$Age)[2],]$age_conc = 1
  BDR_metadata_filter[BDR_metadata_filter$Age >= summary(BDR_metadata_filter$Age)[2] & 
                        BDR_metadata_filter$Age < summary(BDR_metadata_filter$Age)[3],]$age_conc = 2
  BDR_metadata_filter[BDR_metadata_filter$Age >= summary(BDR_metadata_filter$Age)[3] & 
                        BDR_metadata_filter$Age < summary(BDR_metadata_filter$Age)[5],]$age_conc = 3
  BDR_metadata_filter[BDR_metadata_filter$Age >= summary(BDR_metadata_filter$Age)[5],]$age_conc = 4
  
  BDR_metadata_filter$dementia = NA
  BDR_metadata_filter$diag = NA
  BDR_metadata_filter$CDR[BDR_metadata_filter$CDR == -9] = NA
  BDR_metadata_filter$CDR[BDR_metadata_filter$CDR == 18] = NA
  BDR_metadata_filter$CDR[BDR_metadata_filter$CDR == 6] = NA
  BDR_metadata_filter$dementia[BDR_metadata_filter$MMSE < 24] = 1
  BDR_metadata_filter$dementia[BDR_metadata_filter$CDR >= 1] = 1
  BDR_metadata_filter$dementia[BDR_metadata_filter$CDR < 1] = 0
  BDR_metadata_filter$dementia[is.na(BDR_metadata_filter$dementia)] = 2
  BDR_metadata_filter$diag[BDR_metadata_filter$Cerad <= 1 & BDR_metadata_filter$Braak_tangle <= 2 & BDR_metadata_filter$dementia == 0] = "Control" 
  BDR_metadata_filter$diag[BDR_metadata_filter$Cerad == 0 & BDR_metadata_filter$Braak_tangle <= 3 & BDR_metadata_filter$dementia == 0] = "Control" 
  BDR_metadata_filter$diag[BDR_metadata_filter$Cerad >= 1 & BDR_metadata_filter$Braak_tangle >= 3 & BDR_metadata_filter$dementia == 0] = "AsymAD"
  BDR_metadata_filter$diag[BDR_metadata_filter$Cerad >= 2 & BDR_metadata_filter$Braak_tangle >= 3 & BDR_metadata_filter$dementia == 1] = "AD" 
  BDR_metadata_filter$diag[is.na(BDR_metadata_filter$diag)] = "Unassigned"
  
  
##### Curating dataset #####
  BDR_ok = BDR_metadata_filter[BDR_metadata_filter$diag != "Unassigned",]
  
  # Clinical AD
  BDR_unassigned = BDR_metadata_filter[BDR_metadata_filter$diag == "Unassigned",]
  BDR_unassigned[grep("Alzheimer", BDR_unassigned$clinical),]$clinical
  BDR_unassigned[grep("Alzheimer", BDR_unassigned$clinical),]$diag = "AD"
  BDR_ok = rbind(BDR_ok, BDR_unassigned[BDR_unassigned$diag != "Unassigned",])
  BDR_unassigned = BDR_unassigned[BDR_unassigned$diag == "Unassigned",]
  
  
  # Clinical dementia
  BDR_unassigned[grep("dementia", BDR_unassigned$clinical),]$clinical
  BDR_unassigned[grep("dementia", BDR_unassigned$clinical),]$dementia = 1
  
  BDR_unassigned$diag[BDR_unassigned$Cerad <= 1 & BDR_unassigned$Braak_tangle <= 2 & BDR_unassigned$dementia == 0] = "Control" 
  BDR_unassigned$diag[BDR_unassigned$Cerad == 0 & BDR_unassigned$Braak_tangle <= 3 & BDR_unassigned$dementia == 0] = "Control" 
  BDR_unassigned$diag[BDR_unassigned$Cerad >= 1 & BDR_unassigned$Braak_tangle >= 3 & BDR_unassigned$dementia == 0] = "AsymAD"
  BDR_unassigned$diag[BDR_unassigned$Cerad >= 2 & BDR_unassigned$Braak_tangle >= 3 & BDR_unassigned$dementia == 1] = "AD" 
  table(BDR_unassigned$diag)
  BDR_ok = rbind(BDR_ok, BDR_unassigned[BDR_unassigned$diag != "Unassigned",])
  BDR_unassigned = BDR_unassigned[BDR_unassigned$diag == "Unassigned",]
  
  # Braak tangle and Cerad NA's can't be fixed - no data
  # No abnormalities = controls
  BDR_unassigned[grep("No abnormality detected", BDR_unassigned$clinical),]$clinical
  BDR_unassigned[grep("No abnormality detected", BDR_unassigned$clinical),]$diag = "Control"
  BDR_ok = rbind(BDR_ok, BDR_unassigned[BDR_unassigned$diag != "Unassigned",])
  BDR_unassigned = BDR_unassigned[BDR_unassigned$diag == "Unassigned",]
  
  # Rest
  BDR_unassigned$clinical
  BDR_unassigned$clinical_plus
  BDR_unassigned$pathology
  BDR_unassigned$Braak_LB
  
  # Final = We exclude BDR_unassigned 88 remaining samples because not enough information OR other types of Dementia
  # Final = We exclude the 25 AsymAD samples and age <65
  
  
  BDR_ok_age = BDR_ok[BDR_ok$Age >= 65,]
  table(BDR_ok_age$diag)
  BDR_ok_ADCtrl = BDR_ok_age[BDR_ok_age$diag != "AsymAD",]
  
  table(BDR_ok_ADCtrl$diag)
  table(BDR_ok_ADCtrl$dementia)

  table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$MMSE)
  table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$CDR)
  table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$Braak_tangle)
  table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$Cerad)
  
  table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD" & BDR_ok_ADCtrl$CDR < 1,]$MMSE)

  
##### Curating Controls #####
  BDR_ok_ADCtrl = BDR_ok_ADCtrl[order(BDR_ok_ADCtrl$BBNId),]
  BDR_ok_ADCtrl$BBNId
  BDR_filter = BDR_filter[order(BDR_filter$BBNId),]
  BDR_filter$BBNId
  BDR_ok_ADCtrl$Thal_phase = BDR_filter$Thal.stage
  
  BDR_ok_ADCtrl$Thal_number = NA
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == "Thal phase 0 "] = 0
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == "Thal phase 1 "] = 1
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == "Thal phase 2 "] = 2
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == "Thal phase 3 "] = 3
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == "Thal phase 4 "] = 4
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == "Thal phase 5 "] = 5
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == "Thal phase not assessed "] = NA
  BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal_phase == ""] = NA
  table(BDR_ok_ADCtrl$Thal_phase)
  table(BDR_ok_ADCtrl$Thal_number)
  
  BDR_Ctrl = BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "Control",]
  BDR_Ctrl$new_diag = NA
  BDR_Ctrl$Thal_phase
  
  
  #Braak stage 0-2
  table(BDR_Ctrl$BraakTangle_numeric)
  BDR_Ctrl_BTnonOK = BDR_Ctrl[which(BDR_Ctrl$Braak_tangle >= 4),]
  
  # Argyrophilic
  grep("+Ndgen.NonADTau.AGD", BDR_Ctrl$pathology, fixed = TRUE)
  BDR_Ctrl_Argnonok = BDR_Ctrl[grep("+Ndgen.NonADTau.AGD", BDR_Ctrl$pathology, fixed = TRUE),]
  BDR_Ctrl_Argnonok$pathology  

  # Thal and Cerad 0
  table(BDR_Ctrl$Cerad, BDR_Ctrl$Thal_number)
  BDR_Ctrl_amynonok = BDR_Ctrl[which(BDR_Ctrl$Cerad > 0 | BDR_Ctrl$Thal_number > 0),]

  # TDP43
  grep("+Ndgen.TDP.TDPaccomp", BDR_Ctrl$pathology, fixed = TRUE)
  BDR_Ctrl_TDP43 = BDR_Ctrl[grep("+Ndgen.TDP.TDPaccomp", BDR_Ctrl$pathology, fixed = TRUE),]
  BDR_Ctrl_TDP43$pathology
  
  # Alpha synuclein - no info ??
  grep("a-", BDR_Ctrl$pathology, fixed = TRUE)
  BDR_Ctrl_alpha = BDR_Ctrl[grep("a-", BDR_Ctrl$pathology, fixed = TRUE),]
  BDR_Ctrl_alpha$pathology
  
  # Cerebrovascular
  
  grep("vascular", BDR_Ctrl$clinical)
  grep("vascular", BDR_Ctrl$clinical_plus)
  grep("vascular", BDR_Ctrl$pathology)
  
  final_cereb = c(grep("vascular", BDR_Ctrl$clinical, fixed = TRUE),
                  grep("vascular", BDR_Ctrl$clinical_plus, fixed = TRUE),
                  grep("+CVD.Vasc.Athero.Sevathero", BDR_Ctrl$pathology, fixed = TRUE))
  BDR_Ctrl_cereb = BDR_Ctrl[unique(final_cereb),]
  BDR_Ctrl_cereb$clinical
  BDR_Ctrl_cereb$clinical_plus
  BDR_Ctrl_cereb$pathology
  
  # Leukemia - None ?
  grep("LEUKEMIA", BDR_Ctrl$Cause_of_Death)
  
  # Neuropathologic correlate - None ?
  grep("psychia", BDR_Ctrl$clinical_plus)
  
  
  # Final assessment
  final_ctrl_nonok = c(grep("vascular", BDR_Ctrl$clinical, fixed = TRUE),
                       grep("vascular", BDR_Ctrl$clinical_plus, fixed = TRUE),
                       grep("+CVD.Vasc.Athero.Sevathero", BDR_Ctrl$pathology, fixed = TRUE),
                       grep("+Ndgen.TDP.TDPaccomp", BDR_Ctrl$pathology, fixed = TRUE),
                       which(BDR_Ctrl$Cerad > 1 | BDR_Ctrl$Thal_number > 2),
                       grep("+Ndgen.NonADTau.AGD", BDR_Ctrl$pathology, fixed = TRUE),
                       which(BDR_Ctrl$Braak_tangle >= 3))
  final_ctrl_nonok = order(unique(final_ctrl_nonok))
  
  
  BDR_remove = BDR_Ctrl[-final_ctrl_nonok,]
  BDR_remove$Basename
##### Train / Test models #####
  #### Normal sampling #### 
  # training_ID = splitstackshape::stratified(BDR_metadata_AD, c("age_conc","Gender"),size = 0.6)
  # training_ID = training_ID$Basename
  # test_ID = BDR_metadata_AD[!BDR_metadata_AD$Basename %in% training_ID,]$Basename
  #### Kstone sampling ####
  # sample_select = kenStone(X = t(mSet_betas),
  #                          k = round(0.7*ncol(mSet_betas)),
  #                          pc = 20,
  #                          .center = TRUE,
  #                          .scale = TRUE)
  # 
  # training = mSet_betas[,sample_select$model]
  # testing = mSet_betas[,sample_select$test]
  # BDR_metadata_train = BDR_ok_ADCtrl[sample_select$model,]
  # BDR_metadata_test = BDR_ok_ADCtrl[sample_select$test,]
  #### Split AD Ctrl ####
  
  BDR_ok_ADCtrl = BDR_ok_ADCtrl[order(BDR_ok_ADCtrl$Basename),]
  mSet_betas = mSet_betas[,order(colnames(mSet_betas))]
  head(colnames(mSet_betas))
  head(BDR_ok_ADCtrl$Basename)
  table(BDR_ok_ADCtrl$diag)
  
  BDR_ok_ADCtrl$age_conc = NA ### Making age into simpler intervals
  BDR_ok_ADCtrl[BDR_ok_ADCtrl$Age < summary(BDR_ok_ADCtrl$Age)[2],]$age_conc = 1
  BDR_ok_ADCtrl[BDR_ok_ADCtrl$Age >= summary(BDR_ok_ADCtrl$Age)[2] & 
                        BDR_ok_ADCtrl$Age < summary(BDR_ok_ADCtrl$Age)[3],]$age_conc = 2
  BDR_ok_ADCtrl[BDR_ok_ADCtrl$Age >= summary(BDR_ok_ADCtrl$Age)[3] & 
                        BDR_ok_ADCtrl$Age < summary(BDR_ok_ADCtrl$Age)[5],]$age_conc = 3
  BDR_ok_ADCtrl[BDR_ok_ADCtrl$Age >= summary(BDR_ok_ADCtrl$Age)[5],]$age_conc = 4
  
  BDR_metadata_AD = BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]
  BDR_metadata_Ctrl = BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "Control",]
  mSet_betas_AD = mSet_betas[,colnames(mSet_betas) %in% BDR_metadata_AD$Basename]
  mSet_betas_Ctrl = mSet_betas[,colnames(mSet_betas) %in% BDR_metadata_Ctrl$Basename]
  
  set.seed(3)  
  
  # gc()
  # affM_Methyl_all = aff_matrix_calc(mSet_betas)
  # affM_Methyl_noCor = aff_matrix_calc(mSet_betas_noCor)
##### Building datasets #####

  list_datasets = list()
  cols_df = ""
  
  list_datasets[[length(list_datasets)+1]] = mSet_betas_AD
  names(list_datasets)[length(list_datasets)] = "Methyl"
  
  # list_aff_matrix[[length(list_aff_matrix)+1]] = affM_Methyl_all
  # names(list_aff_matrix)[length(list_aff_matrix)] = "Methyl"
  dirname = ""
  dirname = paste(dirname, "Methyl", sep = "_")
  cols_df = paste0(cols_df, "M")
  
  dirname = gsub("^_", "", dirname)
  dir.create("results/")
  dir.create(paste0("results/", dirname))
  dir.create(paste0("results/", dirname, "/Spectrum"))
  dir.create(paste0("results/", dirname, "/HClust"))
  dir.create(paste0("results/", dirname, "/Kmeans"))
  dir.create(paste0("results/", dirname, "/PAM"))

##### Building models #####
  
  multi_Spectrum = Spectrum_modelization(list_datasets, spectrum_method = 2, 
                                         cluster_alg = "GMM", num_clusters = 2)
  
  png(paste0("results/", dirname, "/Spectrum/eigenvector_analysis.png"), width = 800, height = 600)
  plot(multi_Spectrum$eigenvector_analysis$K, multi_Spectrum$eigensystem$values[1:11])
  dev.off()
  
  multi_Spectrum = Spectrum_modelization(list_datasets, spectrum_method = 3, 
                                         cluster_alg = "GMM", num_clusters = 4)
  plot_model(multi_Spectrum, "Spectrum", colnames(mSet_betas_AD), paste0("results/", dirname, "/Spectrum/"))
  df_assignment = data.frame("Spectrum" = multi_Spectrum$assignments,
                             row.names = names(multi_Spectrum$assignments))  

  multi_HC = CC_modelization(list_datasets[[1]],
                             title_model = paste0("results/", dirname, "/HClust/"), 
                             cluster_alg = "hc", distance = "pearson")
  plot_model(multi_HC, "CC", colnames(mSet_betas_AD), paste0("results/", dirname, "/HClust/"),
             numCC_clusters = 3)
  df_assignment$HC = multi_HC[[3]]$consensusClass
  
  multi_KM = CC_modelization(list_datasets[[1]],
                             title_model = paste0("results/", dirname, "/Kmeans/"), 
                             cluster_alg = "km", distance = "euclidean")
  plot_model(multi_KM, "CC", colnames(mSet_betas_AD), paste0("results/", dirname, "/Kmeans/"),
             numCC_clusters = 4)
  df_assignment$KM = multi_KM[[4]]$consensusClass
  
  multi_PAM = CC_modelization(list_datasets[[1]],
                             title_model = paste0("results/", dirname, "/PAM/"), 
                             cluster_alg = "pam", distance = "pearson")
  plot_model(multi_PAM, "CC", colnames(mSet_betas_AD), paste0("results/", dirname, "/PAM/"),
             numCC_clusters = 4)
  df_assignment$PAM = multi_PAM[[4]]$consensusClass
  
  
  df_NMI_models = calcNMI_all(df_assignment)
  models_tables(df_assignment, df_NMI_models, paste0("results/", dirname, "/"))
  colnames(df_assignment) = paste(cols_df, colnames(df_assignment), sep = "_")
  
##### mixOmics #####
  
  BDR_metadata_AD$HC_clusters = multi_HC[[3]]$consensusClass
  BDR_metadata_AD$KM_clusters = multi_KM[[4]]$consensusClass
  BDR_metadata_AD$PAM_clusters = multi_PAM[[4]]$consensusClass
  BDR_metadata_AD$Spec_clusters = multi_Spectrum$assignments
  
  
  diablo_Methyls = list()
  diablo_Methyls[[1]] = mixOmics::splsda(t(mSet_betas_AD), BDR_metadata_AD$HC_clusters, keepX=c(200,200,200,200,200,200), ncomp = 6)
  diablo_Methyls[[2]] = mixOmics::splsda(t(mSet_betas_AD), BDR_metadata_AD$KM_clusters, keepX=c(200,200,200,200), ncomp = 4)
  diablo_Methyls[[3]] = mixOmics::splsda(t(mSet_betas_AD), BDR_metadata_AD$PAM_clusters, keepX=c(200,200,200,200), ncomp = 4)
  diablo_Methyls[[4]] = mixOmics::splsda(t(mSet_betas_AD), BDR_metadata_AD$Spec_clusters, keepX=c(200,200,200,200), ncomp = 4)
  
  for(dMethyl in seq(1:length(diablo_Methyls))){
    if(dMethyl == 1){
      outdir = paste0("results/", dirname, "/HClust/")
    } else if (dMethyl == 2){
      outdir = paste0("results/", dirname, "/Kmeans/")
    } else if (dMethyl == 3){
      outdir = paste0("results/", dirname, "/PAM/")
    } else if (dMethyl == 4){
      outdir = paste0("results/", dirname, "/Spectrum/")
    }
    
    plot_mixOmics_all(model = diablo_Methyls[[dMethyl]], outdir = outdir)
  }
  
  save.image("./current_status.Rdata")
  
##### Shuffle clusters #####  
  
  HC_roc = auroc(diablo_Methyls[[1]])
  KM_roc = auroc(diablo_Methyls[[2]])
  PAM_roc = auroc(diablo_Methyls[[3]])
  Spec_roc = auroc(diablo_Methyls[[4]])
  final_aucs = c(mean(HC_roc$Comp4[,1]),
                 mean(KM_roc$Comp4[,1]),
                 mean(PAM_roc$Comp4[,1]),
                 mean(Spec_roc$Comp4[,1]))
  
  
  HC_var = plotVar(diablo_Methyls[[1]])
  HC_features = HC_var$names
  HC_var = plotVar(diablo_Methyls[[1]])
  HC_features = c(HC_features, HC_var$names)
  
  KM_var = plotVar(diablo_Methyls[[2]])
  KM_features = KM_var$names
  KM_var = plotVar(diablo_Methyls[[2]])
  KM_features = c(KM_features, KM_var$names)
  
  PAM_var = plotVar(diablo_Methyls[[3]])
  PAM_features = PAM_var$names
  PAM_var = plotVar(diablo_Methyls[[3]])
  PAM_features = c(PAM_features, PAM_var$names)
  
  Spec_var = plotVar(diablo_Methyls[[4]])
  Spec_features = Spec_var$names
  Spec_var = plotVar(diablo_Methyls[[4]])
  Spec_features = c(Spec_features, Spec_var$names)
  
  df_shuffle = data.frame(matrix(nrow = 100, ncol = 8))
  colnames(df_shuffle) = c("HC_auc","HC_ratio",
                           "KM_auc","KM_ratio",
                           "PAM_auc","PAM_ratio",
                           "Spec_auc","Spec_ratio")
  
  for(i in seq(1:100)){
    HC_shuffle = sample(BDR_metadata_AD$HC_clusters)
    KM_shuffle = sample(BDR_metadata_AD$KM_clusters)
    PAM_shuffle = sample(BDR_metadata_AD$PAM_clusters)
    Spec_shuffle = sample(BDR_metadata_AD$Spec_clusters)
    
    diablo_shuffles = list()
    diablo_shuffles[[1]] = mixOmics::splsda(t(mSet_betas_AD), HC_shuffle, keepX=c(200,200,200,200), ncomp = 4)
    diablo_shuffles[[2]] = mixOmics::splsda(t(mSet_betas_AD), KM_shuffle, keepX=c(200,200,200,200), ncomp = 4)
    diablo_shuffles[[3]] = mixOmics::splsda(t(mSet_betas_AD), PAM_shuffle, keepX=c(200,200,200,200), ncomp = 4)
    diablo_shuffles[[4]] = mixOmics::splsda(t(mSet_betas_AD), Spec_shuffle, keepX=c(200,200,200,200), ncomp = 4)
    
    HC_roc = auroc(diablo_shuffles[[1]])
    KM_roc = auroc(diablo_shuffles[[2]])
    PAM_roc = auroc(diablo_shuffles[[3]])
    Spec_roc = auroc(diablo_shuffles[[4]])
    
    HC_shuffle_var = plotVar(diablo_shuffles[[1]])
    HC_shuffle_features = HC_shuffle_var$names
    HC_shuffle_var = plotVar(diablo_shuffles[[1]])
    HC_shuffle_features = c(HC_shuffle_features, HC_shuffle_var$names)
    
    KM_shuffle_var = plotVar(diablo_shuffles[[2]])
    KM_shuffle_features = KM_shuffle_var$names
    KM_shuffle_var = plotVar(diablo_shuffles[[2]])
    KM_shuffle_features = c(KM_shuffle_features, KM_shuffle_var$names)
    
    PAM_shuffle_var = plotVar(diablo_shuffles[[3]])
    PAM_shuffle_features = PAM_shuffle_var$names
    PAM_shuffle_var = plotVar(diablo_shuffles[[3]])
    PAM_shuffle_features = c(PAM_shuffle_features, PAM_shuffle_var$names)
    
    Spec_shuffle_var = plotVar(diablo_shuffles[[4]])
    Spec_shuffle_features = Spec_shuffle_var$names
    Spec_shuffle_var = plotVar(diablo_shuffles[[4]])
    Spec_shuffle_features = c(Spec_shuffle_features, Spec_shuffle_var$names)
    
    HC_inter = intersect(HC_features, HC_shuffle_features)
    KM_inter = intersect(KM_features, KM_shuffle_features)
    PAM_inter = intersect(PAM_features, PAM_shuffle_features)
    Spec_inter = intersect(Spec_features, Spec_shuffle_features)
    
    shuffle_res = c(mean(HC_roc$Comp4[,1]),
                    length(HC_inter)/length(HC_features),
                    mean(KM_roc$Comp4[,1]),
                    length(KM_inter)/length(KM_features),
                    mean(PAM_roc$Comp4[,1]),
                    length(PAM_inter)/length(PAM_features),
                    mean(Spec_roc$Comp4[,1]),
                    length(Spec_inter)/length(Spec_features))
    df_shuffle[i,] = shuffle_res
  }
  
  
  colnames(df_shuffle) = c("HC_accuracy","HC_intersect",
                           "KM_accuracy","KM_intersect",
                           "PAM_accuracy","PAM_intersect",
                           "Spec_accuracy","Spec_intersect")
  
  ggplot(data = df_shuffle, aes(x = HC_accuracy, y = HC_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave("./results/Methyl/HClust/HC_acc.jpeg")
  ggplot(data = df_shuffle, aes(x = KM_accuracy, y = KM_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave("./results/Methyl/Kmeans/KM_acc.jpeg")
  ggplot(data = df_shuffle, aes(x = PAM_accuracy, y = PAM_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave("./results/Methyl/PAM//PAM_acc.jpeg")
  ggplot(data = df_shuffle, aes(x = Spec_accuracy, y = Spec_intersect)) +
    geom_point() +
    xlim(0,1) + ylim(0,1) +
    theme_minimal()
  ggsave("./results/Methyl/Spectrum/Spec_acc.jpeg")
  
  
  plot(df_shuffle[,1], df_shuffle[,2])
  plot(df_shuffle[,3], df_shuffle[,4])
  plot(df_shuffle[,5], df_shuffle[,6])
  plot(df_shuffle[,7], df_shuffle[,8])
  
##### Cross-method clustering #####

  table(BDR_metadata_AD$HC_clusters, BDR_metadata_AD$KM_clusters)
  table(BDR_metadata_AD$HC_clusters, BDR_metadata_AD$PAM_clusters)
  table(BDR_metadata_AD$HC_clusters, BDR_metadata_AD$Spec_clusters)
  table(BDR_metadata_AD$KM_clusters, BDR_metadata_AD$PAM_clusters)
  table(BDR_metadata_AD$KM_clusters, BDR_metadata_AD$Spec_clusters)
  table(BDR_metadata_AD$PAM_clusters, BDR_metadata_AD$Spec_clusters)
  
  BDR_metadata_AD$all_clust = NA
  n_profile = 0
  for(i in seq(1:nrow(BDR_metadata_AD))){
    tmp_profile = paste(c(BDR_metadata_AD[i,]$HC_clusters,
                        BDR_metadata_AD[i,]$KM_clusters,
                        BDR_metadata_AD[i,]$PAM_clusters,
                        BDR_metadata_AD[i,]$Spec_clusters),
                        collapse = "-")
    BDR_metadata_AD[i,]$all_clust = tmp_profile
  }
  
  table(BDR_metadata_AD$all_clust)[order(table(BDR_metadata_AD$all_clust))]

  table(BDR_metadata_AD$all_clust) %>% 
    as.data.frame() %>% 
    arrange(Var1)
  
    
##### Clinical NMI #####
  
  
  df_assignment_train = df_assignment[rownames(df_assignment) %in% BDR_metadata_AD$Basename,]
  
  df_NMI_clinical = calcNMI_clinical(df_assignment = df_assignment_train, BDR_metadata = BDR_metadata_AD)
  df_NMI_clinical[is.na(df_NMI_clinical)] = ""
  write.table(df_NMI_clinical, file = "./results/Methyl/clinical_NMI.csv", 
              sep = "\t", quote = F, row.names = F, col.names = T)
  ##### Optimise parameters #####
  
  # perf_splsda <- perf(diablo_Methyls[[1]], validation = "Mfold", 
  #                           folds = 5, nrepeat = 10, # use repeated cross-validation
  #                           progressBar = FALSE, auc = TRUE) # include AUC values
  # 
  # colnames(training)
  # tuned_splsda <- tune.splsda(X = diablo_Methyls[[1]]$X, Y = diablo_Methyls[[1]]$Y, 
  #                                  ncomp = 4, # calculate for first 4 components
  #                                  validation = 'Mfold',
  #                                  already.tested.X = c(100,100),
  #                                  folds = 4, nrepeat = 10, # use repeated cross-validation
  #                                  dist = 'max.dist', # use max.dist measure
  #                                  measure = "BER", # use balanced error rate of dist measure
  #                                  test.keepX = c(1:10,  seq(20, 300, 10)),
  #                                  scale = FALSE)
  # plot(tuned_splsda, col = color.jet(4)) # plot output of variable number tuning
  # tuned_splsda$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
  # tuned_splsda$choice.keepX # what are the optimal values of variables according to tune.splsda()
  
  ##### Prediction onto testing #####
  
  testing = cbind(testing_AD, testing_Ctrl)
  pheno_test = rbind(BDR_metadata_test_AD, BDR_metadata_test_Ctrl)
  
  HC_pred = predict(diablo_Methyls[[1]], newdata = t(testing))
  KM_pred = predict(diablo_Methyls[[2]], newdata = t(testing))
  PAM_pred = predict(diablo_Methyls[[3]], newdata = t(testing))
  Spec_pred = predict(diablo_Methyls[[4]], newdata = t(testing))
  
  HC_roc = auroc(diablo_Methyls[[1]], newdata = t(testing), 
                 outcome.test = HC_pred$MajorityVote$mahalanobis.dist[,6])
  table(HC_pred$MajorityVote$mahalanobis.dist[,6], pheno_test$diag)
  
  ##### Sanity check #####
  
  diablo_Tests = list()
  diablo_Tests[[1]] = mixOmics::splsda(t(testing), HC_pred, keepX=c(200,200,200,200), ncomp = 4)
  diablo_Tests[[2]] = mixOmics::splsda(t(testing), KM_pred, keepX=c(200,200,200,200), ncomp = 4)
  diablo_Tests[[3]] = mixOmics::splsda(t(testing), PAM_pred, keepX=c(200,200,200,200), ncomp = 4)
  diablo_Tests[[4]] = mixOmics::splsda(t(testing), Spec_pred, keepX=c(200,200,200,200), ncomp = 4)
  
  for(dMethyl in seq(1:length(diablo_Tests))){
    if(dMethyl == 1){
      outdir = paste0("results/", dirname, "/HClust_test/")
      dir.create(outdir)
    } else if (dMethyl == 2){
      outdir = paste0("results/", dirname, "/Kmeans_test/")
      dir.create(outdir)
    } else if (dMethyl == 3){
      outdir = paste0("results/", dirname, "/PAM_test/")
      dir.create(outdir)
    } else if (dMethyl == 4){
      outdir = paste0("results/", dirname, "/Spectrum_test/")
      dir.create(outdir)
    }
    
    plot_mixOmics_all(model = diablo_Tests[[dMethyl]], outdir = outdir)
  }
  
  
  
##### Projecting to ROSMAP #####
  
  load("ROSMAP_data/ROSMAP.711.beta.pheno.Rdata")

  ROSMAP_data = betas[rownames(betas) %in% colnames(training),]
  diablo_pred = predict(diablo_Methyls[[3]], newdata = ROSMAP_data)
  
################################## Testing ##################################
  # dir.create(paste0("results/", dirname, "/iClusterBayes"))
  # multi_iCB = ICB_modelization(list_datasets,
  #                              data_distrib = rep("gaussian", length(list_datasets)), K = 3)
  # # png(paste0("results/", dirname, "/iClusterBayes/feature_importance.png"), width = 2000, height = 1600)
  # # plotHMBayes(multi_iCB, list_datasets, type = rep("gaussian", length(list_datasets)),
  # #             sample.order = multi_iCB$clusters, width = 20)
  # # dev.off()
  # 
  # plot_model(multi_iCB, "iCB", colnames(mSet_betas), paste0("results/", dirname, "/iClusterBayes/"))
  # 
  # df_assignment$iCB = multi_iCB$clusters
  # 
  # dir.create(paste0("results/", dirname, "/ANF"))
  # 
  # multi_ANF = ANF_modelization(list_aff_matrix, K = 20, num_clusters = c(2:6))
  # 
  # eigen(multi_ANF[[4]][[1]])
  # 
  # CDF_model = list()
  # for (i in seq(1:length(multi_ANF))){
  #   CDF_model[[i]] = multi_ANF[[i]][[1]]
  # }
  # 
  # 
  # plot_model(multi_ANF, "ANF", colnames(mSet_betas), paste0("results/", dirname, "/ANF/"))
  # 
  # df_assignment$ANF = multi_ANF[[2]]
  # 
  # if(length(list_aff_matrix) >= 2){
  #   dir.create(paste0("results/", dirname, "/SNF"))
  #   multi_SNF = SNF_modelization(list_aff_matrix, K = 20, t = 20, num_clusters = 4)
  #   plot_model(multi_SNF, "SNF", colnames(mSet_betas_cor), paste0("results/", dirname, "/SNF/"))
  #   
  #   df_assignment$SNF = multi_SNF[[2]]
  #   
  #   dir.create(paste0("results/", dirname, "/SNF.CC"))
  #   multi_SNF.CC = SNF.CC_modelization(datasets = list_datasets, num_clusters = 4)
  #   plot_model(multi_SNF.CC, "SNF.CC", colnames(mSet_betas_cor), paste0("results/", dirname, "/SNF.CC/"))
  #   
  # }
  
# multi_JIVE = JIVE_modelization(datasets = list("RNA" = RNA_gene_norm,
#                                                "miRNA" = miRNA_norm, 
#                                                "circRNA" = circRNA_norm, 
#                                                "proteomics" = proteomics_norm))
# 
# JIVE_analysis(jive_model = multi_JIVE)
# 
# 
# X_diablo = list("RNA" = t(RNA_gene_norm),
#                 "miRNA" = t(miRNA_norm), 
#                 "circRNA" = t(circRNA_norm), 
#                 "proteomics" = t(proteomics_norm))
# Y_diablo = BDR_metadata$Case.Control
# keepX_diablo = list("RNA" = c(50, 50),
#                     "miRNA" = c(50, 50),
#                     "circRNA" = c(50, 50),
#                     "proteomics" = c(50, 50))
# multi_DIABLO = block.splsda(X = X_diablo,
#                             Y = Y_diablo,
#                             keepX = keepX_diablo)
# 
# plotIndiv(multi_DIABLO) ## sample plot
# plotVar(multi_DIABLO, var.names = F, legend = TRUE, pch = c(20,20,20,20)) ## variable plot
# 
# plotIndiv(multi_DIABLO, 
#           ind.names = FALSE, 
#           legend=TRUE, cex=c(1,2,3), pch = c(15,18,20),
#           title = 'BRCA with DIABLO')
#  
# Y = multi_DIABLO$Y
# plotDiablo(multi_DIABLO, ncomp = 1)
# 
# circosPlot(multi_DIABLO, cutoff=0.8)
# 
# png("test.png", width = 3000, height = 1600)
# lol = cimDiablo(multi_DIABLO, color.blocks = c('darkorchid', 'brown1', 'lightgreen', "dodgerblue2"), margins = c(5, 22), comp = 1, legend.position = "right")
# dev.off()
# 
# plotLoadings(multi_DIABLO, comp = 1, contrib = "max")
# 
# network(multi_DIABLO, blocks = c(1,2,3,4),
#         color.node = c('darkorchid', 'brown1', 'lightgreen', "dodgerblue2"), 
#         cutoff = 0.6, save = 'png', name.save = 'DIABLOnetwork')
# 
# 
# Y_unsup_DIABLO = as.data.frame(matrix(0, ncol = 3, nrow = ncol(RNA_gene)))
# colnames(Y_unsup_DIABLO) = c("case","intermediate","control")
# 
# for (i in seq(1:nrow(BDR_metadata))){
#   if (BDR_metadata$Case.Control[i] == "case"){
#     Y_unsup_DIABLO[i,1] = 1
#   } else if (BDR_metadata$Case.Control[i] == "intermediate"){
#     Y_unsup_DIABLO[i,2] = 1
#   } else if (BDR_metadata$Case.Control[i] == "control"){
#     Y_unsup_DIABLO[i,3] = 1
#   }
# }
# rownames(Y_unsup_DIABLO) = BDR_metadata$No.
# rownames(Y_unsup_DIABLO) = gsub("10230_BBN_", "BDR_", rownames(Y_unsup_DIABLO))
# 
# Y_unsup_DIABLO = as.matrix(Y_unsup_DIABLO)
# 
# multi_unsup_DIABLO = block.spls(X = X_diablo,
#                                 Y = Y_unsup_DIABLO,
#                                 ncomp = 2,
#                                 keepX = keepX_diablo,
#                                 scale = TRUE,
#                                 tol = 1e-06,
#                                 max.iter = 100,
#                                 near.zero.var = FALSE,
#                                 all.outputs = TRUE)
# 
# plotIndiv(multi_unsup_DIABLO) ## sample plot
# plotVar(multi_unsup_DIABLO, var.names = F, legend = TRUE, pch = c(20,20,20,20)) ## variable plot
# 
# plotIndiv(multi_unsup_DIABLO, 
#           ind.names = FALSE, 
#           legend=TRUE, cex=c(1), pch = c(15,18,20,4,5),
#           title = 'BRCA with DIABLO')
# 
# Y = multi_unsup_DIABLO$Y
# plotDiablo(multi_unsup_DIABLO, ncomp = 1)
# 
# plotLoadings(multi_unsup_DIABLO, comp = 1, contrib = "max")
# 
# network(multi_unsup_DIABLO, blocks = c(1,2,3,4),
#         color.node = c('darkorchid', 'brown1', 'lightgreen', "dodgerblue2"), 
#         cutoff = 0.94, save = 'png', name.save = 'DIABLOnetwork')
# 
# 
# 
# multi_iCB = iClusterBayes(dt1 = t(RNA_gene_norm),
#               dt2 = t(miRNA_norm),
#               dt3 = t(circRNA_norm),
#               dt4 = t(proteomics_norm),
#               type = c("gaussian","gaussian","gaussian","gaussian"),
#               K = 2)
# 
# png("testicb.png")
# plotHeatmap(multi_iCB,
#             datasets = X_diablo,
#             type = c("gaussian","gaussian","gaussian","gaussian"))
# dev.off()
# 
# iCB_sim_matrix(multi_iCB, label = colnames(RNA_gene_norm))
# 
# 
# 
# result=ExecuteSNF.CC(datasets=datasets, clusterNum=4, K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500, 
#                    title = "GBM", plot = "png", finalLinkage ="average")
# 
# 
# 
# group1=result$group
# distanceMatrix1=result$distanceMatrix
# p_value=survAnalysis(mainTitle="GBM1",time,status,group1,
#                      distanceMatrix1,similarity=TRUE)

  