WGCNAblock <- function(){

##### Loading #####
  
setwd("D:/valentin/main/")
source("./scripts/utils.R")
source("./scripts/WGCNA/WGCNA_functions.R")
#Load WGCNA model (per cohort)
load("./WGCNA/PITT/PITT_subtyping.wgcna.network.rdat")
load("WGCNA/PITT/PITT_subtyping.wgcna.softThreshold.rdat")
load("WGCNA/PITT/PITT_subtyping.wgcna1.rdat")

load("betas/ROSMAP_final_set_27-11-23.Rdata")
load("betas/UKBBN_final_set_27-9-2023.Rdata")
load("betas/PITT_final_set_15-11-23.Rdata")
load("phenos_3cohorts.Rdata")
load("D:/valentin/main/pheno_ROSMAP.Rdata")
##### Prepare data #####
  betas_UKBBN = betas_UKBBN[rownames(betas_UKBBN) %in% intersect(rownames(betas_UKBBN), rownames(betas_PITT)),]
  betas_PITT = betas_PITT[rownames(betas_PITT) %in% intersect(rownames(betas_UKBBN), rownames(betas_PITT)),]
  betas_ROSMAP = betas_ROSMAP[rownames(betas_ROSMAP) %in% rownames(betas_PITT),]
  
  betas_UKBBN = betas_UKBBN[,order(colnames(betas_UKBBN))]
  betas_PITT = betas_PITT[,order(colnames(betas_PITT))]
  betas_ROSMAP = betas_ROSMAP[,order(colnames(betas_ROSMAP))]
  
  pheno_UKBBN = pheno_UKBBN[order(pheno_UKBBN$Basename),]
  pheno_PITT = pheno_PITT[order(pheno_PITT$Basename),]
  pheno_ROSMAP = pheno_ROSMAP[order(pheno_ROSMAP$Basename),]
  
  identical(colnames(betas_UKBBN), pheno_UKBBN$Basename)
  identical(colnames(betas_PITT), pheno_PITT$Basename)
  identical(colnames(betas_ROSMAP), pheno_ROSMAP$Basename)
  
  data_matrix_AD_UKBBN = betas_UKBBN[,colnames(betas_UKBBN) %in% pheno_UKBBN[pheno_UKBBN$HC_clusters != "Control",]$Basename]
  data_matrix_AD_PITT = betas_PITT[,colnames(betas_PITT) %in% pheno_PITT[pheno_PITT$HC_clusters != "Control",]$Basename]
  data_matrix_AD_ROSMAP = betas_ROSMAP[,colnames(betas_ROSMAP) %in% pheno_ROSMAP[pheno_ROSMAP$HC_clusters != "Control",]$Basename]
  
  # Generate eigen values for each module and each cohort
  # Which cohort to use in the WGCNA model
  cohort = "ROSMAP"
  WGCNA_labels = net$colors
  inter_450k_EPIC = rownames(data_matrix_AD_ROSMAP)
  
  data_matrix_AD_UKBBN = data_matrix_AD_UKBBN[rownames(data_matrix_AD_UKBBN) %in% names(WGCNA_labels),]
  data_matrix_AD_PITT = data_matrix_AD_PITT[rownames(data_matrix_AD_PITT) %in% names(WGCNA_labels),]

  #Generate the 450k sub-WGCNA model
  WGCNA_labels_450k = WGCNA_labels[names(WGCNA_labels) %in% inter_450k_EPIC]
  table(WGCNA_labels)
  table(WGCNA_labels_450k)
  
  #Calculate eigengenes values EPIC
  output = gen_eigengenes(t(data_matrix_AD_UKBBN),
                          t(data_matrix_AD_PITT),
                          NULL,
                          WGCNA_labels)

  UKBBN_eigen_EPIC = output[[1]]
  PITT_eigen_EPIC = output[[2]]

  #Calculate eigengenes values 450k
  matrix_450k_UKBBN = data_matrix_AD_UKBBN[rownames(data_matrix_AD_UKBBN) %in% inter_450k_EPIC,]
  matrix_450k_PITT = data_matrix_AD_PITT[rownames(data_matrix_AD_PITT) %in% inter_450k_EPIC,]

  output = gen_eigengenes(t(matrix_450k_UKBBN),
                          t(matrix_450k_PITT),
                          t(data_matrix_AD_ROSMAP),
                          WGCNA_labels_450k)

  UKBBN_eigen_450k = output[[1]]
  PITT_eigen_450k = output[[2]]
  ROSMAP_eigen_450k = output[[3]]

  #Generate AD-only pheno files
  cohort = "PITT"
  clustering = "HC_clusters"
  dir.create(paste0("./WGCNA/", cohort, "/", clustering))
  output = gen_phenos(pheno_1 = pheno_UKBBN,
                      pheno_2 = pheno_PITT,
                      pheno_3 = pheno_ROSMAP,
                      clustering = clustering)
  
  WGCNA_phenoUKBBN = output[[1]]
  WGCNA_phenoPITT = output[[2]]
  WGCNA_phenoROSMAP = output[[3]]
  
  #Generate standard names for datasets and future function
  if(cohort == "PITT"){
    main_pheno = WGCNA_phenoPITT
    pheno_disc_1 = WGCNA_phenoUKBBN
    pheno_disc_3 = WGCNA_phenoROSMAP
    
    main_expr_EPIC = data_matrix_AD_PITT
    main_expr_450k = matrix_450k_PITT
    
    main_eEPIC = PITT_eigen_EPIC
    main_e450k = PITT_eigen_450k
    
    disc_expr_EPIC_1 = data_matrix_AD_UKBBN
    disc_expr_450k_1 = matrix_450k_UKBBN
    
    disc_eEPIC_1 = UKBBN_eigen_EPIC
    disc_e450k_1 = UKBBN_eigen_450k
    
    disc_expr_450k_3 = data_matrix_AD_ROSMAP
    
    pref_1 = "PITT"
    pref_2 = "UKBBN"
    pref_4 = "ROSMAP"
  } else if(cohort == "UKBBN"){

    main_pheno = WGCNA_phenoUKBBN
    pheno_disc_1 = WGCNA_phenoPITT
    pheno_disc_3 = WGCNA_phenoROSMAP
    
    main_expr_EPIC = data_matrix_AD_UKBBN
    main_expr_450k = matrix_450k_UKBBN
    
    main_eEPIC = UKBBN_eigen_EPIC
    main_e450k = UKBBN_eigen_450k
    
    disc_expr_EPIC_1 = data_matrix_AD_PITT
    disc_expr_450k_1 = matrix_450k_PITT
    
    disc_eEPIC_1 = PITT_eigen_EPIC
    disc_e450k_1 = PITT_eigen_450k
    
    disc_expr_450k_3 = data_matrix_AD_ROSMAP
    
    pref_1 = "UKBBN"
    pref_2 = "PITT"
    pref_4 = "ROSMAP"
  } else if(cohort == "ROSMAP"){
    main_pheno = WGCNA_phenoROSMAP
    pheno_disc_1 = WGCNA_phenoUKBBN
    pheno_disc_3 = WGCNA_phenoPITT
    
    main_expr_EPIC = data_matrix_AD_ROSMAP
    main_expr_450k = data_matrix_AD_ROSMAP
    
    disc_expr_EPIC_1 = matrix_450k_UKBBN
    disc_expr_450k_1 = matrix_450k_UKBBN
    
    disc_expr_EPIC_3 = matrix_450k_PITT
    disc_expr_450k_3 = matrix_450k_PITT
    
    pref_1 = "ROSMAP"
    pref_2 = "UKBBN"
    pref_4 = "PITT"
  } 
 
   
##### WGCNA outputs #####
  # Correlation of WGCNA with clusters
  trait_UKBBN = gen_subtype_matrix(WGCNA_phenoUKBBN, clustering = clustering, nsubtypes = length(unique(WGCNA_phenoUKBBN[[clustering]])))
  heatmap_module_to_subtype(eigengenes = UKBBN_eigen_EPIC, pheno_trait = trait_UKBBN, 
                            paste0("./WGCNA/", cohort, "/", clustering, "/corUKBBN_heatmap.pdf"))
  
  trait_PITT = gen_subtype_matrix(WGCNA_phenoPITT, clustering = clustering, nsubtypes = 3)
  heatmap_module_to_subtype(eigengenes = PITT_eigen_EPIC, pheno_trait = trait_PITT, 
                            paste0("./WGCNA/", cohort, "/", clustering, "/corPITT_heatmap.pdf"))
  
  trait_ROSMAP = gen_subtype_matrix(WGCNA_phenoROSMAP, clustering = clustering, nsubtypes = 3)
  heatmap_module_to_subtype(eigengenes = ROSMAP_eigen_450k, pheno_trait = trait_ROSMAP, 
                            paste0("./WGCNA/", cohort, "/", clustering, "/corROSMAP_heatmap.pdf"))
  
  # Calculating anova/tukey tests for clusters and modules
  gen_anova_tukey_tests(pheno = WGCNA_phenoUKBBN, clustering = clustering, 
                        eigengenes = UKBBN_eigen_EPIC, filename = paste0("./WGCNA/", cohort, "/", clustering, "/anovatukeyU.pdf"))
  gen_anova_tukey_tests(pheno = WGCNA_phenoPITT, clustering = clustering, 
                        eigengenes = PITT_eigen_EPIC, filename = paste0("./WGCNA/", cohort, "/", clustering, "/anovatukeyP.pdf"))
  gen_anova_tukey_tests(pheno = WGCNA_phenoROSMAP, clustering = clustering, 
                        eigengenes = ROSMAP_eigen_450k, filename = paste0("./WGCNA/", cohort, "/", clustering, "/anovatukeyR.pdf"))
  
  # Comparing either EPIC/450K or two cohorts to one another (ANOVA/Tukey)
  double_anovas(pheno_1 = pheno_disc_1, eigen_1 = disc_eEPIC_1, pref_1 = "UKBBN_EPIC",
                pheno_2 = pheno_disc_1, eigen_2 = disc_e450k_1, pref_2 = "UKBBN_450k",
                clustering = clustering, filename = paste0("./WGCNA/", cohort, "/", clustering, "/anovas2cohort_UE45_UKBBN.pdf"))
  double_anovas(pheno_1 = main_pheno, eigen_1 = main_eEPIC, pref_1 = "PITT_EPIC",
                pheno_2 = main_pheno, eigen_2 = main_e450k, pref_2 = "PITT_450k",
                clustering = clustering, filename = paste0("./WGCNA/", cohort, "/", clustering, "/anovas2cohort_UE45_PITT.pdf"))
  double_anovas(pheno_1 = main_pheno, eigen_1 = main_eEPIC, pref_1 = pref_1,
                pheno_2 = pheno_disc_2, eigen_2 = disc_eEPIC_2, pref_2 = pref_3,
                clustering = clustering, paste0("./WGCNA/", cohort, "/", clustering, "/anovas2cohort_UP.pdf"))
  
  # Comparing all cohorts to one another (ANOVA/Tukey)
  triple_anovas(pheno_1 = pheno_disc_1, eigen_1 = disc_eEPIC_1, pref_1 = pref_2,
                pheno_2 = main_pheno, eigen_2 = main_eEPIC, pref_2 = pref_1,
                pheno_3 = pheno_disc_2, eigen_3 = disc_eEPIC_2, pref_3 = pref_3,
                pheno_4 = pheno_disc_3, eigen_4 = disc_e450k_3, pref_4 = pref_4,
                clustering = clustering, paste0("./WGCNA/", cohort, "/", clustering, "/anovas3cohort_noBDR.pdf"))
  
  # Get loadings of modules 
  gen_comploading_plots(pheno_1 = pheno_disc_1, eigen_1 = disc_eEPIC_1, cohort_1 = pref_1,
                        pheno_2 = main_pheno, eigen_2 = main_eEPIC, cohort_2 = pref_2,
                        pheno_3 = pheno_disc_2, eigen_3 = disc_eEPIC_2, cohort_3 = pref_3,
                        pheno_4 = pheno_disc_3, eigen_4 = disc_e450k_3, cohort_4 = pref_4,
                        clustering = clustering, title = "Loadings_UKBBNnet",
                        filename = paste0("./WGCNA/", cohort, "/", clustering, "/loadings.pdf"))
  
  # Bootstrap randomization of clusters per module (specificity check)
  gen_random_anovas(pheno_1 = main_pheno, eigen_1 = main_eEPIC, pref_1 = pref_1,
                    pheno_2 = pheno_disc_1, eigen_2 = disc_eEPIC_1, pref_2 = pref_2,
                    pheno_3 = pheno_disc_2, eigen_3 = disc_eEPIC_2, pref_3 = pref_3,
                    pheno_4 = pheno_disc_3, eigen_4 = disc_e450k_3, pref_4 = pref_4,
                    clustering = clustering,
                    filename = paste0("./WGCNA/", cohort, "/", clustering, "/randomized_noBDR.pdf"))
  
##### Subtypes #####
  # Mapping one replication cohort on the discovery
  pred_mapping_1pred(data_disc = t(main_expr_EPIC), pheno_disc = main_pheno,
                     data_pred = t(disc_expr_EPIC_1), pheno_pred = pheno_disc_1, cohort = cohort, disc_cohort = pref_2,
                     clustering = clustering, filename = paste0("./WGCNA/", cohort, "/", clustering, "/mapping_to_", pref_2, "EPIC.pdf"))
  pred_mapping_1pred(data_disc = t(main_expr_450k), pheno_disc = main_pheno,
                     data_pred = t(disc_expr_450k_3), pheno_pred = pheno_disc_3, cohort = cohort, disc_cohort = pref_4,
                     clustering = clustering, filename = paste0("./WGCNA/", cohort, "/", clustering, "/mapping_to_", pref_4, "450k.pdf"))
  
  # Testing over-representation of the replication in the discovery for either EPIC or 450K sets
  df_all_tests_fisher = data.frame(matrix(nrow = 0, ncol = 5))
  for(subtype_hull in unique(main_pheno[[clustering]])){
    tmp_fisher = test_chulls(data_disc = t(main_expr_EPIC), pheno_disc = main_pheno,
                       data_pred = t(disc_expr_EPIC_1), pheno_pred = pheno_disc_1, 
                       clustering = clustering, disc = pref_2)
    df_all_tests_fisher = rbind(df_all_tests_fisher, tmp_fisher)

  }
  for(subtype_hull in unique(main_pheno[[clustering]])){
  tmp_fisher = test_chulls(data_disc = t(main_expr_450k), pheno_disc = main_pheno,
                           data_pred = t(disc_expr_450k_3), pheno_pred = pheno_disc_3, 
                           clustering = clustering, disc = pref_4)
  df_all_tests_fisher = rbind(df_all_tests_fisher, tmp_fisher)
  }
  
  df_all_tests_fisher$P_value_Fisher = signif(df_all_tests_fisher$P_value_Fisher, digits = 3)
  
  # Output table
  pdf(paste0("./WGCNA/", cohort, "/", clustering, "/table_all_fisher.pdf"), width = 10, height = 9)
  grid.arrange(tableGrob(df_all_tests_fisher))
  dev.off()
  
  # Mapping all replication cohorts on the discovery
  pred_mapping_multpred(data_disc = t(main_expr_EPIC), pheno_disc = main_pheno,
                     data_pred_1 = t(disc_expr_EPIC_1), pheno_pred_1 = pheno_disc_1, 
                     data_pred_2 = t(disc_expr_EPIC_2), pheno_pred_2 = pheno_disc_2, 
                     clustering = clustering, filename = paste0("./WGCNA/", cohort, "/", clustering, "/mapping_to_cohortsEPIC.pdf"))
  pred_mapping_multpred(data_disc = t(main_expr_450k), pheno_disc = main_pheno,
                     data_pred_1 = t(disc_expr_450k_1), pheno_pred_1 = pheno_disc_1, 
                     data_pred_2 = t(disc_expr_450k_2), pheno_pred_2 = pheno_disc_2, 
                     data_pred_3 = t(disc_expr_450k_3), pheno_pred_3 = pheno_disc_3, 
                     clustering = clustering, filename = paste0("./WGCNA/", cohort, "/", clustering, "/mapping_to_cohorts450k.pdf"))
  
  
  pdf(paste0("./WGCNA/", cohort, "/", clustering, "/median_all_cohorts.pdf"), width = 9, height = 9)

##### Correlations Red/Blue #####
  
  # Adding subtype information in the pheno files
  red_methods_U = WGCNA_phenoUKBBN[WGCNA_phenoUKBBN$HC_clusters == 2 & WGCNA_phenoUKBBN$KM_clusters == 2,]
  red_methods_P = WGCNA_phenoPITT[WGCNA_phenoPITT$HC_clusters %in% c(3) & WGCNA_phenoPITT$KM_clusters %in% c(3),]
  red_methods_R = WGCNA_phenoROSMAP[WGCNA_phenoROSMAP$HC_clusters %in% 2 & WGCNA_phenoROSMAP$KM_clusters %in% 2,]
  red_samples = c(red_methods_P$Basename,
                  red_methods_R$Basename,
                  red_methods_U$Basename)
  
  blue_methods_U = WGCNA_phenoUKBBN[WGCNA_phenoUKBBN$HC_clusters == "1" & WGCNA_phenoUKBBN$KM_clusters == "1",]
  blue_methods_P = WGCNA_phenoPITT[WGCNA_phenoPITT$HC_clusters %in% 2 & WGCNA_phenoPITT$KM_clusters %in% 2,]
  blue_methods_R = WGCNA_phenoROSMAP[WGCNA_phenoROSMAP$HC_clusters %in% c(1) & WGCNA_phenoROSMAP$KM_clusters %in% c(1),]
  blue_samples = c(blue_methods_P$Basename,
                  blue_methods_R$Basename,
                  blue_methods_U$Basename)
  
  red_methods_U$subtype = "red"
  red_methods_R$subtype = "red"
  red_methods_P$subtype = "red"
  red_methods_U$cohort = "UKBBN"
  red_methods_R$cohort = "ROSMAP"
  red_methods_P$cohort = "PITT"
  
  blue_methods_U$subtype = "blue"
  blue_methods_R$subtype = "blue"
  blue_methods_P$subtype = "blue"
  blue_methods_U$cohort = "UKBBN"
  blue_methods_R$cohort = "ROSMAP"
  blue_methods_P$cohort = "PITT"
  
  table_redblue = data.frame("Basename" = c(red_methods_B$Basename,red_methods_R$Basename,red_methods_U$Basename,red_methods_P$Basename,
                                            blue_methods_B$Basename,blue_methods_R$Basename,blue_methods_U$Basename,blue_methods_P$Basename),
                             "subtype" = c(red_methods_B$subtype,red_methods_R$subtype,red_methods_U$subtype,red_methods_P$subtype,
                                            blue_methods_B$subtype,blue_methods_R$subtype,blue_methods_U$subtype,blue_methods_P$subtype),
                             "cohort" = c(red_methods_B$cohort,red_methods_R$cohort,red_methods_U$cohort,red_methods_P$cohort,
                                           blue_methods_B$cohort,blue_methods_R$cohort,blue_methods_U$cohort,blue_methods_P$cohort))
  
  
  # Making a global matrix
  data_matrix_all = matrix(ncol = 0, nrow = 148951)  
  data_matrix_all = cbind(data_matrix_all, matrix_450k_PITT[,colnames(matrix_450k_PITT) %in% table_redblue$Basename])
  data_matrix_all = cbind(data_matrix_all, matrix_450k_UKBBN[,colnames(matrix_450k_UKBBN) %in% table_redblue$Basename])
  data_matrix_all = cbind(data_matrix_all, data_matrix_AD_ROSMAP[,colnames(data_matrix_AD_ROSMAP) %in% table_redblue$Basename])

  WGCNA_phenoUKBBN$subtype = NA
  WGCNA_phenoUKBBN[WGCNA_phenoUKBBN$Basename %in% red_methods_U$Basename,]$subtype = "red"
  WGCNA_phenoUKBBN[WGCNA_phenoUKBBN$Basename %in% blue_methods_U$Basename,]$subtype = "blue"
  WGCNA_phenoUKBBN[is.na(WGCNA_phenoUKBBN$subtype),]$subtype = "grey90"
  
  WGCNA_phenoPITT$subtype = NA
  WGCNA_phenoPITT[WGCNA_phenoPITT$Basename %in% red_methods_P$Basename,]$subtype = "red"
  WGCNA_phenoPITT[WGCNA_phenoPITT$Basename %in% blue_methods_P$Basename,]$subtype = "blue"
  WGCNA_phenoPITT[is.na(WGCNA_phenoPITT$subtype),]$subtype = "grey90"
  
  WGCNA_phenoROSMAP$subtype = NA
  WGCNA_phenoROSMAP[WGCNA_phenoROSMAP$Basename %in% red_methods_R$Basename,]$subtype = "red"
  WGCNA_phenoROSMAP[WGCNA_phenoROSMAP$Basename %in% blue_methods_R$Basename,]$subtype = "blue"
  WGCNA_phenoROSMAP[is.na(WGCNA_phenoROSMAP$subtype),]$subtype = "grey90"
  
  table(WGCNA_phenoUKBBN$subtype)
  table(WGCNA_phenoPITT$subtype)
  table(WGCNA_phenoROSMAP$subtype)
  
  # Output all the cohorts in an R object
  pheno_PITT$subtype = NA
  pheno_PITT[!is.na(match(pheno_PITT$Basename, WGCNA_phenoPITT$Basename)),]$subtype = WGCNA_phenoPITT$subtype
  pheno_PITT[is.na(pheno_PITT$subtype),]$subtype = "Control"
  pheno_UKBBN$subtype = NA
  pheno_UKBBN[!is.na(match(pheno_UKBBN$Basename, WGCNA_phenoUKBBN$Basename)),]$subtype = WGCNA_phenoUKBBN$subtype
  pheno_UKBBN[is.na(pheno_UKBBN$subtype),]$subtype = "Control"
  pheno_ROSMAP$subtype = NA
  pheno_ROSMAP[!is.na(match(pheno_ROSMAP$Basename, WGCNA_phenoROSMAP$Basename)),]$subtype = WGCNA_phenoROSMAP$subtype
  pheno_ROSMAP[is.na(pheno_ROSMAP$subtype),]$subtype = "Control"
  
  save(pheno_PITT, pheno_ROSMAP, pheno_UKBBN, file = "phenos_meth.Rdata")
  
  grey_phenoUKBBN = WGCNA_phenoUKBBN[WGCNA_phenoUKBBN$subtype == "grey90",]
  grey_phenoPITT = WGCNA_phenoPITT[WGCNA_phenoPITT$subtype == "grey90",]
  grey_phenoROSMAP = WGCNA_phenoROSMAP[WGCNA_phenoROSMAP$subtype == "grey90",]
  
  #####
  
  # Generating median corrplots (This is the pretty output, after looking at the correlations a first time)  
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  splsda_redblue = mixOmics::splsda(t(main_expr_450k), main_pheno$HC_clusters, ncomp = 6)
  
  cor_median = median_samples(pheno_1 = main_pheno, data_1 = main_expr_450k, pref_1 = pref_1,
                              pheno_2 = pheno_disc_1, data_2 = disc_expr_450k_1, pref_2 = pref_2,
                              pheno_3 = pheno_disc_3, data_3 = disc_expr_450k_3, pref_3 = pref_4,
                              clustering = "HC_clusters")
  
  names_tosub = c("blue", "blue", "blue", "red",
                  "red", "red", "grey", "grey", "grey")
  pdf(file = paste0("WGCNA/",cohort,"/testKM_corrplot_redblue.pdf"), width = 10, height = 10)
  cplot = corrplot(cor(t(cor_median)), method = "color", order = "hclust", hclust.method = "ward.D2",
                   tl.col = "black", addCoef.col = "black", insig = "blank", type = "upper")
  dev.off()
  pred_1 = predict(splsda_redblue, newdata = t(data_matrix_AD_ROSMAP))
  pred_2 = predict(splsda_redblue, newdata = t(data_matrix_AD_ROSMAP)) 
  pred_3 = predict(splsda_redblue, newdata = t(matrix_450k_PITT))
  pred_median = predict(splsda_redblue, newdata = cor_median)
  
  col_PITT = WGCNA_phenoPITT$subtype
  col_PITT[col_PITT == "blue"] = "steelblue1"
  col_PITT[col_PITT == "red"] = "brown2"
  
  col_UKBBN = WGCNA_phenoUKBBN$subtype
  col_UKBBN[col_UKBBN == "blue"] = "steelblue1"
  col_UKBBN[col_UKBBN == "red"] = "brown2"
  
  col_ROSMAP = WGCNA_phenoROSMAP$subtype
  col_ROSMAP[col_ROSMAP == "blue"] = "steelblue1"
  col_ROSMAP[col_ROSMAP == "red"] = "brown2"
  
  plot_redblue = plotIndiv(splsda_redblue, comp = 1:2, rep.space = "X-variate",
            style="graphics",ind.names=FALSE, legend = FALSE, col = c("white","white","white"),
            title = paste0("KMeans ", cohort, " discovery map"))
  
  pdf(file = paste0("WGCNA/",cohort,"/HC_indmap_redblue_noBDR.pdf"), width = 9, height = 9)
  plotIndiv(splsda_redblue, comp = 1:2, rep.space = "X-variate",
            style="graphics",ind.names=FALSE, legend = FALSE, col = c("white","white","white"),
            title = paste0("HClust ", cohort, " discovery map"),
            xlim = c(min(plot_redblue$df$x)-100,max(plot_redblue$df$x)+100))
  points(pred_1$variates[, 1], pred_1$variates[, 2], 
         pch = 18, cex = 1.4, col = col_ROSMAP)
  points(pred_2$variates[, 1], pred_2$variates[, 2], 
         pch = 18, cex = 1.4, col = col_ROSMAP)
  points(splsda_redblue$variates$X[,1], splsda_redblue$variates$X[,2],
         pch = 16, cex = 1, col = col_UKBBN)
  points(pred_3$variates[, 1], pred_3$variates[, 2], 
         pch = 17, cex = 1, col = col_PITT)
  legend("topleft",title = "Cohort",legend = c("UKBBN","PITT","ROSMAP"), pch = c(16,17,18), cex = 1.1)
  colors_med = c("#000080", "grey50", "#000080", "#420D09", "#000080", "#420D09",
                 "grey50", "#420D09", "grey50")
  text(pred_median$variates[, 1], pred_median$variates[, 2], 
       pch = 1, cex = 1, col = colors_med, labels = rownames(cor_median))
  dev.off()

#### Test PCA ####
identical(rownames(matrix_450k_PITT), rownames(matrix_450k_UKBBN))
identical(rownames(matrix_450k_PITT), rownames(main_expr_450k))
  
all_data_betas = cbind(disc_expr_450k_1, disc_expr_450k_3, main_expr_450k)
pheno_ROSMAP_topca = data.frame(WGCNA_phenoROSMAP$Basename, WGCNA_phenoROSMAP$diag, WGCNA_phenoROSMAP$subtype, WGCNA_phenoROSMAP$HC_clusters, WGCNA_phenoROSMAP$KM_clusters)
pheno_ROSMAP_topca$cohort = "ROSMAP"
colnames(pheno_ROSMAP_topca) = c("Basename", "diag", "subtype", "HC_clusters", "KM_clusters", "cohort")
pheno_UKBBN_topca = data.frame(WGCNA_phenoUKBBN$Basename, WGCNA_phenoUKBBN$diag, WGCNA_phenoUKBBN$subtype, WGCNA_phenoUKBBN$HC_clusters, WGCNA_phenoUKBBN$KM_clusters)
pheno_UKBBN_topca$cohort = "UKBBN"
colnames(pheno_UKBBN_topca) = c("Basename", "diag", "subtype", "HC_clusters", "KM_clusters", "cohort")
pheno_PITT_topca = data.frame(WGCNA_phenoPITT$Basename, WGCNA_phenoPITT$diag, WGCNA_phenoPITT$subtype, WGCNA_phenoPITT$HC_clusters, WGCNA_phenoPITT$KM_clusters)
pheno_PITT_topca$cohort = "PITT"
colnames(pheno_PITT_topca) = c("Basename", "diag", "subtype", "HC_clusters", "KM_clusters", "cohort")

pheno_all_topca = rbind(pheno_PITT_topca, pheno_ROSMAP_topca, pheno_UKBBN_topca)
pheno_all_topca$subtypenum = pheno_all_topca$subtype
pheno_all_topca[pheno_all_topca$subtypenum == "grey90",]$subtypenum = "3"
pheno_all_topca[pheno_all_topca$subtypenum == "blue",]$subtypenum = "1"
pheno_all_topca[pheno_all_topca$subtypenum == "red",]$subtypenum = "2"

identical(pheno_all_topca$Basename, colnames(all_data_betas))
set.seed(3)
pca_all_data = prcomp(t(all_data_betas), scale. = TRUE)

out_pca = data.frame(pca_all_data$x, subtype = pheno_all_topca$subtype)

ggplot(out_pca, aes(x = PC1, y = PC2, col = subtype)) +
  geom_point(size=2) +
  scale_color_manual(values = c("steelblue2","grey90","brown2"))+ 
  theme_classic()

}
