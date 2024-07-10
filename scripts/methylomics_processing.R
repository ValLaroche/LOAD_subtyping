library(wateRmelon)
library(sva)
library(matrixStats)
library(CETYGO)
setwd("D:/valentin/main/")
load("betas/ROSMAP_betas_only.Rdata")
pheno_ADC = read.table("pheno_ROSMAP_full.txt", sep = "\t", header = T)

#Doublecheck outliers from BDR
outliers <- outlyx(mSet, plot=FALSE)

print(outliers) # Keep this
# 1 outlier in the 454 set, non existant in the 383 set
# No outliers in the 240 UKBBN set
# 5 outliers in the 226 PITT set
# No outliers in 201 MSBB set

msetEPIC_outliers = mSet[,mSet@phenoData@data$barcode %in% rownames(outliers[outliers$outliers == FALSE,])]
msetEPIC_outliers = msetEPIC_outliers[,msetEPIC_outliers@phenoData@data$barcode %in% pheno$Basename]


bsc <- bscon(msetEPIC_outliers)
#Remove anything below 0.8
#everything around .95
#MSBB 202242410182_R05C01 = 65%

msetEPIC_pfilter <- pfilter(msetEPIC_outliers)
#Keep report of pfilter
#1 samples having 1 % of sites with a detection p-value greater than 0.05 were removed
# Samples removed: 202410160041_R04C01
# 844 sites were removed as beadcount <3 in 5 % of samples
# 13988 sites having 1 % of samples with a detection p-value greater than 0.05 were removed

#PITT
# 0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed
#Samples removed:
#1525 sites were removed as beadcount <3 in 5 % of samples
#4124 sites having 1 % of samples with a detection p-value greater than 0.05 were removed

#MSBB
# 0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed
# Samples removed:
#   965 sites were removed as beadcount <3 in 5 % of samples
# 4033 sites having 1 % of samples with a detection p-value greater than 0.05 were removed

estimatesex = wateRmelon::estimateSex(betas(msetEPIC_pfilter), do_plot=FALSE)

#If mismatch remove samples
# > table(pred_sex$predicted_sex)
# 
# 47,XXY Female   Male
# 2    186    193
# > table(pred_sex$truesex)
# 
# female   male
# 188    195
#2 males are outlier + pfilter removed
#2 females as 47 xxy

#PITT 1 sample 203968030048_R07C01

#MSBB 202242410223_R02C01 / 202262730156_R08C01 males as Females


# msetEPIC_ecc = estimateCellCounts.wmln(msetEPIC_pfilter, referencePlatform = "IlluminaHumanMethylationEPIC",
#                         compositeCellType = "DLPFC", cellTypes = c("NeuN_neg","NeuN_pos"))
# New estimator for 4 different celltypes

predPropAll<-list()

for(method in names(modelBrainCoef)){
  for(j in 1:length(modelBrainCoef[[method]])){
    if(!is.null(modelBrainCoef[[method]][[j]])){
      predPropAll[[method]][[j]]<-projectCellTypeWithError(betas(msetEPIC_pfilter), modelBrainCoef[[method]][[j]])
    }
  } 
}

predPropBest<-cbind(predPropAll[["ANOVA"]][[1]][,c("NeuNNeg_SOX10Neg", "NeuNPos")], predPropAll[["IDOL"]][[5]][,c("NeuNNeg_Sox10Neg_IRF8Pos","NeuNNeg_Sox10Neg_IRF8Neg")], predPropAll[["ANOVA"]][[3]][,c("SATB2Neg", "SATB2Pos")], predPropAll[["ANOVA"]][[6]][,c("NeuNNeg", "NeuNPos_SOX6Pos", "NeuNPos_SOX6Neg")], predPropAll[["IDOL"]][[4]][,"NeuNNeg_SOX10Pos"])
colnames(predPropBest)[10]<-"NeuNNeg_SOX10Pos"
# Panel 8 looks like the most relevant (Neun/Exc | Neun/Inh | Microglia | Astrocyte | Oligodendrocytes)


msetEPIC_BMIQ <- BMIQ(msetEPIC_pfilter)

msetEPIC_betas <- betas(msetEPIC_BMIQ)

pwod_bet <- pwod(msetEPIC_betas)
#Remove those if any


# pheno_PFC = pheno_PFC[pheno_PFC$prop >= 0.3,]


msetEPIC_betas = msetEPIC_betas[-grep("ch.", rownames(msetEPIC_betas)),]
dim(msetEPIC_betas)

# msetEPIC_betas = ComBat(dat=msetEPIC_betas, batch=as.character(pheno$Institute), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
msetEPIC_betas = ComBat(dat=msetEPIC_betas, batch=as.character(pheno$Plate), mod=NULL, par.prior=TRUE, prior.plots=FALSE)




# New resid loop with imputation of missing values
mset_betas_lm = betas_cohort
mset_betas_lm = mset_betas_lm[complete.cases(mset_betas_lm),]

for(i in seq(1:nrow(mset_betas_lm))){
  mset_betas_lm[i,] = lm(mset_betas_lm[i,] ~ pheno_ADC$age_death + 
                           pheno_ADC$msex + 
                           pheno_ADC$Sample_Plate +
                           pheno_ADC$Neuronal_Inhibitory + 
                           pheno_ADC$Neuronal_Excitatory + 
                           pheno_ADC$Microglia + 
                           pheno_ADC$Astrocyte + 
                           pheno_ADC$Oligodendrocyte, na.action = na.exclude)$residuals
}

SNPall = read.table("data/SNPProbes_McCartney.txt", header = T)
SNPallEURAF = SNPall[which(SNPall$EUR_AF >= 0.05 & SNPall$EUR_AF <= 0.95),]

crosshyb = read.table("data/CrossHydridisingProbes_McCartney.txt")

manifest_EPIC = read.csv("data/MethylationEPIC_v-1-0_B4.csv", header = T, skip = 7)
manifest_EPIC_XY = manifest_EPIC[manifest_EPIC$CHR %in% c("X", "Y"),]

to_remove_cpgs = unique(c(SNPallEURAF$IlmnID, manifest_EPIC_XY$IlmnID, crosshyb$V1))
#72240 BDR
mset_betas_lm = mset_betas_lm[-(which(rownames(mset_betas_lm) %in% to_remove_cpgs)),]
mset_betas_lm = mset_betas_lm[-grep("rs", rownames(mset_betas_lm)),]

mset_betas_lmdf = as.data.frame(mset_betas_lm)

df_values = data.frame(matrix(ncol = 2, nrow = nrow(mset_betas_lmdf)))
### Use 80-90% quantile absolute deviation ?
for(i in 1:nrow(mset_betas_lmdf)){
  tmp_cpg = as.numeric(mset_betas_lmdf[i,])
  value = c(mad(tmp_cpg), mad(tmp_cpg, center = quantile(tmp_cpg, probs = 0.9)))
  df_values[i,] = value
}

rownames(df_values) = rownames(mset_betas_lmdf)
colnames(df_values) = c("mad","90pc")

# hist(df_values$mad)
# hist(df_values$"90pc")

df_values = df_values[order(df_values$"90pc", decreasing = T),]
perc = 0.50
nfeatures = round(perc*length(df_values$"90pc"))

df_values_perc = df_values[seq(1:nfeatures),]

betas_ROSMAP = mset_betas_lmdf[rownames(mset_betas_lmdf) %in% rownames(df_values_perc),]

# hist(df_values$"90pc")

gc()
save.image(file = "BDR_preproc_final.rdata")


#####
tmset_betas_lm = as.data.frame(t(mset_betas_lm))
identical(rownames(tmset_betas_lm), pheno$Basename)
tmset_betas_lm$Basename = pheno$Basename
tmset_betas_lm$Cell_type

list_all_lme = list(unique(tmset_betas_lm$Cell_type))

for(celltype in unique(pheno$Cell_Type)){
  pheno$Cell_type_allv1 = pheno$Cell_Type
  pheno[pheno$Cell_type_allv1 == celltype,]$Cell_type_allv1 = celltype
  pheno[pheno$Cell_type_allv1 != celltype,]$Cell_type_allv1 = "Other"
  
  tmset_betas_lm$Cell_type = pheno$Cell_type_allv1
  
  
  for(cpg in colnames(tmset_betas_lm)[1:(ncol(tmset_betas_lm)-2)]){
    tmp_df = data.frame(tmset_betas_lm[[cpg]], tmset_betas_lm$Cell_type, tmset_betas_lm$Basename)
    colnames(tmp_df) = c("cpg", "Cell_type", "Basename")
    cpg_val = tmset_betas_lm[[cpg]]
    Cell_type = tmset_betas_lm$Cell_type
    Basename = tmset_betas_lm$Basename
    cpg_lme = nlme::lme(cpg_val~Cell_type, random=~1|Basename,  control = lmeControl(opt = "optim"))
    list_all_lme[[celltype]][[cpg]] = summary(cpg_lme)$tTable
  }
    
}



