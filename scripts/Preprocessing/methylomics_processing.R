library(wateRmelon)
library(sva)
library(matrixStats)
library(CETYGO)
setwd("D:/valentin/main/")
load("betas/ROSMAP_betas_only.Rdata")
pheno_ADC = read.table("pheno_ROSMAP_full.txt", sep = "\t", header = T)

#Doublecheck outliers 
outliers <- outlyx(mSet, plot=FALSE)

print(outliers) # Keep this

#Filter for outliers
msetEPIC_outliers = mSet[,mSet@phenoData@data$barcode %in% rownames(outliers[outliers$outliers == FALSE,])]
msetEPIC_outliers = msetEPIC_outliers[,msetEPIC_outliers@phenoData@data$barcode %in% pheno$Basename]

#Bisulfite conversion rates to filter low quality samples
bsc <- bscon(msetEPIC_outliers)
#Remove anything below 0.8

#Bead count and detection p-value filter for QC
msetEPIC_pfilter <- pfilter(msetEPIC_outliers)

#Methylation genders estimation
estimatesex = wateRmelon::estimateSex(betas(msetEPIC_pfilter), do_plot=FALSE)
#If mismatch remove samples

#CETYGO brain cell-type deconvolution
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

#BMIQ normalisation
msetEPIC_BMIQ <- BMIQ(msetEPIC_pfilter)

#Extract beta values
msetEPIC_betas <- betas(msetEPIC_BMIQ)

#Probe-wise outlier detection
pwod_bet <- pwod(msetEPIC_betas)
#Remove those if any

dim(msetEPIC_betas)

#Combat for batch effect correction
# msetEPIC_betas = ComBat(dat=msetEPIC_betas, batch=as.character(pheno$Institute), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
msetEPIC_betas = ComBat(dat=msetEPIC_betas, batch=as.character(pheno$Plate), mod=NULL, par.prior=TRUE, prior.plots=FALSE)


#Regression of covariates (age/sex/cell type)
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

#Removal of unwanted probes
SNPall = read.table("data/SNPProbes_McCartney.txt", header = T)
SNPallEURAF = SNPall[which(SNPall$EUR_AF >= 0.05 & SNPall$EUR_AF <= 0.95),]

crosshyb = read.table("data/CrossHydridisingProbes_McCartney.txt")

manifest_EPIC = read.csv("data/MethylationEPIC_v-1-0_B4.csv", header = T, skip = 7)
manifest_EPIC_XY = manifest_EPIC[manifest_EPIC$CHR %in% c("X", "Y"),]

to_remove_cpgs = unique(c(SNPallEURAF$IlmnID, manifest_EPIC_XY$IlmnID, crosshyb$V1))
#72240 probes

#remove non-cpg probes
mset_betas_lm = mset_betas_lm[-(which(rownames(mset_betas_lm) %in% to_remove_cpgs)),]
mset_betas_lm = mset_betas_lm[-grep("rs", rownames(mset_betas_lm)),]

#90 percentile absolute deviation for selection of most variant probes
mset_betas_lmdf = as.data.frame(mset_betas_lm)

df_values = data.frame(matrix(ncol = 2, nrow = nrow(mset_betas_lmdf)))
for(i in 1:nrow(mset_betas_lmdf)){
  tmp_cpg = as.numeric(mset_betas_lmdf[i,])
  value = c(mad(tmp_cpg), mad(tmp_cpg, center = quantile(tmp_cpg, probs = 0.9)))
  df_values[i,] = value
}

rownames(df_values) = rownames(mset_betas_lmdf)
colnames(df_values) = c("mad","90pc")

df_values = df_values[order(df_values$"90pc", decreasing = T),]
perc = 0.50
nfeatures = round(perc*length(df_values$"90pc"))

df_values_perc = df_values[seq(1:nfeatures),]

betas_ROSMAP = mset_betas_lmdf[rownames(mset_betas_lmdf) %in% rownames(df_values_perc),]

gc()
save.image(file = "ROSMAP_preproc_final.rdata")
