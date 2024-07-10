library(wateRmelon)
library(sva)
library(matrixStats)
library(CETYGO)
library(nlme)
library(foreach)
library(doParallel)
library(pbapply)
setwd("/lustre/projects/Research_Project-191391/valentin/Subtyping/FACS-Celltype/")
load("Celltype_study.rdata")
args<-commandArgs(TRUE)

celltype <- args[1]

cores=detectCores()
cl <- makeCluster(cores[1]-1)

mset_betas_tolm = msetEPIC_betas
mset_betas_tolm = mset_betas_tolm[complete.cases(mset_betas_tolm),]
pheno$Age = as.numeric(pheno$Age)

# Control for covariates
resid <- function(row, age, sex){
  fit <- try(
    lm( row ~ age + sex),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
}

age = pheno$Age
sex = pheno$Sex

clusterExport(cl, "sex")
clusterExport(cl, "age")
clusterExport(cl, "resid")

mset_betas_lm = t(pbapply(mset_betas_tolm, 1, resid, age, sex, cl = cl))

# Remove unwanted probes
SNPall = read.table("SNPProbes_McCartney.txt", header = T)
SNPallEURAF = SNPall[which(SNPall$EUR_AF >= 0.05 & SNPall$EUR_AF <= 0.95),]

crosshyb = read.table("CrossHydridisingProbes_McCartney.txt")

manifest_EPIC = read.csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv", header = T, skip = 7)
manifest_EPIC_XY = manifest_EPIC[manifest_EPIC$CHR %in% c("X", "Y"),]

to_remove_cpgs = unique(c(SNPallEURAF$IlmnID, manifest_EPIC_XY$IlmnID, crosshyb$V1))

mset_betas_lm = mset_betas_lm[-(which(rownames(mset_betas_lm) %in% to_remove_cpgs)),]
mset_betas_lm = mset_betas_lm[-grep("rs", rownames(mset_betas_lm)),]

print(celltype)
if(celltype == "Trip"){
  celltype = "Trip neg"
}

# Put selected cell type against the rest
pheno$Cell_type_allv1 = pheno$Cell_Type
pheno[pheno$Cell_type_allv1 == celltype,]$Cell_type_allv1 = celltype
pheno[pheno$Cell_type_allv1 != celltype,]$Cell_type_allv1 = "Other"

#tmset_betas_lm$Cell_type = pheno$Cell_type_allv1


Cell_type = pheno$Cell_type_allv1
Basename = pheno$Basename

# Perform LME, accounting for cell-type and individuals
f_lme <- function(row, Cell_type, Basename){
  
  inputdata  = data.frame(row = row, Cell_type = Cell_type, Basename = Basename)
  cpg_lme = nlme::lme(row~Cell_type, random=~1|Basename,  control = nlme::lmeControl(opt = "optim"))
  sum_cpg = c(summary(cpg_lme)$tTable[2,])
  return(sum_cpg)
  
}

clusterExport(cl, "Cell_type")
clusterExport(cl, "Basename")
clusterExport(cl, "f_lme")

tmset_betas_lm = t(mset_betas_lm)

# Output the results
AllLME = pbapply(mset_betas_lm,1,f_lme, Cell_type, Basename, cl = cl)

outdir = paste0("/lustre/projects/Research_Project-191391/valentin/Subtyping/FACS-Celltype/", celltype, "_allcpg.rdata")
save(AllLME, file = outdir)
stopCluster(cl)