setwd("D:/valentin/main/")
library(wateRmelon)
library(sva)
library(matrixStats)
library(CETYGO)
sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}
load("FANS_UKBBN/current.Rdata")
pheno_UKBBN_FANS = read.csv("FANS_UKBBN/sampleSheet.csv")
betas_FANS = readEPIC(idatPath = "FANS_UKBBN/", barcodes = pheno_UKBBN_FANS$X)

outliers <- outlyx(betas_FANS, plot=FALSE)
betas_FANS = betas_FANS[,betas_FANS@phenoData@data$barcode %in% rownames(outliers[outliers$outliers == FALSE,])]
print(outliers) # Keep this

bsc <- bscon(betas_FANS)

msetEPIC_pfilter <- pfilter(betas_FANS)
# 14 samples having 1 % of sites with a detection p-value greater than 0.05 were removed 
# Samples removed: 205854150085_R07C01 205854150085_R08C01 205854150048_R02C01 205854150068_R03C01 205854150068_R07C01 205854150068_R08C01 205854150121_R07C01 205854150103_R06C01 205854150103_R08C01 205854150080_R03C01 205854150076_R04C01 206513950171_R01C01 
# 2607 sites were removed as beadcount <3 in 5 % of samples 
# 18026 sites having 1 % of samples with a detection p-value greater than 0.05 were removed
estimatesex = wateRmelon::estimateSex(betas(msetEPIC_pfilter), do_plot=FALSE)

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

msetEPIC_BMIQ <- BMIQ(msetEPIC_pfilter)

msetEPIC_betas <- betas(msetEPIC_BMIQ)

pwod_bet <- pwod(msetEPIC_betas)

msetEPIC_betas = msetEPIC_betas[-grep("ch.", rownames(msetEPIC_betas)),]
dim(msetEPIC_betas)

pheno_UKBBN_FANS = pheno_UKBBN_FANS[pheno_UKBBN_FANS$X %in% colnames(mset_betas_lm),]
pheno_UKBBN_FANS$Brain.Bank = NA
pheno_UKBBN_FANS$Brain.Bank = pheno_UKBBN[match(pheno_UKBBN_FANS$Individual, pheno_UKBBN$BBNId),]$Brain.Bank
msetEPIC_betas = ComBat(dat=msetEPIC_betas, batch=as.character(pheno_UKBBN_FANS$Brain.Bank), mod=NULL, par.prior=TRUE, prior.plots=FALSE)

# New resid loop with imputation of missing values
mset_betas_lm = msetEPIC_betas
mset_betas_lm = mset_betas_lm[complete.cases(mset_betas_lm),]
for(i in seq(1:nrow(mset_betas_lm))){
  mset_betas_lm[i,] = lm(mset_betas_lm[i,] ~ pheno_UKBBN_FANS$Age + 
                           pheno_UKBBN_FANS$Sex)$residuals
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
mset_betas_lmdf = mset_betas_lmdf[rownames(mset_betas_lmdf) %in% inter_450k_EPIC,]

pheno_UKBBN_FANS$subtype = pheno_UKBBN[match(pheno_UKBBN_FANS$Individual, pheno_UKBBN$BBNId),]$subtype

list_EWAS = list()
for(i in list(c("red","Control"), c("blue", "Control"))){
  for(j in unique(pheno_UKBBN_FANS$Cell_Type)){
    tmp_pheno = pheno_UKBBN_FANS[pheno_UKBBN_FANS$Cell_Type == j,]
    tmp_pheno = tmp_pheno[tmp_pheno$subtype %in% i,]
    tmp_betas = mset_betas_lmdf[,colnames(mset_betas_lmdf) %in% tmp_pheno$X]
    
    tmp_EWAS = limma_DA_FANS(beta_values = tmp_betas,tmp_pheno, 
                    subtype1 = i[1], subtype2 = i[2])
    name = paste0(paste(i, collapse = ""),j)
    list_EWAS[[name]] = tmp_EWAS
  }
}
#Save here
for(i in seq(1:length(list_EWAS))){
  list_EWAS[[i]] = as.data.frame(list_EWAS[[i]])
  list_EWAS[[i]] = list_EWAS[[i]][list_EWAS[[i]]$P.t1 < 5e-2,]
  list_EWAS[[i]] = list_EWAS[[i]][order(list_EWAS[[i]]$P.t1),]
  list_EWAS[[i]]$gene = manifest_EPIC[match(rownames(list_EWAS[[i]]), manifest_EPIC$IlmnID),]$UCSC_RefGene_Name
  print(head(list_EWAS[[i]]))
}

for(i in seq(1:length(list_EWAS))){
  for(j in seq(1:length(list_EWAS))){
    if(i != j){
      print(c(names(list_EWAS)[i],names(list_EWAS)[j]))
      print(length(intersect(rownames(list_EWAS[[i]]),rownames(list_EWAS[[j]]))))
    }
  }
}
inter_mic = intersect(rownames(list_EWAS[["redControlIRF8+"]]),rownames(list_EWAS[["blueControlIRF8+"]]))
list_EWAS[["redControlIRF8+"]][rownames(list_EWAS[["redControlIRF8+"]]) %in% inter_mic,]
list_EWAS[["blueControlIRF8+"]][rownames(list_EWAS[["blueControlIRF8+"]]) %in% inter_mic,]



list_Gometh = list()
for(i in names(list_EWAS)){
  tmp_EWAS = list_EWAS[[i]][list_EWAS[[i]]$P.t1 < 1e-3,]
  list_Gometh[[i]] = gometh(sig.cpg = rownames(tmp_EWAS),
                            all.cpg = rownames(mset_betas_lmdf),
                            collection = c("GO", "KEGG"),
                            array.type = c("EPIC"),
                            sig.genes = FALSE)
  list_Gometh[[i]] = list_Gometh[[i]][list_Gometh[[i]]$P.DE < 0.05,]
  list_Gometh[[i]] = list_Gometh[[i]][order(list_Gometh[[i]]$P.DE),]
}

sig_redC = read.csv("SUBTYPING-PAPER-TABLEFIG/sig_methyl_450k_redC.csv")
sig_blueC = read.csv("SUBTYPING-PAPER-TABLEFIG/sig_methyl_450k_blueC.csv")

intersect(sig_redC$X, rownames(list_EWAS[["redControlNeuN+"]]))
intersect(sig_redC$X, rownames(list_EWAS[["redControlSOX10+"]]))
intersect(sig_redC$X, rownames(list_EWAS[["redControlIRF8+"]]))
intersect(sig_redC$X, rownames(list_EWAS[["redControlNeg"]]))

intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeuN+"]]))
intersect(sig_blueC$X, rownames(list_EWAS[["blueControlSOX10+"]]))
intersect(sig_blueC$X, rownames(list_EWAS[["blueControlIRF8+"]]))
intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeg"]]))

intersect(intersect(sig_redC$X, rownames(list_EWAS[["redControlNeuN+"]])),
          intersect(sig_redC$X, rownames(list_EWAS[["redControlSOX10+"]])))
intersect(intersect(sig_redC$X, rownames(list_EWAS[["redControlNeuN+"]])),
          intersect(sig_redC$X, rownames(list_EWAS[["redControlIRF8+"]])))
intersect(intersect(sig_redC$X, rownames(list_EWAS[["redControlNeuN+"]])),
          intersect(sig_redC$X, rownames(list_EWAS[["redControlNeg"]])))
intersect(intersect(sig_redC$X, rownames(list_EWAS[["redControlSOX10+"]])),
          intersect(sig_redC$X, rownames(list_EWAS[["redControlIRF8+"]])))
intersect(intersect(sig_redC$X, rownames(list_EWAS[["redControlSOX10+"]])),
          intersect(sig_redC$X, rownames(list_EWAS[["redControlNeg"]])))
intersect(intersect(sig_redC$X, rownames(list_EWAS[["redControlIRF8+"]])),
          intersect(sig_redC$X, rownames(list_EWAS[["redControlNeg"]])))

intersect(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeuN+"]])),
          intersect(sig_blueC$X, rownames(list_EWAS[["blueControlSOX10+"]])))
intersect(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeuN+"]])),
          intersect(sig_blueC$X, rownames(list_EWAS[["blueControlIRF8+"]])))
intersect(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeuN+"]])),
          intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeg"]])))
intersect(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlSOX10+"]])),
          intersect(sig_blueC$X, rownames(list_EWAS[["blueControlIRF8+"]])))
intersect(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlSOX10+"]])),
          intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeg"]])))
intersect(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlIRF8+"]])),
          intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeg"]])))

load("450k_cpgs.Rdata")

inter_cpg = length(intersect(sig_redC$X, rownames(list_EWAS[["redControlNeuN+"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_redC), nrow(list_EWAS[["redControlNeuN+"]])) - inter_cpg)
neun_RvC = matrix(c(inter_cpg, #inter 
                    nrow(sig_redC) - inter_cpg, #subtype left 
                    nrow(list_EWAS[["redControlNeuN+"]]) - inter_cpg, #celltype left
                    total_cpg + inter_cpg), #all
                  nrow = 2, ncol = 2)
fisher.test(neun_RvC, alternative = "greater")

inter_cpg = length(intersect(sig_redC$X, rownames(list_EWAS[["redControlSOX10+"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_redC), nrow(list_EWAS[["redControlSOX10+"]])) - inter_cpg)
olig_RvC = matrix(c(inter_cpg, #inter 
                    nrow(sig_redC) - inter_cpg, #subtype left 
                    nrow(list_EWAS[["redControlSOX10+"]]) - inter_cpg, #celltype left
                    total_cpg + inter_cpg), #all
                  nrow = 2, ncol = 2)
fisher.test(olig_RvC, alternative = "greater")

inter_cpg = length(intersect(sig_redC$X, rownames(list_EWAS[["redControlIRF8+"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_redC), nrow(list_EWAS[["redControlIRF8+"]])))
mic_RvC = matrix(c(inter_cpg, #inter 
                   nrow(sig_redC) - inter_cpg, #subtype left 
                   nrow(list_EWAS[["redControlIRF8+"]]) - inter_cpg, #celltype left
                   total_cpg + inter_cpg), #all
                 nrow = 2, ncol = 2)
fisher.test(mic_RvC, alternative = "greater")


inter_cpg = length(intersect(sig_redC$X, rownames(list_EWAS[["redControlNeg"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_redC), nrow(list_EWAS[["redControlNeg"]])) - inter_cpg)
ast_RvC = matrix(c(inter_cpg, #inter 
                   nrow(sig_redC) - inter_cpg, #subtype left 
                   nrow(list_EWAS[["redControlNeg"]]) - inter_cpg, #celltype left
                   total_cpg + inter_cpg), #all
                 nrow = 2, ncol = 2)
fisher.test(ast_RvC, alternative = "greater")

inter_cpg = length(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeuN+"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_blueC), nrow(list_EWAS[["blueControlNeuN+"]])) - inter_cpg)
neun_BvC = matrix(c(inter_cpg, #inter 
                    nrow(sig_blueC) - inter_cpg, #subtype left 
                    nrow(list_EWAS[["blueControlNeuN+"]]) - inter_cpg, #celltype left
                    total_cpg + inter_cpg), #all
                  nrow = 2, ncol = 2)
fisher.test(neun_BvC, alternative = "greater")

inter_cpg = length(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlSOX10+"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_blueC), nrow(list_EWAS[["blueControlSOX10+"]])) - inter_cpg)
olig_BvC = matrix(c(inter_cpg, #inter 
                    nrow(sig_blueC) - inter_cpg, #subtype left 
                    nrow(list_EWAS[["blueControlSOX10+"]]) - inter_cpg, #celltype left
                    total_cpg + inter_cpg), #all
                  nrow = 2, ncol = 2)
fisher.test(olig_BvC, alternative = "greater")

inter_cpg = length(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlIRF8+"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_blueC), nrow(list_EWAS[["blueControlIRF8+"]])) - inter_cpg)
mic_BvC = matrix(c(inter_cpg, #inter 
                   nrow(sig_blueC) - inter_cpg, #subtype left 
                   nrow(list_EWAS[["blueControlIRF8+"]]) - inter_cpg, #celltype left
                   total_cpg + inter_cpg), #all
                 nrow = 2, ncol = 2)
fisher.test(mic_BvC, alternative = "greater")

inter_cpg = length(intersect(sig_blueC$X, rownames(list_EWAS[["blueControlNeg"]])))
total_cpg = length(inter_450k_EPIC) - (sum(nrow(sig_blueC), nrow(list_EWAS[["blueControlNeg"]])) - inter_cpg)
ast_BvC = matrix(c(inter_cpg, #inter 
                   nrow(sig_blueC) - inter_cpg, #subtype left 
                   nrow(list_EWAS[["blueControlNeg"]]) - inter_cpg, #celltype left
                   total_cpg + inter_cpg), #all
                 nrow = 2, ncol = 2)
fisher.test(ast_BvC)

