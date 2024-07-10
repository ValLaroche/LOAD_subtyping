setwd("D:/valentin/main/")
sources = dir("./scripts/",full.names=TRUE)
for(s in sources){
  source(s)
}
load("FANS_UKBBN/current.Rdata")

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

