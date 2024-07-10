WGCNArest <- function(){
  setwd("D:/valentin/main/")
# load("./WGCNA/UKBBN/UKBBN_subtyping.wgcna.network.rdat")
# load("WGCNA/UKBBN/UKBBN_subtyping.wgcna.softThreshold.rdat")
# load("WGCNA/UKBBN/UKBBN_subtyping.wgcna1.rdat")
# load("WGCNA/PITT/PITT_toBDR_preservation_5000parallel_ModulePreservation.Rdata")
# load("WGCNA/to_WGCNA.Rdata")
# load("savepitt.Rdata")
load("WGCNA/current_WGCNA_analysis.Rdata")  

WGCNA_network = net
table(WGCNA_network$colors) ## label 0 is reserved for genes outside of all modules

WGCNA_labels = WGCNA_network$colors

table(WGCNA_labels)
options(stringsAsFactors = FALSE)

##### BDR #####

WGCNA_eigengenesBDR = moduleEigengenes(t(data_matrix_AD_BDR), WGCNA_labels)$eigengenes

WGCNA_phenoBDR = pheno_BDR[pheno_BDR$HC_clusters != "Control",]

#####

WGCNA_traitsHC = as.data.frame(cbind(WGCNA_phenoBDR$HC_clusters,WGCNA_phenoBDR$HC_clusters,WGCNA_phenoBDR$HC_clusters))
rownames(WGCNA_traitsHC) = WGCNA_phenoBDR$Basename
colnames(WGCNA_traitsHC) = c("sub1","sub2","sub3")

WGCNA_traitsHC[WGCNA_traitsHC$sub1 == 1,]$sub1 = "1"
WGCNA_traitsHC[WGCNA_traitsHC$sub1 != 1,]$sub1 = "0"
WGCNA_traitsHC[WGCNA_traitsHC$sub2 != 2,]$sub2 = "0"
WGCNA_traitsHC[WGCNA_traitsHC$sub2 == 2,]$sub2 = "1"
WGCNA_traitsHC[WGCNA_traitsHC$sub3 != 3,]$sub3 = "0"
WGCNA_traitsHC[WGCNA_traitsHC$sub3 == 3,]$sub3 = "1"

#####

WGCNA_traitsKM = as.data.frame(cbind(WGCNA_phenoBDR$KM_clusters,WGCNA_phenoBDR$KM_clusters,WGCNA_phenoBDR$KM_clusters))
rownames(WGCNA_traitsKM) = WGCNA_phenoBDR$Basename
colnames(WGCNA_traitsKM) = c("sub1","sub2","sub3")

WGCNA_traitsKM[WGCNA_traitsKM$sub1 == 1,]$sub1 = "1"
WGCNA_traitsKM[WGCNA_traitsKM$sub1 != 1,]$sub1 = "0"
WGCNA_traitsKM[WGCNA_traitsKM$sub2 != 2,]$sub2 = "0"
WGCNA_traitsKM[WGCNA_traitsKM$sub2 == 2,]$sub2 = "1"
WGCNA_traitsKM[WGCNA_traitsKM$sub3 != 3,]$sub3 = "0"
WGCNA_traitsKM[WGCNA_traitsKM$sub3 == 3,]$sub3 = "1"

#####
WGCNA_traitsPAM = as.data.frame(cbind(WGCNA_phenoBDR$PAM_clusters,WGCNA_phenoBDR$PAM_clusters,WGCNA_phenoBDR$PAM_clusters))
rownames(WGCNA_traitsPAM) = WGCNA_phenoBDR$Basename
colnames(WGCNA_traitsPAM) = c("sub1","sub2","sub3")

WGCNA_traitsPAM[WGCNA_traitsPAM$sub1 == 1,]$sub1 = "1"
WGCNA_traitsPAM[WGCNA_traitsPAM$sub1 != 1,]$sub1 = "0"
WGCNA_traitsPAM[WGCNA_traitsPAM$sub2 != 2,]$sub2 = "0"
WGCNA_traitsPAM[WGCNA_traitsPAM$sub2 == 2,]$sub2 = "1"
WGCNA_traitsPAM[WGCNA_traitsPAM$sub3 != 3,]$sub3 = "0"
WGCNA_traitsPAM[WGCNA_traitsPAM$sub3 == 3,]$sub3 = "1"

#####
WGCNA_traitsSpec = as.data.frame(cbind(WGCNA_phenoBDR$Spec_clusters,WGCNA_phenoBDR$Spec_clusters,WGCNA_phenoBDR$Spec_clusters))
rownames(WGCNA_traitsSpec) = WGCNA_phenoBDR$Basename
colnames(WGCNA_traitsSpec) = c("sub1","sub2","sub3")

WGCNA_traitsSpec[WGCNA_traitsSpec$sub1 == 1,]$sub1 = "1"
WGCNA_traitsSpec[WGCNA_traitsSpec$sub1 != 1,]$sub1 = "0"
WGCNA_traitsSpec[WGCNA_traitsSpec$sub2 != 2,]$sub2 = "0"
WGCNA_traitsSpec[WGCNA_traitsSpec$sub2 == 2,]$sub2 = "1"
WGCNA_traitsSpec[WGCNA_traitsSpec$sub3 != 3,]$sub3 = "0"
WGCNA_traitsSpec[WGCNA_traitsSpec$sub3 == 3,]$sub3 = "1"

#####

WGCNA_traitsCM = as.data.frame(cbind(WGCNA_phenoBDR$CM_clusters,WGCNA_phenoBDR$CM_clusters,WGCNA_phenoBDR$CM_clusters,WGCNA_phenoBDR$CM_clusters))
rownames(WGCNA_traitsCM) = WGCNA_phenoBDR$Basename
colnames(WGCNA_traitsCM) = c("sub1","sub2","sub3","excluded")

WGCNA_traitsCM[WGCNA_traitsCM$sub1 != "Cluster 1",]$sub1 = "0"
WGCNA_traitsCM[WGCNA_traitsCM$sub1 == "Cluster 1",]$sub1 = "1"
WGCNA_traitsCM[WGCNA_traitsCM$sub2 != "Cluster 2",]$sub2 = "0"
WGCNA_traitsCM[WGCNA_traitsCM$sub2 == "Cluster 2",]$sub2 = "1"
WGCNA_traitsCM[WGCNA_traitsCM$sub3 != "Cluster 3",]$sub3 = "0"
WGCNA_traitsCM[WGCNA_traitsCM$sub3 == "Cluster 3",]$sub3 = "1"
WGCNA_traitsCM[WGCNA_traitsCM$excluded != "Excluded",]$excluded = "0"
WGCNA_traitsCM[WGCNA_traitsCM$excluded == "Excluded",]$excluded = "1"
#####

WGCNA_MEorder = orderMEs(WGCNA_eigengenesBDR)
moduleTraitCor <- cor(WGCNA_MEorder, WGCNA_traitsCM, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(WGCNA_traitsHC))

## visualize module-trait assocations in color-coded table (color=coded by correlation value):
#sizeGrWindow(10,6)
## will display correlations and their p values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
#
pdf(file = "WGCNA/CM_BDR_modulecor.pdf", width = 15, height = 10)
par(mar = c(6, 8.5, 3, 3))
## display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#####

WGCNA_toanovaBDR = cbind(WGCNA_phenoBDR$HC_clusters, WGCNA_eigengenesBDR)
colnames(WGCNA_toanovaBDR)[1] = "clusters"

WGCNA_toanovaBDR = WGCNA_toanovaBDR[WGCNA_toanovaBDR[,1] != "Excluded",]

tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaBDR)))){
  tmptukey_table <- aov(WGCNA_toanovaBDR[,i] ~ clusters, data = WGCNA_toanovaBDR) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanovaBDR[,-1])*1.75, 
                                          max(WGCNA_toanovaBDR[,-1])*2, 
                                          max(WGCNA_toanovaBDR[,-1])*1.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaBDR)[c(2:17)], each = 3)

WGCNA_toanovaBDR$Basename = rownames(WGCNA_toanovaBDR)
WGCNA_toanovaBDRmelt = melt(WGCNA_toanovaBDR, id.vars = c("Basename","clusters"))

pdf(file = "WGCNA/BDR/HC_UKBBNnet_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanovaBDRmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free") +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanovaBDR[, i]*1.2, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj", tip.length = 0.03, bracket.shorten = 0.05) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside")
dev.off()

##### UKBBN #####

WGCNA_eigengenesUKBBN = moduleEigengenes(t(data_matrix_AD_UKBBN), WGCNA_labels)$eigengenes
WGCNA_MEorderUKBBN = orderMEs(WGCNA_eigengenesUKBBN)

WGCNA_phenoUKBBN = pheno_UKBBN[pheno_UKBBN$HC_clusters != "Control",]



#####
WGCNA_traitsUKBBNHC = as.data.frame(cbind(WGCNA_phenoUKBBN$HC_clusters,WGCNA_phenoUKBBN$HC_clusters,WGCNA_phenoUKBBN$HC_clusters))
rownames(WGCNA_traitsUKBBNHC) = WGCNA_phenoUKBBN$Basename
colnames(WGCNA_traitsUKBBNHC) = c("sub1","sub2","sub3")

WGCNA_traitsUKBBNHC[WGCNA_traitsUKBBNHC$sub1 != "1",]$sub1 = "0"
WGCNA_traitsUKBBNHC[WGCNA_traitsUKBBNHC$sub1 == "1",]$sub1 = "1"
WGCNA_traitsUKBBNHC[WGCNA_traitsUKBBNHC$sub2 != "2",]$sub2 = "0"
WGCNA_traitsUKBBNHC[WGCNA_traitsUKBBNHC$sub2 == "2",]$sub2 = "1"
WGCNA_traitsUKBBNHC[WGCNA_traitsUKBBNHC$sub3 != "3",]$sub3 = "0"
WGCNA_traitsUKBBNHC[WGCNA_traitsUKBBNHC$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsUKBBNKM = as.data.frame(cbind(WGCNA_phenoUKBBN$KM_clusters,WGCNA_phenoUKBBN$KM_clusters,WGCNA_phenoUKBBN$KM_clusters))
rownames(WGCNA_traitsUKBBNKM) = WGCNA_phenoUKBBN$Basename
colnames(WGCNA_traitsUKBBNKM) = c("sub1","sub2","sub3")

WGCNA_traitsUKBBNKM[WGCNA_traitsUKBBNKM$sub1 != "1",]$sub1 = "0"
WGCNA_traitsUKBBNKM[WGCNA_traitsUKBBNKM$sub1 == "1",]$sub1 = "1"
WGCNA_traitsUKBBNKM[WGCNA_traitsUKBBNKM$sub2 != "2",]$sub2 = "0"
WGCNA_traitsUKBBNKM[WGCNA_traitsUKBBNKM$sub2 == "2",]$sub2 = "1"
WGCNA_traitsUKBBNKM[WGCNA_traitsUKBBNKM$sub3 != "3",]$sub3 = "0"
WGCNA_traitsUKBBNKM[WGCNA_traitsUKBBNKM$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsUKBBNPAM = as.data.frame(cbind(WGCNA_phenoUKBBN$PAM_clusters,WGCNA_phenoUKBBN$PAM_clusters,WGCNA_phenoUKBBN$PAM_clusters))
rownames(WGCNA_traitsUKBBNPAM) = WGCNA_phenoUKBBN$Basename
colnames(WGCNA_traitsUKBBNPAM) = c("sub1","sub2","sub3")

WGCNA_traitsUKBBNPAM[WGCNA_traitsUKBBNPAM$sub1 != "1",]$sub1 = "0"
WGCNA_traitsUKBBNPAM[WGCNA_traitsUKBBNPAM$sub1 == "1",]$sub1 = "1"
WGCNA_traitsUKBBNPAM[WGCNA_traitsUKBBNPAM$sub2 != "2",]$sub2 = "0"
WGCNA_traitsUKBBNPAM[WGCNA_traitsUKBBNPAM$sub2 == "2",]$sub2 = "1"
WGCNA_traitsUKBBNPAM[WGCNA_traitsUKBBNPAM$sub3 != "3",]$sub3 = "0"
WGCNA_traitsUKBBNPAM[WGCNA_traitsUKBBNPAM$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsUKBBNSpec = as.data.frame(cbind(WGCNA_phenoUKBBN$Spec_clusters,WGCNA_phenoUKBBN$Spec_clusters,WGCNA_phenoUKBBN$Spec_clusters))
rownames(WGCNA_traitsUKBBNSpec) = WGCNA_phenoUKBBN$Basename
colnames(WGCNA_traitsUKBBNSpec) = c("sub1","sub2","sub3")

WGCNA_traitsUKBBNSpec[WGCNA_traitsUKBBNSpec$sub1 != "1",]$sub1 = "0"
WGCNA_traitsUKBBNSpec[WGCNA_traitsUKBBNSpec$sub1 == "1",]$sub1 = "1"
WGCNA_traitsUKBBNSpec[WGCNA_traitsUKBBNSpec$sub2 != "2",]$sub2 = "0"
WGCNA_traitsUKBBNSpec[WGCNA_traitsUKBBNSpec$sub2 == "2",]$sub2 = "1"
WGCNA_traitsUKBBNSpec[WGCNA_traitsUKBBNSpec$sub3 != "3",]$sub3 = "0"
WGCNA_traitsUKBBNSpec[WGCNA_traitsUKBBNSpec$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsUKBBNCM = as.data.frame(cbind(WGCNA_phenoUKBBN$CM_clusters,WGCNA_phenoUKBBN$CM_clusters,WGCNA_phenoUKBBN$CM_clusters, WGCNA_phenoUKBBN$CM_clusters))
rownames(WGCNA_traitsUKBBNCM) = WGCNA_phenoUKBBN$Basename
colnames(WGCNA_traitsUKBBNCM) = c("sub1","sub2","sub3", "excluded")

WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$sub1 != "Cluster 1",]$sub1 = "0"
WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$sub1 == "Cluster 1",]$sub1 = "1"
WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$sub2 != "Cluster 2",]$sub2 = "0"
WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$sub2 == "Cluster 2",]$sub2 = "1"
WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$sub3 != "Cluster 3",]$sub3 = "0"
WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$sub3 == "Cluster 3",]$sub3 = "1"
WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$excluded != "Excluded",]$excluded = "0"
WGCNA_traitsUKBBNCM[WGCNA_traitsUKBBNCM$excluded == "Excluded",]$excluded = "1"

#####

moduleTraitCor <- cor(WGCNA_MEorderUKBBN, WGCNA_traitsUKBBNHC, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(WGCNA_traitsUKBBNHC))

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
#
pdf(file = "WGCNA/HC_UKBBN_modulecor.pdf", width = 15, height = 10)
par(mar = c(6, 8.5, 3, 3))
## display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#####

WGCNA_toanovaUKBBN = cbind(WGCNA_phenoUKBBN$HC_clusters, WGCNA_eigengenesUKBBN)
colnames(WGCNA_toanovaUKBBN)[1] = "clusters"

WGCNA_toanovaUKBBN = WGCNA_toanovaUKBBN[WGCNA_toanovaUKBBN[,1] != "Excluded",]


tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaUKBBN)))){
  tmptukey_table <- aov(WGCNA_toanovaUKBBN[,i] ~ clusters, data = WGCNA_toanovaUKBBN) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanovaUKBBN[,-1])*1.75, 
                                          max(WGCNA_toanovaUKBBN[,-1])*2, 
                                          max(WGCNA_toanovaUKBBN[,-1])*1.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaUKBBN)[c(2:8)], each = 3)

WGCNA_toanovaUKBBN$Basename = rownames(WGCNA_toanovaUKBBN)
WGCNA_toanovaUKBBNmelt = melt(WGCNA_toanovaUKBBN, id.vars = c("Basename","clusters"))

labels = paste0(rownames(stats), "  ", signif(stats$Zsummary.pres, digits = 4))
names(labels) = paste0("ME",rownames(stats))

pdf(file = "WGCNA/UKBBN/HC_PITTnet_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanovaUKBBNmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free", 
             labeller = labeller(variable = labels)) +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanovaUKBBN[, i]*1.2, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj", tip.length = 0.03, bracket.shorten = 0.05) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside")
dev.off()


#####

WGCNA_toanovaBDR = cbind(WGCNA_phenoBDR$HC_clusters, WGCNA_eigengenesBDR)
colnames(WGCNA_toanovaBDR)[1] = "clusters"
WGCNA_toanovaBDR$clusters = paste0("BDR", WGCNA_toanovaBDR$clusters)


WGCNA_toanovaUKBBN = cbind(WGCNA_phenoUKBBN$HC_clusters, WGCNA_eigengenesUKBBN)
colnames(WGCNA_toanovaUKBBN)[1] = "clusters"
WGCNA_toanovaUKBBN$clusters = paste0("UKBBN", WGCNA_toanovaUKBBN$clusters)

WGCNA_toanova_all = rbind(WGCNA_toanovaBDR, WGCNA_toanovaUKBBN)


tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaUKBBN)))){
  tmptukey_table <- aov(WGCNA_toanova_all[,i] ~ clusters, data = WGCNA_toanova_all) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*2, 
                                          max(WGCNA_toanova_all[,-1])*2.1, 
                                          max(WGCNA_toanova_all[,-1])*2.2,
                                          max(WGCNA_toanova_all[,-1])*2.3, 
                                          max(WGCNA_toanova_all[,-1])*2.4, 
                                          max(WGCNA_toanova_all[,-1])*2.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaUKBBN)[c(2:14)], each = 15)


labels = paste0(rownames(stats), "  MP=", signif(stats$Zsummary.pres, digits = 4))
names(labels) = paste0("ME",rownames(stats))

WGCNA_toanova_allmelt = melt(WGCNA_toanova_all, id.vars = "clusters")
WGCNA_toanova_allmelt$clusters = factor(WGCNA_toanova_allmelt$clusters,
                                        levels = c("BDR1","UKBBN1","BDR2","UKBBN2","BDR3","UKBBN3"))

pdf(file = "WGCNA/UKBBN/HC_BDRUKBBN_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanova_allmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free", 
             labeller = labeller(variable = labels)) +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanova_allmelt$value*1, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
        axis.text.x=element_text(size=7))
dev.off()


##### PITT #####

WGCNA_eigengenesPITT = moduleEigengenes(t(data_matrix_AD_PITT), WGCNA_labels)$eigengenes
WGCNA_MEorderPITT = orderMEs(WGCNA_eigengenesPITT)

WGCNA_phenoPITT = pheno_PITT[pheno_PITT$HC_clusters != "Control",]

#####
WGCNA_traitsPITTHC = as.data.frame(cbind(WGCNA_phenoPITT$HC_clusters,WGCNA_phenoPITT$HC_clusters,WGCNA_phenoPITT$HC_clusters))
rownames(WGCNA_traitsPITTHC) = WGCNA_phenoPITT$Basename
colnames(WGCNA_traitsPITTHC) = c("sub1","sub2","sub3")

WGCNA_traitsPITTHC[WGCNA_traitsPITTHC$sub1 != "1",]$sub1 = "0"
WGCNA_traitsPITTHC[WGCNA_traitsPITTHC$sub1 == "1",]$sub1 = "1"
WGCNA_traitsPITTHC[WGCNA_traitsPITTHC$sub2 != "2",]$sub2 = "0"
WGCNA_traitsPITTHC[WGCNA_traitsPITTHC$sub2 == "2",]$sub2 = "1"
WGCNA_traitsPITTHC[WGCNA_traitsPITTHC$sub3 != "3",]$sub3 = "0"
WGCNA_traitsPITTHC[WGCNA_traitsPITTHC$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsPITTKM = as.data.frame(cbind(WGCNA_phenoPITT$KM_clusters,WGCNA_phenoPITT$KM_clusters,WGCNA_phenoPITT$KM_clusters))
rownames(WGCNA_traitsPITTKM) = WGCNA_phenoPITT$Basename
colnames(WGCNA_traitsPITTKM) = c("sub1","sub2","sub3")

WGCNA_traitsPITTKM[WGCNA_traitsPITTKM$sub1 != "1",]$sub1 = "0"
WGCNA_traitsPITTKM[WGCNA_traitsPITTKM$sub1 == "1",]$sub1 = "1"
WGCNA_traitsPITTKM[WGCNA_traitsPITTKM$sub2 != "2",]$sub2 = "0"
WGCNA_traitsPITTKM[WGCNA_traitsPITTKM$sub2 == "2",]$sub2 = "1"
WGCNA_traitsPITTKM[WGCNA_traitsPITTKM$sub3 != "3",]$sub3 = "0"
WGCNA_traitsPITTKM[WGCNA_traitsPITTKM$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsPITTPAM = as.data.frame(cbind(WGCNA_phenoPITT$PAM_clusters,WGCNA_phenoPITT$PAM_clusters,WGCNA_phenoPITT$PAM_clusters))
rownames(WGCNA_traitsPITTPAM) = WGCNA_phenoPITT$Basename
colnames(WGCNA_traitsPITTPAM) = c("sub1","sub2","sub3")

WGCNA_traitsPITTPAM[WGCNA_traitsPITTPAM$sub1 != "1",]$sub1 = "0"
WGCNA_traitsPITTPAM[WGCNA_traitsPITTPAM$sub1 == "1",]$sub1 = "1"
WGCNA_traitsPITTPAM[WGCNA_traitsPITTPAM$sub2 != "2",]$sub2 = "0"
WGCNA_traitsPITTPAM[WGCNA_traitsPITTPAM$sub2 == "2",]$sub2 = "1"
WGCNA_traitsPITTPAM[WGCNA_traitsPITTPAM$sub3 != "3",]$sub3 = "0"
WGCNA_traitsPITTPAM[WGCNA_traitsPITTPAM$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsPITTSpec = as.data.frame(cbind(WGCNA_phenoPITT$Spec_clusters,WGCNA_phenoPITT$Spec_clusters,WGCNA_phenoPITT$Spec_clusters))
rownames(WGCNA_traitsPITTSpec) = WGCNA_phenoPITT$Basename
colnames(WGCNA_traitsPITTSpec) = c("sub1","sub2","sub3")

WGCNA_traitsPITTSpec[WGCNA_traitsPITTSpec$sub1 != "1",]$sub1 = "0"
WGCNA_traitsPITTSpec[WGCNA_traitsPITTSpec$sub1 == "1",]$sub1 = "1"
WGCNA_traitsPITTSpec[WGCNA_traitsPITTSpec$sub2 != "2",]$sub2 = "0"
WGCNA_traitsPITTSpec[WGCNA_traitsPITTSpec$sub2 == "2",]$sub2 = "1"
WGCNA_traitsPITTSpec[WGCNA_traitsPITTSpec$sub3 != "3",]$sub3 = "0"
WGCNA_traitsPITTSpec[WGCNA_traitsPITTSpec$sub3 == "3",]$sub3 = "1"

#####

WGCNA_traitsPITTCM = as.data.frame(cbind(WGCNA_phenoPITT$CM_clusters,WGCNA_phenoPITT$CM_clusters,WGCNA_phenoPITT$CM_clusters, WGCNA_phenoPITT$CM_clusters))
rownames(WGCNA_traitsPITTCM) = WGCNA_phenoPITT$Basename
colnames(WGCNA_traitsPITTCM) = c("sub1","sub2","sub3", "excluded")

WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$sub1 != "Cluster 1",]$sub1 = "0"
WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$sub1 == "Cluster 1",]$sub1 = "1"
WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$sub2 != "Cluster 2",]$sub2 = "0"
WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$sub2 == "Cluster 2",]$sub2 = "1"
WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$sub3 != "Cluster 3",]$sub3 = "0"
WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$sub3 == "Cluster 3",]$sub3 = "1"
WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$excluded != "Excluded",]$excluded = "0"
WGCNA_traitsPITTCM[WGCNA_traitsPITTCM$excluded == "Excluded",]$excluded = "1"

#####

moduleTraitCor <- cor(WGCNA_MEorderPITT, WGCNA_traitsPITTHC, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(WGCNA_traitsPITTHC))

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
#
pdf(file = "WGCNA/HC_PITT_modulecor.pdf", width = 15, height = 10)
par(mar = c(6, 8.5, 3, 3))
## display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#####

WGCNA_toanovaPITT = cbind(WGCNA_phenoPITT$HC_clusters, WGCNA_eigengenesPITT)
colnames(WGCNA_toanovaPITT)[1] = "clusters"

WGCNA_toanovaPITT = WGCNA_toanovaPITT[WGCNA_toanovaPITT[,1] != "Excluded",]

tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaPITT)))){
  tmptukey_table <- aov(WGCNA_toanovaPITT[,i] ~ clusters, data = WGCNA_toanovaPITT) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanovaPITT[,-1])*1.75, 
                                          max(WGCNA_toanovaPITT[,-1])*2, 
                                          max(WGCNA_toanovaPITT[,-1])*1.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaPITT)[c(2:17)], each = 3)

WGCNA_toanovaPITT$Basename = rownames(WGCNA_toanovaPITT)
WGCNA_toanovaPITTmelt = melt(WGCNA_toanovaPITT, id.vars = c("Basename","clusters"))

pdf(file = "WGCNA/PITT/HC_UKBBNnet_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanovaPITTmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free") +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanovaPITT[, i]*1.2, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj", tip.length = 0.03, bracket.shorten = 0.05) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside")
dev.off()


#####

WGCNA_toanovaBDR = cbind(WGCNA_phenoBDR$HC_clusters, WGCNA_eigengenesBDR)
colnames(WGCNA_toanovaBDR)[1] = "clusters"
WGCNA_toanovaBDR$clusters = paste0("BDR", WGCNA_toanovaBDR$clusters)


WGCNA_toanovaUKBBN = cbind(WGCNA_phenoUKBBN$HC_clusters, WGCNA_eigengenesUKBBN)
colnames(WGCNA_toanovaUKBBN)[1] = "clusters"
WGCNA_toanovaUKBBN$clusters = paste0("UKBBN", WGCNA_toanovaUKBBN$clusters)

WGCNA_toanova_all = rbind(WGCNA_toanovaBDR, WGCNA_toanovaUKBBN)


tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaUKBBN)))){
  tmptukey_table <- aov(WGCNA_toanova_all[,i] ~ clusters, data = WGCNA_toanova_all) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*2, 
                                          max(WGCNA_toanova_all[,-1])*2.1, 
                                          max(WGCNA_toanova_all[,-1])*2.2,
                                          max(WGCNA_toanova_all[,-1])*2.3, 
                                          max(WGCNA_toanova_all[,-1])*2.4, 
                                          max(WGCNA_toanova_all[,-1])*2.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaUKBBN)[c(2:14)], each = 15)


labels = paste0(rownames(stats), "  MP=", signif(stats$Zsummary.pres, digits = 4))
names(labels) = paste0("ME",rownames(stats))

WGCNA_toanova_allmelt = melt(WGCNA_toanova_all, id.vars = "clusters")
WGCNA_toanova_allmelt$clusters = factor(WGCNA_toanova_allmelt$clusters,
                                        levels = c("BDR1","UKBBN1","BDR2","UKBBN2","BDR3","UKBBN3"))

pdf(file = "WGCNA/UKBBN/HC_BDRUKBBN_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanova_allmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free", 
             labeller = labeller(variable = labels)) +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanova_allmelt$value*1, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
        axis.text.x=element_text(size=7))
dev.off()


##### Loadings #####
# HC

models_network_UKBBN = list(list(),list(),list(),list(),list())

splsdaBDR_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesBDR, WGCNA_phenoBDR$HC_clusters, ncomp = 6)
# plotIndiv(splsdaBDR_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadBDR_UKBBNnet = data.frame(plotLoadings(splsdaBDR_UKBBNnet, comp = 1, contrib = "max", legend = FALSE))
loadBDR_UKBBNnet$module = as.factor(rownames(loadBDR_UKBBNnet))

loadBDR_UKBBNnet$cohort = "BDR"

models_network_UKBBN[[1]][[1]] = splsdaBDR_UKBBNnet

splsdaUKBBN_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesUKBBN, WGCNA_phenoUKBBN$HC_clusters, ncomp = 6)
# plotIndiv(splsdaUKBBN_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadUKBBN_UKBBNnet = data.frame(plotLoadings(splsdaUKBBN_UKBBNnet, comp = 1, contrib = "max", size.legend = 0.7, title = "HCUKBBN"))
loadUKBBN_UKBBNnet$module = as.factor(rownames(loadUKBBN_UKBBNnet))

loadUKBBN_UKBBNnet$cohort = "UKBBN"

models_network_UKBBN[[1]][[2]] = splsdaUKBBN_UKBBNnet


merge_loadings_UKBBN = rbind(loadBDR_UKBBNnet, loadUKBBN_UKBBNnet)

ggplot(merge_loadings_UKBBN, aes(x = X.importance, y = module, fill = X.GroupContrib)) +
  ggtitle("HC_BDRtoUKBBN") +
  geom_bar(stat = "identity") +
  xlab("Loadings") + ylab("Module") +
  facet_wrap(~cohort) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))


# KM
splsdaBDR_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesBDR, WGCNA_phenoBDR$KM_clusters, ncomp = 6)
# plotIndiv(splsdaBDR_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadBDR_UKBBNnet = data.frame(plotLoadings(splsdaBDR_UKBBNnet, comp = 1, contrib = "max", legend = FALSE))
loadBDR_UKBBNnet$module = as.factor(rownames(loadBDR_UKBBNnet))

loadBDR_UKBBNnet$cohort = "BDR"

models_network_UKBBN[[2]][[1]] = splsdaBDR_UKBBNnet

splsdaUKBBN_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesUKBBN, WGCNA_phenoUKBBN$KM_clusters, ncomp = 6)
# plotIndiv(splsdaUKBBN_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadUKBBN_UKBBNnet = data.frame(plotLoadings(splsdaUKBBN_UKBBNnet, comp = 1, contrib = "max", size.legend = 0.7, title = "KMUKBBN"))
loadUKBBN_UKBBNnet$module = as.factor(rownames(loadUKBBN_UKBBNnet))

loadUKBBN_UKBBNnet$cohort = "UKBBN"

models_network_UKBBN[[2]][[2]] = splsdaUKBBN_UKBBNnet


merge_loadings_UKBBN = rbind(loadBDR_UKBBNnet, loadUKBBN_UKBBNnet)

ggplot(merge_loadings_UKBBN, aes(x = X.importance, y = module, fill = X.GroupContrib)) +
  ggtitle("KM_BDRtoUKBBN") +
  geom_bar(stat = "identity") +
  xlab("Loadings") + ylab("Module") +
  facet_wrap(~cohort) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

#PAM

splsdaBDR_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesBDR, WGCNA_phenoBDR$PAM_clusters, ncomp = 6)
# plotIndiv(splsdaBDR_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadBDR_UKBBNnet = data.frame(plotLoadings(splsdaBDR_UKBBNnet, comp = 1, contrib = "max", legend = FALSE))
loadBDR_UKBBNnet$module = as.factor(rownames(loadBDR_UKBBNnet))

loadBDR_UKBBNnet$cohort = "BDR"

models_network_UKBBN[[3]][[1]] = splsdaBDR_UKBBNnet

splsdaUKBBN_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesUKBBN, WGCNA_phenoUKBBN$PAM_clusters, ncomp = 6)
# plotIndiv(splsdaUKBBN_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadUKBBN_UKBBNnet = data.frame(plotLoadings(splsdaUKBBN_UKBBNnet, comp = 1, contrib = "max", size.legend = 0.7, title = "PAMUKBBN"))
loadUKBBN_UKBBNnet$module = as.factor(rownames(loadUKBBN_UKBBNnet))

loadUKBBN_UKBBNnet$cohort = "UKBBN"

models_network_UKBBN[[3]][[2]] = splsdaUKBBN_UKBBNnet


merge_loadings_UKBBN = rbind(loadBDR_UKBBNnet, loadUKBBN_UKBBNnet)

ggplot(merge_loadings_UKBBN, aes(x = X.importance, y = module, fill = X.GroupContrib)) +
  ggtitle("PAM_BDRtoUKBBN") +
  geom_bar(stat = "identity") +
  xlab("Loadings") + ylab("Module") +
  facet_wrap(~cohort) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# Spec

splsdaBDR_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesBDR, WGCNA_phenoBDR$Spec_clusters, ncomp = 6)
# plotIndiv(splsdaBDR_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadBDR_UKBBNnet = data.frame(plotLoadings(splsdaBDR_UKBBNnet, comp = 1, contrib = "max", legend = FALSE))
loadBDR_UKBBNnet$module = as.factor(rownames(loadBDR_UKBBNnet))

loadBDR_UKBBNnet$cohort = "BDR"

models_network_UKBBN[[4]][[1]] = splsdaBDR_UKBBNnet

splsdaUKBBN_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesUKBBN, WGCNA_phenoUKBBN$Spec_clusters, ncomp = 6)
# plotIndiv(splsdaUKBBN_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadUKBBN_UKBBNnet = data.frame(plotLoadings(splsdaUKBBN_UKBBNnet, comp = 1, contrib = "max", size.legend = 0.7, title = "SpecUKBBN"))
loadUKBBN_UKBBNnet$module = as.factor(rownames(loadUKBBN_UKBBNnet))

loadUKBBN_UKBBNnet$cohort = "UKBBN"

models_network_UKBBN[[4]][[2]] = splsdaUKBBN_UKBBNnet


merge_loadings_UKBBN = rbind(loadBDR_UKBBNnet, loadUKBBN_UKBBNnet)

ggplot(merge_loadings_UKBBN, aes(x = X.importance, y = module, fill = X.GroupContrib)) +
  ggtitle("Spec_BDRtoUKBBN") +
  geom_bar(stat = "identity") +
  xlab("Loadings") + ylab("Module") +
  facet_wrap(~cohort) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# CM

splsdaBDR_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesBDR, WGCNA_phenoBDR$CM_clusters, ncomp = 6)
# plotIndiv(splsdaBDR_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadBDR_UKBBNnet = data.frame(plotLoadings(splsdaBDR_UKBBNnet, comp = 1, contrib = "max", legend = FALSE))
loadBDR_UKBBNnet$module = as.factor(rownames(loadBDR_UKBBNnet))

loadBDR_UKBBNnet$cohort = "BDR"

models_network_UKBBN[[5]][[1]] = splsdaBDR_UKBBNnet

splsdaUKBBN_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesUKBBN, WGCNA_phenoUKBBN$CM_clusters, ncomp = 6)
# plotIndiv(splsdaUKBBN_UKBBNnet,ind.names = FALSE, comp = c(1,2), 
#           legend=TRUE, ellipse = TRUE) 
loadUKBBN_UKBBNnet = data.frame(plotLoadings(splsdaUKBBN_UKBBNnet, comp = 1, contrib = "max", size.legend = 0.7, title = "CMUKBBN"))
loadUKBBN_UKBBNnet$module = as.factor(rownames(loadUKBBN_UKBBNnet))

loadUKBBN_UKBBNnet$cohort = "UKBBN"

models_network_UKBBN[[5]][[2]] = splsdaUKBBN_UKBBNnet


merge_loadings_UKBBN = rbind(loadBDR_UKBBNnet, loadUKBBN_UKBBNnet)

ggplot(merge_loadings_UKBBN, aes(x = X.importance, y = module, fill = X.GroupContrib)) +
  ggtitle("CM_BDRtoUKBBN") +
  geom_bar(stat = "identity") +
  xlab("Loadings") + ylab("Module") +
  facet_wrap(~cohort) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))


##### Module membership #####

setwd(in_dir)

print("Calculating Module Memberships...")
trait_ = as.data.frame(WGCNA_traitsUKBBNHC$sub1)
names(trait_) = "sub1"

MEs0 = net$MEs
MEs = orderMEs(MEs0)
modNames = substring(names(MEs), 3)
nSamples = ncol(betas)
geneModuleMembership = as.data.frame(cor(t(betas), MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(t(betas), as.numeric(as.factor(trait_[,trait])), use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(trait_), sep="")
names(GSPvalue) = paste("p.GS.", names(trait_), sep="")
modules = modNames

out_pref = "WGCNA/UKBBN/modulemembership/HC/sub1/UKBBN_membership"
dir.create("WGCNA/UKBBN/modulemembership")
dir.create("WGCNA/UKBBN/modulemembership/HC")
dir.create("WGCNA/UKBBN/modulemembership/HC/sub1")

for (i in 1:length(modules)) {
  module = modules[i]
  column = match(module, modNames)
  moduleGenes = net$colors==module
  
  print(paste("Generating Plot for",module,"module..."))
  sizeGrWindow(7, 7)
  par(mfrow = c(1,1))
  pdf(paste0(out_pref,".wgcna.ModuleMembership.",module,".pdf"))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for sub1",
                     main = paste("Module membership vs. gene significance\n"),
                     abline = T,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
  dev.off()
  
  print(paste("Saving",module,"modules nodes..."))
  write.csv(names(net$colors[moduleGenes]),file = paste0(out_pref,".wgcna.genes.",module,".csv"))
}


##### comparing BDR/UKBBN modules #####
load("./WGCNA/UKBBN/UKBBN_subtyping.wgcna.network.rdat")
UKBBN_net = net
load("./WGCNA/BDR/BDR_subtyping.wgcna.network.rdat")
BDR_net = net

colors_BDR = unique(BDR_net$colors)
colors_UKBBN = unique(UKBBN_net$colors)

module_BDR = names(BDR_net$colors[BDR_net$colors == colors_BDR[1]])
module_UKBBN = names(UKBBN_net$colors[UKBBN_net$colors == colors_UKBBN[1]])

for(tmp_module_BDR in colors_BDR){
  for(tmp_module_UKBBN in colors_UKBBN){
    module_BDR = names(BDR_net$colors[BDR_net$colors == tmp_module_BDR])
    module_UKBBN = names(UKBBN_net$colors[UKBBN_net$colors == tmp_module_UKBBN])
    
    print(paste0("BDR module ", tmp_module_BDR, " / UKBBN module ", tmp_module_UKBBN, " / ",
                 length(intersect(module_BDR, module_UKBBN))))
  }
}

table(BDR_net$colors)
table(UKBBN_net$colors)


print(paste0("BDR module ", colors_BDR, " / UKBBN module ", colors_UKBBN, " / ",
             length(intersect(module_BDR, module_UKBBN))))




##### all cohorts #####

WGCNA_toanovaBDR = cbind(WGCNA_phenoBDR$HC_clusters, WGCNA_eigengenesBDR)
WGCNA_toanovaBDR = WGCNA_toanovaBDR[WGCNA_toanovaBDR[,1] != "Excluded",]
colnames(WGCNA_toanovaBDR)[1] = "clusters"
WGCNA_toanovaBDR$clusters = paste0("BDR", WGCNA_toanovaBDR$clusters)

WGCNA_toanovaUKBBN = cbind(WGCNA_phenoUKBBN$HC_clusters, WGCNA_eigengenesUKBBN)
WGCNA_toanovaUKBBN = WGCNA_toanovaUKBBN[WGCNA_toanovaUKBBN[,1] != "Excluded",]
colnames(WGCNA_toanovaUKBBN)[1] = "clusters"
WGCNA_toanovaUKBBN$clusters = paste0("UKBBN", WGCNA_toanovaUKBBN$clusters)

WGCNA_toanovaPITT = cbind(WGCNA_phenoPITT$HC_clusters, WGCNA_eigengenesPITT)
WGCNA_toanovaPITT = WGCNA_toanovaPITT[WGCNA_toanovaPITT[,1] != "Excluded",]
colnames(WGCNA_toanovaPITT)[1] = "clusters"
WGCNA_toanovaPITT$clusters = paste0("PITT", WGCNA_toanovaPITT$clusters)

WGCNA_toanova_all = rbind(WGCNA_toanovaBDR, WGCNA_toanovaUKBBN,WGCNA_toanovaPITT)

tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaUKBBN)))){
  tmptukey_table <- aov(WGCNA_toanova_all[,i] ~ clusters, data = WGCNA_toanova_all) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*2))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaUKBBN)[c(2:17)], each = 36)

# labels = paste0(rownames(stats), "  MP=", signif(stats$Zsummary.pres, digits = 4))
# names(labels) = paste0("ME",rownames(stats))

WGCNA_toanova_allmelt = melt(WGCNA_toanova_all, id.vars = "clusters")

pdf(file = "WGCNA/HC_UKBBNallcohorts_anovatukey.pdf", width = 18, height = 10)

print(ggplot(WGCNA_toanova_allmelt, aes(x = clusters, y = value, color = clusters)) +
        geom_boxplot() +
        geom_jitter(size = .5) +
        facet_wrap(~ variable, nrow = 3, scales = "free")+ 
        # labeller = labeller(variable = labels)) +
        stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                           method = "anova", label.y = max(WGCNA_toanova_allmelt$value*1, na.rm = TRUE)) +
        stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
        theme_bw() +
        theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
              axis.text.x=element_text(size=7)))
dev.off()



##### random permutations #####

pdf(file = "WGCNA/UKBBN/HC_random_anovatukey.pdf", width = 18, height = 10)

for(i in seq(1:100)){
  
  BDR_random = sample(WGCNA_phenoBDR$HC_clusters)
  
  WGCNA_toanovaBDR = cbind(BDR_random, WGCNA_eigengenesBDR)
  WGCNA_toanovaBDR = WGCNA_toanovaBDR[WGCNA_toanovaBDR[,1] != "Excluded",]
  colnames(WGCNA_toanovaBDR)[1] = "clusters"
  WGCNA_toanovaBDR$clusters = paste0("BDR", WGCNA_toanovaBDR$clusters)
  
  UKBBN_random = sample(WGCNA_phenoUKBBN$HC_clusters)
  
  WGCNA_toanovaUKBBN = cbind(UKBBN_random, WGCNA_eigengenesUKBBN)
  WGCNA_toanovaUKBBN = WGCNA_toanovaUKBBN[WGCNA_toanovaUKBBN[,1] != "Excluded",]
  colnames(WGCNA_toanovaUKBBN)[1] = "clusters"
  WGCNA_toanovaUKBBN$clusters = paste0("UKBBN", WGCNA_toanovaUKBBN$clusters)
  
  PITT_random = sample(WGCNA_phenoPITT$HC_clusters)
  
  WGCNA_toanovaPITT = cbind(PITT_random, WGCNA_eigengenesPITT)
  WGCNA_toanovaPITT = WGCNA_toanovaPITT[WGCNA_toanovaPITT[,1] != "Excluded",]
  colnames(WGCNA_toanovaPITT)[1] = "clusters"
  WGCNA_toanovaPITT$clusters = paste0("PITT", WGCNA_toanovaPITT$clusters)
  
  WGCNA_toanova_all = rbind(WGCNA_toanovaBDR, WGCNA_toanovaUKBBN,WGCNA_toanovaPITT)
  
  
  tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
  for(i in c(2:(ncol(WGCNA_toanovaUKBBN)))){
    tmptukey_table <- aov(WGCNA_toanova_all[,i] ~ clusters, data = WGCNA_toanova_all) %>%
      tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanova_all[,-1])*1.1, 
                                            max(WGCNA_toanova_all[,-1])*1.2, 
                                            max(WGCNA_toanova_all[,-1])*1.3,
                                            max(WGCNA_toanova_all[,-1])*1.4, 
                                            max(WGCNA_toanova_all[,-1])*1.5, 
                                            max(WGCNA_toanova_all[,-1])*1.6,
                                            max(WGCNA_toanova_all[,-1])*1.7, 
                                            max(WGCNA_toanova_all[,-1])*1.8, 
                                            max(WGCNA_toanova_all[,-1])*1.9,
                                            max(WGCNA_toanova_all[,-1])*1.1, 
                                            max(WGCNA_toanova_all[,-1])*1.2, 
                                            max(WGCNA_toanova_all[,-1])*1.3,
                                            max(WGCNA_toanova_all[,-1])*1.4, 
                                            max(WGCNA_toanova_all[,-1])*1.5, 
                                            max(WGCNA_toanova_all[,-1])*1.6,
                                            max(WGCNA_toanova_all[,-1])*1.7, 
                                            max(WGCNA_toanova_all[,-1])*1.8, 
                                            max(WGCNA_toanova_all[,-1])*1.9,
                                            max(WGCNA_toanova_all[,-1])*1.1, 
                                            max(WGCNA_toanova_all[,-1])*1.2, 
                                            max(WGCNA_toanova_all[,-1])*1.3,
                                            max(WGCNA_toanova_all[,-1])*1.4, 
                                            max(WGCNA_toanova_all[,-1])*1.5, 
                                            max(WGCNA_toanova_all[,-1])*1.6,
                                            max(WGCNA_toanova_all[,-1])*1.7, 
                                            max(WGCNA_toanova_all[,-1])*1.9,
                                            max(WGCNA_toanova_all[,-1])*1.1, 
                                            max(WGCNA_toanova_all[,-1])*1.2, 
                                            max(WGCNA_toanova_all[,-1])*1.3,
                                            max(WGCNA_toanova_all[,-1])*1.4, 
                                            max(WGCNA_toanova_all[,-1])*1.5, 
                                            max(WGCNA_toanova_all[,-1])*1.6,
                                            max(WGCNA_toanova_all[,-1])*1.7, 
                                            max(WGCNA_toanova_all[,-1])*1.8, 
                                            max(WGCNA_toanova_all[,-1])*1.9,
                                            max(WGCNA_toanova_all[,-1])*2))
    
    tukey_table = rbind(tukey_table, tmptukey_table)
  }
  tukey_table$variable = rep(colnames(WGCNA_toanovaUKBBN)[c(2:17)], each = 36)
  
  # labels = paste0(rownames(stats), "  MP=", signif(stats$Zsummary.pres, digits = 4))
  # names(labels) = paste0("ME",rownames(stats))
  
  WGCNA_toanova_allmelt = melt(WGCNA_toanova_all, id.vars = "clusters")
  
  
  print(ggplot(WGCNA_toanova_allmelt, aes(x = clusters, y = value, color = clusters)) +
          geom_boxplot() +
          geom_jitter(size = .5) +
          facet_wrap(~ variable, nrow = 3, scales = "free")+ 
          # labeller = labeller(variable = labels)) +
          stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                             method = "anova", label.y = max(WGCNA_toanova_allmelt$value*1, na.rm = TRUE)) +
          stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
          theme_bw() +
          theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
                axis.text.x=element_text(size=7)))
  
}
dev.off()

##### prediction of samples based on eigen values #####

# From BDR
splsdaBDR_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesBDR, WGCNA_phenoBDR$HC_clusters, ncomp = 6)

BDRtoUKBBN_pred = predict(splsdaBDR_UKBBNnet, newdata = WGCNA_eigengenesUKBBN) 

# colors = WGCNA_phenoBDR$CM_clusters
# colors[colors == "Cluster 1"] = "red"
# colors[colors == "Cluster 2"] = "blue"
# colors[colors == "Cluster 3"] = "green"
# colors[colors == "Excluded"] = "grey"

pdf(paste0("./WGCNA/BDR/HC2_BDRtoUKBBN_projection_indiv.pdf"), width = 15, height = 15)
plotIndiv(splsdaBDR_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(BDRtoUKBBN_pred$variates[, 1], BDRtoUKBBN_pred$variates[, 2], 
       pch = 19, cex = 1.2, col = WGCNA_phenoUKBBN$HC_clusters)
dev.off()


##### projection from UKBBN #####
#### eigen ####
splsdaUKBBN_UKBBNnet = mixOmics::splsda(WGCNA_eigengenesUKBBN, WGCNA_phenoUKBBN$HC_clusters, ncomp = 6)

UKBBNtoBDR_pred = predict(splsdaUKBBN_UKBBNnet, newdata = WGCNA_eigengenesBDR) 
UKBBNtoPITT_pred = predict(splsdaUKBBN_UKBBNnet, newdata = WGCNA_eigengenesPITT) 

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoBDR_projection_indiv.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1.2, col = WGCNA_phenoBDR$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoPITT_projection_indiv.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1.2, col = WGCNA_phenoPITT$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoboth_projection_indiv.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1, col = WGCNA_phenoPITT$HC_clusters)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1, col = WGCNA_phenoBDR$HC_clusters)

dev.off()

#### cpg ####
splsdaUKBBN_UKBBNnet = mixOmics::splsda(t(data_matrix_AD_UKBBN), WGCNA_phenoUKBBN$HC_clusters, ncomp = 6)

UKBBNtoBDR_pred = predict(splsdaUKBBN_UKBBNnet, newdata = t(data_matrix_AD_BDR))
UKBBNtoPITT_pred = predict(splsdaUKBBN_UKBBNnet, newdata = t(data_matrix_AD_PITT))

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoBDR_projectionCPG_indiv.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1.2, col = WGCNA_phenoBDR$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoPITT_projectionCPG_indiv.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1.2, col = WGCNA_phenoPITT$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoboth_projectionCPG_indiv.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1, col = WGCNA_phenoPITT$HC_clusters)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1, col = WGCNA_phenoBDR$HC_clusters)

dev.off()


## CM

# From UKBBN

UKBBN_CM = WGCNA_phenoUKBBN[WGCNA_phenoUKBBN$CM_clusters != "Excluded",]
BDR_CM = WGCNA_phenoBDR[WGCNA_phenoBDR$CM_clusters != "Excluded",]

UKBBNeigen = WGCNA_eigengenesUKBBN[rownames(WGCNA_eigengenesUKBBN) %in% UKBBN_CM$Row.names,]
BDReigen = WGCNA_eigengenesBDR[rownames(WGCNA_eigengenesBDR) %in% BDR_CM$Basename,]

splsdaUKBBN_UKBBNnet = mixOmics::splsda(UKBBNeigen, UKBBN_CM$CM_clusters, ncomp = 6)

UKBBNtoBDR_pred = predict(splsdaUKBBN_UKBBNnet, newdata = BDReigen) 

colors = BDR_CM$CM_clusters
colors[colors == "Cluster 1"] = "1"
colors[colors == "Cluster 2"] = "2"
colors[colors == "Cluster 3"] = "3"

pdf(paste0("./WGCNA/UKBBN/CM_UKBBNtoBDR_projection_indiv.pdf"), width = 7, height = 7)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",
          style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1.2, col =colors)
dev.off()


##### 450k matching #####

cpg_list_450k = read.table("450k_cpg_list.txt")
colnames(cpg_list_450k) = "450k"

cpg_list_EPIC = rownames(data_matrix_AD_UKBBN)

inter_450k_EPIC = intersect(cpg_list_EPIC, cpg_list_450k$`450k`)

# WGCNA_labels
# black        blue       brown       green greenyellow        grey     magenta        pink 
# 3145        9248        8980        3280         519      241151        1380        2847 
# purple         red         tan   turquoise      yellow 
# 1368        3219         474       32099        7318 

WGCNA_labels_450k = WGCNA_labels[names(WGCNA_labels) %in% inter_450k_EPIC]
# WGCNA_labels_450k
# black        blue       brown       green greenyellow        grey     magenta        pink 
# 1175        3870        4147        1142         387      116155         654        1434 
# purple         red         tan   turquoise      yellow 
# 779        1335         276       11468        2704 

table(WGCNA_labels)
table(WGCNA_labels_450k)

WGCNA_PITT = pheno_PITT[pheno_PITT$HC_clusters != "Control",]

matrix_450k_BDR = data_matrix_AD_BDR[rownames(data_matrix_AD_BDR) %in% inter_450k_EPIC,]
matrix_450k_UKBBN = data_matrix_AD_UKBBN[rownames(data_matrix_AD_UKBBN) %in% inter_450k_EPIC,]
matrix_450k_PITT = data_matrix_AD_PITT[rownames(data_matrix_AD_PITT) %in% inter_450k_EPIC,]

BDR_eigen_450k = moduleEigengenes(t(matrix_450k_BDR), WGCNA_labels_450k)$eigengenes
UKBBN_eigen_450k = moduleEigengenes(t(matrix_450k_UKBBN), WGCNA_labels_450k)$eigengenes
PITT_eigen_450k = moduleEigengenes(t(matrix_450k_PITT), WGCNA_labels_450k)$eigengenes

##### 450k match #####
#### UKBBN ####

WGCNA_toanovaEPIC = cbind(WGCNA_phenoUKBBN$HC_clusters, WGCNA_eigengenesUKBBN)
colnames(WGCNA_toanovaEPIC)[1] = "clusters"
WGCNA_toanovaEPIC$clusters = paste0("EPIC", WGCNA_toanovaEPIC$clusters)

WGCNA_toanova450k = cbind(WGCNA_phenoUKBBN$HC_clusters, UKBBN_eigen_450k)
colnames(WGCNA_toanova450k)[1] = "clusters"
WGCNA_toanova450k$clusters = paste0("450k", WGCNA_toanova450k$clusters)

WGCNA_toanova_all = rbind(WGCNA_toanovaEPIC, WGCNA_toanova450k)

tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaEPIC)))){
  tmptukey_table <- aov(WGCNA_toanova_all[,i] ~ clusters, data = WGCNA_toanova_all) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*2, 
                                          max(WGCNA_toanova_all[,-1])*2.1, 
                                          max(WGCNA_toanova_all[,-1])*2.2,
                                          max(WGCNA_toanova_all[,-1])*2.3, 
                                          max(WGCNA_toanova_all[,-1])*2.4, 
                                          max(WGCNA_toanova_all[,-1])*2.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaEPIC)[c(2:17)], each = 15)


# labels = paste0(rownames(stats), "  MP=", signif(stats$Zsummary.pres, digits = 4))
# names(labels) = paste0("ME",rownames(stats))

WGCNA_toanova_allmelt = melt(WGCNA_toanova_all, id.vars = "clusters")
WGCNA_toanova_allmelt$clusters = factor(WGCNA_toanova_allmelt$clusters,
                                        levels = c("EPIC1","450k1","EPIC2","450k2","EPIC3","450k3"))

pdf(file = "WGCNA/HC_UKBBNUKBBNnet_EPIC450k_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanova_allmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free")+ 
  # labeller = labeller(variable = labels)) +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanova_allmelt$value*1, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
        axis.text.x=element_text(size=7))
dev.off()

#### BDR ####

WGCNA_toanovaEPIC = cbind(WGCNA_phenoBDR$HC_clusters, WGCNA_eigengenesBDR)
colnames(WGCNA_toanovaEPIC)[1] = "clusters"
WGCNA_toanovaEPIC$clusters = paste0("EPIC", WGCNA_toanovaEPIC$clusters)

WGCNA_toanova450k = cbind(WGCNA_phenoBDR$HC_clusters, BDR_eigen_450k)
colnames(WGCNA_toanova450k)[1] = "clusters"
WGCNA_toanova450k$clusters = paste0("450k", WGCNA_toanova450k$clusters)

WGCNA_toanova_all = rbind(WGCNA_toanovaEPIC, WGCNA_toanova450k)

tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaEPIC)))){
  tmptukey_table <- aov(WGCNA_toanova_all[,i] ~ clusters, data = WGCNA_toanova_all) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*2, 
                                          max(WGCNA_toanova_all[,-1])*2.1, 
                                          max(WGCNA_toanova_all[,-1])*2.2,
                                          max(WGCNA_toanova_all[,-1])*2.3, 
                                          max(WGCNA_toanova_all[,-1])*2.4, 
                                          max(WGCNA_toanova_all[,-1])*2.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaEPIC)[c(2:17)], each = 15)


# labels = paste0(rownames(stats), "  MP=", signif(stats$Zsummary.pres, digits = 4))
# names(labels) = paste0("ME",rownames(stats))

WGCNA_toanova_allmelt = melt(WGCNA_toanova_all, id.vars = "clusters")
WGCNA_toanova_allmelt$clusters = factor(WGCNA_toanova_allmelt$clusters,
                                        levels = c("EPIC1","450k1","EPIC2","450k2","EPIC3","450k3"))

pdf(file = "WGCNA/HC_BDRUKBBN_EPIC450k_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanova_allmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free")+ 
  # labeller = labeller(variable = labels)) +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanova_allmelt$value*1, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
        axis.text.x=element_text(size=7))
dev.off()


#### PITT ####

WGCNA_toanovaEPIC = cbind(WGCNA_phenoPITT$HC_clusters, WGCNA_eigengenesPITT)
colnames(WGCNA_toanovaEPIC)[1] = "clusters"
WGCNA_toanovaEPIC$clusters = paste0("EPIC", WGCNA_toanovaEPIC$clusters)

WGCNA_toanova450k = cbind(WGCNA_phenoPITT$HC_clusters, PITT_eigen_450k)
colnames(WGCNA_toanova450k)[1] = "clusters"
WGCNA_toanova450k$clusters = paste0("450k", WGCNA_toanova450k$clusters)

WGCNA_toanova_all = rbind(WGCNA_toanovaEPIC, WGCNA_toanova450k)

tukey_table = data.frame(matrix(nrow = 0,ncol = 9))
for(i in c(2:(ncol(WGCNA_toanovaEPIC)))){
  tmptukey_table <- aov(WGCNA_toanova_all[,i] ~ clusters, data = WGCNA_toanova_all) %>%
    tukey_hsd() %>% mutate(y.position = c(max(WGCNA_toanova_all[,-1])*1.1, 
                                          max(WGCNA_toanova_all[,-1])*1.2, 
                                          max(WGCNA_toanova_all[,-1])*1.3,
                                          max(WGCNA_toanova_all[,-1])*1.4, 
                                          max(WGCNA_toanova_all[,-1])*1.5, 
                                          max(WGCNA_toanova_all[,-1])*1.6,
                                          max(WGCNA_toanova_all[,-1])*1.7, 
                                          max(WGCNA_toanova_all[,-1])*1.8, 
                                          max(WGCNA_toanova_all[,-1])*1.9,
                                          max(WGCNA_toanova_all[,-1])*2, 
                                          max(WGCNA_toanova_all[,-1])*2.1, 
                                          max(WGCNA_toanova_all[,-1])*2.2,
                                          max(WGCNA_toanova_all[,-1])*2.3, 
                                          max(WGCNA_toanova_all[,-1])*2.4, 
                                          max(WGCNA_toanova_all[,-1])*2.5))
  
  tukey_table = rbind(tukey_table, tmptukey_table)
}
tukey_table$variable = rep(colnames(WGCNA_toanovaEPIC)[c(2:17)], each = 15)


# labels = paste0(rownames(stats), "  MP=", signif(stats$Zsummary.pres, digits = 4))
# names(labels) = paste0("ME",rownames(stats))

WGCNA_toanova_allmelt = melt(WGCNA_toanova_all, id.vars = "clusters")
WGCNA_toanova_allmelt$clusters = factor(WGCNA_toanova_allmelt$clusters,
                                        levels = c("EPIC1","450k1","EPIC2","450k2","EPIC3","450k3"))

pdf(file = "WGCNA/HC_PITTUKBBN_EPIC450k_anovatukey.pdf", width = 18, height = 10)
ggplot(WGCNA_toanova_allmelt, aes(x = clusters, y = value, color = clusters)) +
  geom_boxplot() +
  geom_jitter(size = .5) +
  facet_wrap(~ variable, nrow = 3, scales = "free")+ 
  # labeller = labeller(variable = labels)) +
  stat_compare_means(size = 3,aes(label = paste0(after_stat(method), ", p-value = ", after_stat(p.format))),
                     method = "anova", label.y = max(WGCNA_toanova_allmelt$value*1, na.rm = TRUE)) +
  stat_pvalue_manual(tukey_table, size = 3, label = "p.adj.signif", tip.length = 0.03, bracket.shorten = 0.02, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
        axis.text.x=element_text(size=7))
dev.off()


##### projection from UKBBN in 450k array #####
#### eigen ####
splsdaUKBBN_UKBBNnet = mixOmics::splsda(UKBBN_eigen_450k, WGCNA_phenoUKBBN$HC_clusters, ncomp = 6)

UKBBNtoBDR_pred = predict(splsdaUKBBN_UKBBNnet, newdata = BDR_eigen_450k) 
UKBBNtoPITT_pred = predict(splsdaUKBBN_UKBBNnet, newdata = PITT_eigen_450k) 

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoBDR_projection_indiv_450k.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1.2, col = WGCNA_phenoBDR$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoPITT_projection_indiv_450k.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1.2, col = WGCNA_phenoPITT$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoboth_projection_indiv_450k.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1, col = WGCNA_phenoPITT$HC_clusters)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1, col = WGCNA_phenoBDR$HC_clusters)

dev.off()

#### cpg ####
splsdaUKBBN_UKBBNnet = mixOmics::splsda(t(matrix_450k_UKBBN), WGCNA_phenoUKBBN$HC_clusters, ncomp = 6)

UKBBNtoBDR_pred = predict(splsdaUKBBN_UKBBNnet, newdata = t(matrix_450k_BDR))
UKBBNtoPITT_pred = predict(splsdaUKBBN_UKBBNnet, newdata = t(matrix_450k_PITT))

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoBDR_projectionCPG_indiv_450k.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1.2, col = WGCNA_phenoBDR$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoPITT_projectionCPG_indiv_450k.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1.2, col = WGCNA_phenoPITT$HC_clusters)
dev.off()

pdf(paste0("WGCNA/UKBBN/HC_UKBBNtoboth_projectionCPG_indiv_450k.pdf"), width = 10, height = 10)
plotIndiv(splsdaUKBBN_UKBBNnet, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE, legend = TRUE)
points(UKBBNtoPITT_pred$variates[, 1], UKBBNtoPITT_pred$variates[, 2], 
       pch = 3, cex = 1, col = WGCNA_phenoPITT$HC_clusters)
points(UKBBNtoBDR_pred$variates[, 1], UKBBNtoBDR_pred$variates[, 2], 
       pch = 19, cex = 1, col = WGCNA_phenoBDR$HC_clusters)

dev.off()
}