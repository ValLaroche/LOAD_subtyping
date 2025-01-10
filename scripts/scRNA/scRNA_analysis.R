##### scRNA #####
setwd("D:/valentin/main/")
load("pheno_ROSMAP.Rdata")
source("scripts/utils.R")
#Load packages
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)
library(speckle)
library(GeneOverlap)
library(clusterProfiler)
library(scales)
load("new_scRNA_MIC_res.Rdata")
#### Load data ####
pheno_ROSMAP_scrna = readRDS("./scRNA/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")

meta_ROSMAP_scRNA = read.csv("scRNA/MIT_ROSMAP_Multiomics_individual_metadata.csv")
meta_ROSMAP_scRNA = meta_ROSMAP_scRNA[meta_ROSMAP_scRNA$individualID %in% pheno_ROSMAP$individualID,]
meta_ROSMAP_scRNA = unique(meta_ROSMAP_scRNA)

pheno_ROSMAP = pheno_ROSMAP[pheno_ROSMAP$individualID %in% meta_ROSMAP_scRNA$individualID,]

pheno_ROSMAP_scrna = pheno_ROSMAP_scrna[pheno_ROSMAP_scrna$brainRegion == "PFC",]
pheno_ROSMAP_scrna = pheno_ROSMAP_scrna[pheno_ROSMAP_scrna$subject %in% meta_ROSMAP_scRNA$subject,]  
pheno_ROSMAP_scrna = pheno_ROSMAP_scrna[!pheno_ROSMAP_scrna$seurat_clusters %in% c(9,13,14,15),]

meta_ROSMAP_scRNA$AD3types = pheno_ROSMAP_scrna[match(meta_ROSMAP_scRNA$subject, pheno_ROSMAP_scrna$subject),]$ADdiag3types
pheno_ROSMAP$AD3type = meta_ROSMAP_scRNA[match(pheno_ROSMAP$individualID, meta_ROSMAP_scRNA$individualID),]$AD3type

pheno_ROSMAP$individualID
meta_ROSMAP_scRNA$individualID

scRNA_MIC = readRDS("scRNA/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds")
scRNA_MIC = scRNA_MIC[,colnames(scRNA_MIC) %in% rownames(pheno_ROSMAP_scrna)]

scRNA_counts = as.matrix(scRNA_MIC)

scRNA_counts = scRNA_counts[rowSums(scRNA_counts) > 0,]
gc()
# Perform QC
qc = perCellQCMetrics(scRNA_MIC)

# Detect outlisers
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)

scRNA_MIC <- scRNA_MIC[, !ol]
scRNA_counts <- scRNA_counts[, !ol]

# Keep expressed genes
scRNA_counts <- scRNA_counts[rowSums(scRNA_counts) >= 10, ]
gc()

pheno_ROSMAP_scrna = pheno_ROSMAP_scrna[rownames(pheno_ROSMAP_scrna) %in% colnames(scRNA_counts),]
pheno_ROSMAP = pheno_ROSMAP[order(pheno_ROSMAP$individualID),]
meta_ROSMAP_scRNA = meta_ROSMAP_scRNA[order(meta_ROSMAP_scRNA$individualID),]
pheno_ROSMAP = cbind(pheno_ROSMAP, meta_ROSMAP_scRNA)
pheno_ROSMAP_scrna$subtype = pheno_ROSMAP[match(pheno_ROSMAP_scrna$subject, pheno_ROSMAP$subject),]$subtype
pheno_ROSMAP_scrna$celltype = paste0("MIC_", pheno_ROSMAP_scrna$seurat_clusters)
unique(pheno_ROSMAP_scrna$celltype)

#### Filter for group ####

pheno_ROSMAP_scrna_RvC = pheno_ROSMAP_scrna[pheno_ROSMAP_scrna$subtype %in% c("red","Control"),]
scRNA_counts_RvC = scRNA_counts[,colnames(scRNA_counts) %in% rownames(pheno_ROSMAP_scrna_RvC)]
pheno_ROSMAP_scrna_BvC = pheno_ROSMAP_scrna[pheno_ROSMAP_scrna$subtype %in% c("blue","Control"),]
scRNA_counts_BvC = scRNA_counts[,colnames(scRNA_counts) %in% rownames(pheno_ROSMAP_scrna_BvC)]
pheno_ROSMAP_scrna_RvB = pheno_ROSMAP_scrna[pheno_ROSMAP_scrna$subtype %in% c("red","blue"),]
scRNA_counts_RvB = scRNA_counts[,colnames(scRNA_counts) %in% rownames(pheno_ROSMAP_scrna_RvB)]

#### Analysis ####
scRNA_sce = SingleCellExperiment(scRNA_counts)
names(assays(scRNA_sce)) = "counts"
scRNA_sce <- computeLibraryFactors(scRNA_sce)
scRNA_sce <- logNormCounts(scRNA_sce)

scRNA_sce_RvC = SingleCellExperiment(scRNA_counts_RvC)
names(assays(scRNA_sce_RvC)) = "counts"
scRNA_sce_RvC <- computeLibraryFactors(scRNA_sce_RvC)
scRNA_sce_RvC <- logNormCounts(scRNA_sce_RvC)

scRNA_sce_BvC = SingleCellExperiment(scRNA_counts_BvC)
names(assays(scRNA_sce_BvC)) = "counts"
scRNA_sce_BvC <- computeLibraryFactors(scRNA_sce_BvC)
scRNA_sce_BvC <- logNormCounts(scRNA_sce_BvC)

scRNA_sce_RvB = SingleCellExperiment(scRNA_counts_RvB)
names(assays(scRNA_sce_RvB)) = "counts"
scRNA_sce_RvB <- computeLibraryFactors(scRNA_sce_RvB)
scRNA_sce_RvB <- logNormCounts(scRNA_sce_RvB)

gc()

#### UMAP ####

scRNA_sce@colData$cluster_id = pheno_ROSMAP_scrna$celltype
scRNA_sce@colData$subject = pheno_ROSMAP_scrna$subject
scRNA_sce@colData$group_id = pheno_ROSMAP[match(scRNA_sce@colData$subject, pheno_ROSMAP$subject),]$diag
scRNA_sce@colData$gender = pheno_ROSMAP_scrna$msex
scRNA_sce@colData$pmi = pheno_ROSMAP_scrna$pmi
scRNA_sce@colData$age_death = pheno_ROSMAP_scrna$age_death
scRNA_sce@colData$sample_id = paste(scRNA_sce$group_id, scRNA_sce$subject, sep = "_")

scRNA_sce_RvC@colData$cluster_id = pheno_ROSMAP_scrna_RvC$celltype
scRNA_sce_RvC@colData$group_id = pheno_ROSMAP_scrna_RvC$subtype
scRNA_sce_RvC@colData$sample_id = paste(pheno_ROSMAP_scrna_RvC$subtype, pheno_ROSMAP_scrna_RvC$subject, sep = "_")
scRNA_sce_RvC@colData$gender = pheno_ROSMAP_scrna_RvC$msex
scRNA_sce_RvC@colData$pmi = pheno_ROSMAP_scrna_RvC$pmi
scRNA_sce_RvC@colData$age_death = pheno_ROSMAP_scrna_RvC$age_death

scRNA_sce_BvC@colData$cluster_id = pheno_ROSMAP_scrna_BvC$celltype
scRNA_sce_BvC@colData$group_id = pheno_ROSMAP_scrna_BvC$subtype
scRNA_sce_BvC@colData$sample_id = paste(pheno_ROSMAP_scrna_BvC$subtype, pheno_ROSMAP_scrna_BvC$subject, sep = "_")
scRNA_sce_BvC@colData$gender = pheno_ROSMAP_scrna_BvC$msex
scRNA_sce_BvC@colData$pmi = pheno_ROSMAP_scrna_BvC$pmi
scRNA_sce_BvC@colData$age_death = pheno_ROSMAP_scrna_BvC$age_death

scRNA_sce_RvB@colData$cluster_id = pheno_ROSMAP_scrna_RvB$celltype
scRNA_sce_RvB@colData$group_id = pheno_ROSMAP_scrna_RvB$subtype
scRNA_sce_RvB@colData$sample_id = paste(pheno_ROSMAP_scrna_RvB$subtype, pheno_ROSMAP_scrna_RvB$subject, sep = "_")
scRNA_sce_RvB@colData$gender = pheno_ROSMAP_scrna_RvB$msex
scRNA_sce_RvB@colData$pmi = pheno_ROSMAP_scrna_RvB$pmi
scRNA_sce_RvB@colData$age_death = pheno_ROSMAP_scrna_RvB$age_death



table_MIC_cells_RvC = as.data.frame(t(table(scRNA_sce_RvC$cluster_id, scRNA_sce_RvC$sample_id)))
table_MIC_cells_RvC = dcast(data = table_MIC_cells_RvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

table_MIC_cells_BvC = as.data.frame(t(table(scRNA_sce_BvC$cluster_id, scRNA_sce_BvC$sample_id)))
table_MIC_cells_BvC = dcast(data = table_MIC_cells_BvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

table_MIC_cells_ADvC = as.data.frame(t(table(scRNA_sce$cluster_id, scRNA_sce$sample_id)))
table_MIC_cells_ADvC = dcast(data = table_MIC_cells_ADvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

table_MIC_cells_ADvC = table_MIC_cells_ADvC %>% separate(Var1, c('diag', 'id'), sep = "_")
table_MIC_cells_RvC = table_MIC_cells_RvC %>% separate(Var1, c('diag', 'id'), sep = "_")
table_MIC_cells_BvC = table_MIC_cells_BvC %>% separate(Var1, c('diag', 'id'), sep = "_")

table_all_MIC_cells = rbind(table_MIC_cells_ADvC, 
                            table_MIC_cells_RvC[table_MIC_cells_RvC$diag == "red",],
                            table_MIC_cells_BvC[table_MIC_cells_BvC$diag == "blue",])

scRNA_sce_RvC = prepSCE(scRNA_sce_RvC)
t(table(scRNA_sce_RvC$cluster_id, scRNA_sce_RvC$sample_id))

# scRNA_sce_RvC <- runUMAP(scRNA_sce_RvC, pca = 20)
# 
# plot_dr <- function(sce, dr, col){
#   plotReducedDim(sce, dimred = dr, colour_by = col) +
#   guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#   theme_minimal() + theme(aspect.ratio = 1)
# }
# 
# cs_by_k <- split(colnames(scRNA_sce_RvC), scRNA_sce_RvC$cluster_id)
# cs100 <- unlist(sapply(cs_by_k, function(u) sample(u, min(length(u), 100))))
# 
# # plot t-SNE & UMAP colored by cluster & group ID
# for (dr in c("UMAP")){ 
#   for (col in c("cluster_id", "group_id")){ 
#     plot(plot_dr(scRNA_sce_RvC[, cs100], dr, col))
#   }
# }
# ###

scRNA_sce_BvC = prepSCE(scRNA_sce_BvC)
t(table(scRNA_sce_BvC$cluster_id, scRNA_sce_BvC$sample_id))

# scRNA_sce_BvC <- runUMAP(scRNA_sce_BvC, pca = 20)
# 
# cs_by_k <- split(colnames(scRNA_sce_BvC), scRNA_sce_BvC$cluster_id)
# cs100 <- unlist(sapply(cs_by_k, function(u) sample(u, min(length(u), 100))))
# 
# # plot t-SNE & UMAP colored by cluster & group ID
# for (dr in c("UMAP")){ 
#   for (col in c("cluster_id", "group_id")){ 
#     plot(plot_dr(scRNA_sce_BvC[, cs100], dr, col))
#   }
# }
# 
# #### Pseudobulk ####
pb_RvC <- aggregateData(scRNA_sce_RvC,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb_RvC)
t(head(assay(pb_RvC)))

# pb_mds_RvC <- pbMDS(pb_RvC)
# plot(pb_mds_RvC)

#Filter by expr at the pseudobulk level
#Filter gene is by expression, filter samples is outlier detection
formula = ~group_id+gender+pmi+age_death
cd = as.data.frame(colData(pb_RvC))
cd2=cd[,c("group_id","gender","pmi","age_death")]
design=model.matrix(formula,cd2)
pb_RvC = pb_RvC[,rownames(design)]

res_RvC <- pbDS(pb_RvC, design = design, min_cells = 10, filter = "both", verbose = TRUE, coef=as.list(2:5))
tbl_RvC <- res_RvC$table[[1]]
# 
# vapply(res_RvC$table$red, function(u)
#     sum(u$p_val < 0.01), numeric(1))
# 
# 
# tbl_RvC = lapply(tbl_RvC, function(u)
#   filter(u, abs(logFC) > log2(1.2)))
tbl_RvC = lapply(tbl_RvC, function(u)
  filter(u, p_val < 0.05))
tbl_RvC
# 
# gc()
# 
pb_BvC <- aggregateData(scRNA_sce_BvC,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
pb_BvC$group_id = factor(pb_BvC$group_id, levels = c("Control","blue"))
# one sheet per subpopulation

assayNames(pb_BvC)
t(head(assay(pb_BvC)))

# pb_mds_BvC <- pbMDS(pb_BvC)
# plot(pb_mds_BvC)
# 
# metadata(pb_BvC)

#Filter gene is by expression, filter samples is outlier detection
formula = ~group_id+gender+pmi+age_death
cd = as.data.frame(colData(pb_BvC))
cd2=cd[,c("group_id","gender","pmi","age_death")]
design=model.matrix(formula,cd2)

res_BvC <- pbDS(pb_BvC, design = design, min_cells = 10, filter = "none", verbose = TRUE, coef=as.list(2:5))
tbl_BvC <- res_BvC$table[[1]]

# vapply(res_BvC$table$blue, function(u)
#   sum(u$p_val < 0.01), numeric(1))
# 
# # 
# tbl_BvC = lapply(tbl_BvC, function(u)
#   filter(u, abs(logFC) > log2(1.2)))
tbl_BvC = lapply(tbl_BvC, function(u)
  filter(u, p_val < 0.05))
tbl_BvC
# 

# RvB
scRNA_sce_RvB = prepSCE(scRNA_sce_RvB)
t(table(scRNA_sce_RvB$cluster_id, scRNA_sce_RvB$sample_id))
pb_RvB <- aggregateData(scRNA_sce_RvB,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb_RvB)
t(head(assay(pb_RvB)))

# pb_mds_RvB <- pbMDS(pb_RvB)
# plot(pb_mds_RvB)

#Filter by expr at the pseudobulk level
#Filter gene is by expression, filter samples is outlier detection
formula = ~group_id+gender+pmi+age_death
cd = as.data.frame(colData(pb_RvB))
cd2=cd[,c("group_id","gender","pmi","age_death")]
design=model.matrix(formula,cd2)
pb_RvB = pb_RvB[,rownames(design)]

res_RvB <- pbDS(pb_RvB, design = design, min_cells = 10, filter = "both", verbose = TRUE, coef=as.list(2:5))
tbl_RvB <- res_RvB$table[[1]]
# 
# vapply(res_RvB$table$red, function(u)
#     sum(u$p_val < 0.01), numeric(1))
# 
# 
# tbl_RvB = lapply(tbl_RvB, function(u)
#   filter(u, abs(logFC) > log2(1.2)))
tbl_RvB = lapply(tbl_RvB, function(u)
  filter(u, p_val < 0.05))
tbl_RvB





# Check overlap between red and blue microglia state DEA results
overlap_MIC_states = c()
for(i in assayNames(pb_BvC)){
  if(i %in% names(tbl_BvC) & i %in% names(tbl_RvC)){
    print(i)
    print(intersect(tbl_BvC[[i]]$gene, tbl_RvC[[i]]$gene))
  }
}

# Get results from microglia states study and overlapping them with our DEGs
res_DEA_MICstate = read.table("scRNA/ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt")
overlap_MICstate_RED = list()
for(i in names(tbl_RvC)){
  print(i)
  tmp_micstate = res_DEA_MICstate[res_DEA_MICstate$cluster == strsplit(i, split = "_")[[1]][2],]
  print(tbl_RvC[[i]][tbl_RvC[[i]]$gene %in% intersect(tbl_RvC[[i]]$gene, tmp_micstate$gene),])
  overlap_MICstate_RED[[i]] = tbl_RvC[[i]][tbl_RvC[[i]]$gene %in% intersect(tbl_RvC[[i]]$gene, tmp_micstate$gene),]
}
overlap_MICstate_BLUE = list()
for(i in names(tbl_BvC)){
  print(i)
  tmp_micstate = res_DEA_MICstate[res_DEA_MICstate$cluster == strsplit(i, split = "_")[[1]][2],]
  print(tbl_BvC[[i]][tbl_BvC[[i]]$gene %in% intersect(tbl_BvC[[i]]$gene, tmp_micstate$gene),])
  overlap_MICstate_BLUE[[i]] = tbl_BvC[[i]][tbl_BvC[[i]]$gene %in% intersect(tbl_BvC[[i]]$gene, tmp_micstate$gene),] 
}

# Testing for overrepresentation
all_tables_RvC = list()
all_fishers_RvC = list()
for(i in names(tbl_RvC)){
  tmp_micstate = res_DEA_MICstate[res_DEA_MICstate$cluster == strsplit(i, split = "_")[[1]][2],]
  table_red_MICstate = matrix(c(nrow(overlap_MICstate_RED[[i]]), #overlap study + RvC
                                nrow(tmp_micstate), # total RvC
                                nrow(tbl_RvC[[i]]), # total study
                                dim(scRNA_sce_RvC@assays)[1]), #total genes in RNA array after QC
                           nrow = 2)
  all_tables_RvC[[i]] = table_red_MICstate
  fisher_red_MIC = fisher.test(table_red_MICstate, alternative = "greater")
  all_fishers_RvC[[i]] = fisher_red_MIC
}

all_tables_BvC = list()
all_fishers_BvC = list()
for(i in names(tbl_BvC)){
  tmp_micstate = res_DEA_MICstate[res_DEA_MICstate$cluster == strsplit(i, split = "_")[[1]][2],]
  table_blue_MICstate = matrix(c(nrow(overlap_MICstate_BLUE[[i]]), #overlap study + BvC
                                 nrow(tbl_BvC[[i]]), # total BvC
                                 nrow(tmp_micstate), # total study
                                 dim(scRNA_sce_BvC@assays)[1]), #total genes in RNA array after QC
                               nrow = 2)
  all_tables_BvC[[i]] = table_blue_MICstate
  fisher_blue_MIC = fisher.test(table_blue_MICstate, alternative = "greater")
  all_fishers_BvC[[i]] = fisher_blue_MIC
}

# Generate the complete object and the plot
table_res = data.frame(matrix(ncol = 8, nrow = 12))
colnames(table_res) = c("total_RNA_genes","MIC_state_sig", "BLUE_sig", "BLUE_overlap", "BLUE_fishers",
                        "RED_sig", "RED_overlap", "RED_fishers")
rownames(table_res) = c("MIC_0", "MIC_1","MIC_2","MIC_3","MIC_4","MIC_5",
                        "MIC_6","MIC_7","MIC_8","MIC_10","MIC_11","MIC_12")
for(i in rownames(table_res)){
  tmp_micstate = res_DEA_MICstate[res_DEA_MICstate$cluster == strsplit(i, split = "_")[[1]][2],]
  table_res[i,] = c(dim(scRNA_sce_BvC@assays)[1],
                    nrow(tmp_micstate),
                    if(!is.null(nrow(tbl_BvC[[i]]))){nrow(tbl_BvC[[i]])}else{0},
                    if(!is.null(nrow(overlap_MICstate_BLUE[[i]]))){nrow(overlap_MICstate_BLUE[[i]])}else{0},
                    if(!is.null(all_fishers_BvC[[i]]$p.value)){all_fishers_BvC[[i]]$p.value}else{1},
                    if(!is.null(nrow(tbl_RvC[[i]]))){nrow(tbl_RvC[[i]])}else{0},
                    if(!is.null(nrow(overlap_MICstate_RED[[i]]))){nrow(overlap_MICstate_RED[[i]])}else{0},
                    if(!is.null(all_fishers_RvC[[i]]$p.value)){all_fishers_RvC[[i]]$p.value}else{1})
                    
                    
}


table_res_melt = melt(table_res)
table_res_melt = table_res_melt[-c(1:12, 49:60,85:96),]

table_res_melt$group = str_split_fixed(as.character(table_res_melt$variable), "_",2)[,1]
table_res_melt$cell_type = rep(rownames(table_res), 5)
table_res_melt$bin_overlap = str_split_fixed(as.character(table_res_melt$variable), "_",2)[,2]
table_res_melt[table_res_melt$bin_overlap == "state_sig",]$bin_overlap = "sig"

table_res_melt$pvalue = c(rep(1, 24), table_res$BLUE_fishers, rep(1, 12), table_res$RED_fishers) 

table_res_melt$pvalue = as.numeric(format(table_res_melt$pvalue, digits = 3, scientific = TRUE))

table_res_melt$cell_type = gsub("_", " ",table_res_melt$cell_type)

table_res_melt$cell_type = factor(table_res_melt$cell_type, levels = gsub("_", " ",rownames(table_res)))

table_res_melt[table_res_melt$group == "MIC",]$group = "Microglia state"
table_res_melt[table_res_melt$group == "RED",]$group = "RED subtype"
table_res_melt[table_res_melt$group == "BLUE",]$group = "BLUE subtype"

table_res_melt = table_res_melt[!table_res_melt$cell_type == "MIC 12",]
table_res_melt$group = factor(table_res_melt$group, levels = c("Microglia state", "BLUE subtype", "RED subtype"))
table_res_melt$log_value = log(table_res_melt$value + 1)

table_res_melt = table_res_melt[!table_res_melt$variable %in% c("BLUE_sig", "RED_sig"),]

ggplot(table_res_melt, aes(x = cell_type, y = log_value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5) +
  geom_text(aes(label = if_else(pvalue < 1e-2,"*","")), position = position_dodge(width = 0.9), vjust = -0.5, size = 8) +
  geom_text(aes(y = -0.4, label =group ), position = position_dodge(width = 0.9), vjust = 0.5, angle = 90, size = 3) +
  scale_y_continuous(breaks = c(-2,0,1,2,3,4,5,6)) + ylim(-0.5,6) +
  xlab("Microglia state") + ylab("Number of significant genes") +
  # geom_hline(yintercept=1.05) + geom_hline(yintercept=1.15) +
  scale_fill_manual(values = c("black", "steelblue1", "brown2"), name = "Microglia enrichment",
                      labels = c("Microglia state", "BLUE state specific", "RED state specific")) + 
  theme_bw() +
  theme(axis.text.x = element_text(size=10),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

for(i in names(tbl_RvC)){
  print(tbl_RvC[[i]]$gene)
}

# save.image("scRNA/current_scRNA.Rdata")

gc()

#### Table gen ####
table_globalBvC = data.frame(matrix(nrow = 0, ncol = 9))
for(i in names(tbl_BvC)){
  table_globalBvC = rbind(table_globalBvC, tbl_BvC[[i]][order(tbl_BvC[[i]]$p_val),])
}
table_globalRvC = data.frame(matrix(nrow = 0, ncol = 9))
for(i in names(tbl_RvC)){
  table_globalRvC = rbind(table_globalRvC, tbl_RvC[[i]][order(tbl_RvC[[i]]$p_val),])
}
# write.csv(table_globalBvC, file = "scRNA/table_globalBvC.csv")
# write.csv(table_globalRvC, file = "scRNA/table_globalRvC.csv")

#Overlap
table_overlapBvC = data.frame(matrix(nrow = 0, ncol = 9))
for(i in names(overlap_MICstate_BLUE)){
  table_overlapBvC = rbind(table_overlapBvC, overlap_MICstate_BLUE[[i]][order(overlap_MICstate_BLUE[[i]]$p_val),])
}
table_overlapRvC = data.frame(matrix(nrow = 0, ncol = 9))
for(i in names(overlap_MICstate_RED)){
  table_overlapRvC = rbind(table_overlapRvC, overlap_MICstate_RED[[i]][order(overlap_MICstate_RED[[i]]$p_val),])
}
# write.csv(table_overlapBvC, file = "scRNA/table_overlapBvC.csv")
# write.csv(table_overlapRvC, file = "scRNA/table_overlapRvC.csv")

intersect_RB = data.frame(matrix(nrow = 0, ncol = 9))
for(i in intersect(names(tbl_BvC), names(tbl_RvC))){
  inter_genes = intersect(tbl_BvC[[i]]$gene, tbl_RvC[[i]]$gene)
  print(inter_genes)
  inter_R = tbl_RvC[[i]][tbl_RvC[[i]]$gene %in% inter_genes,]
  inter_R = inter_R[order(inter_R$p_val),]
  inter_R$subtype = "red"
  inter_B = tbl_BvC[[i]][tbl_BvC[[i]]$gene %in% inter_genes,]
  inter_B = inter_B[order(inter_B$p_val),]
  inter_B$subtype = "blue"
  intersect_RB = rbind(intersect_RB, inter_B, inter_R)
  intersect_RB = intersect_RB[order(intersect_RB$cluster_id, intersect_RB$gene),]
}

intersect_overlapRB = data.frame(matrix(nrow = 0, ncol = 9))
for(i in intersect(names(overlap_MICstate_BLUE), names(overlap_MICstate_RED))){
  inter_genes = intersect(overlap_MICstate_BLUE[[i]]$gene, overlap_MICstate_RED[[i]]$gene)
  print(inter_genes)
  inter_R = overlap_MICstate_RED[[i]][overlap_MICstate_RED[[i]]$gene %in% inter_genes,]
  inter_R = inter_R[order(inter_R$p_val),]
  inter_R$subtype = "red"
  inter_B = overlap_MICstate_BLUE[[i]][overlap_MICstate_BLUE[[i]]$gene %in% inter_genes,]
  inter_B = inter_B[order(inter_B$p_val),]
  inter_B$subtype = "blue"
  intersect_overlapRB = rbind(intersect_overlapRB, inter_B, inter_R)
  intersect_overlapRB = intersect_overlapRB[order(intersect_overlapRB$cluster_id, intersect_overlapRB$gene),]
}

scRNA_sce = prepSCE(scRNA_sce)
pb_ADvC <- aggregateData(scRNA_sce,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
pb_ADvC$group_id = factor(pb_ADvC$group_id, levels = c("Control","AD"))
# one sheet per subpopulation

assayNames(pb_ADvC)
t(head(assay(pb_ADvC)))

#Filter gene is by expression, filter samples is outlier detection
formula = ~group_id+gender+pmi+age_death
cd = as.data.frame(colData(pb_ADvC))
cd2=cd[,c("group_id","gender","pmi","age_death")]
design=model.matrix(formula,cd2)

res_ADvC <- pbDS(pb_ADvC, design = design, min_cells = 10, filter = "both", verbose = TRUE, coef=as.list(2:5))
tbl_ADvC <- res_ADvC$table[[1]]

res_propeller_ADvC = propeller(clusters = scRNA_sce@colData$cluster_id, sample = scRNA_sce@colData$sample_id, 
                               group = scRNA_sce@colData$group_id)
res_propeller_RvC = propeller(clusters = scRNA_sce_RvC@colData$cluster_id, sample = scRNA_sce_RvC@colData$sample_id, 
                               group = scRNA_sce_RvC@colData$group_id)
res_propeller_BvC = propeller(clusters = scRNA_sce_BvC@colData$cluster_id, sample = scRNA_sce_BvC@colData$sample_id, 
                               group = scRNA_sce_BvC@colData$group_id)
res_propeller_RvB = propeller(clusters = scRNA_sce_RvB@colData$cluster_id, sample = scRNA_sce_RvB@colData$sample_id, 
                              group = scRNA_sce_RvB@colData$group_id)


# Plot cell type proportions
plotCellTypeProps(clusters = scRNA_sce@colData$cluster_id, sample = scRNA_sce@colData$sample_id)
plotCellTypeProps(clusters = scRNA_sce_RvC@colData$cluster_id, sample = scRNA_sce_RvC@colData$sample_id)
plotCellTypeProps(clusters = scRNA_sce_BvC@colData$cluster_id, sample = scRNA_sce_BvC@colData$sample_id)


# Differential expression subtypes
for(mg_state in names(tbl_ADvC)){
  print(mg_state)
  if(mg_state %in% names(tbl_BvC)){
    tS1 = tbl_BvC[[mg_state]] 
  }
  if(mg_state %in% names(tbl_RvC)){
    tS2 = tbl_RvC[[mg_state]]
  }
    
    
    
  if(mg_state %in% names(tbl_RvC) & mg_state %in% names(tbl_BvC)){
    backg = min(nrow(tS1), nrow(tS2))
    
    tS1 = tS1[tS1$p_val < 0.05 & abs(tS1$logFC) > log2(1.2),]
    tS2 = tS2[tS2$p_val < 0.05 & abs(tS2$logFC) > log2(1.2),]
    to_compare = newGeneOverlap(tS1$gene, tS2$gene, backg, spec = "hg19.gene")
    print(testGeneOverlap(to_compare))
  }
    
  if(mg_state %in% names(tbl_BvC)){
    tS1 = tS1[tS1$p_val < 0.05 & abs(tS1$logFC) > log2(1.2),]
    tS1 = tS1[order(tS1$p_val),]
    tS1_glist = sort(tS1$gene, decreasing = T)
    res_MIC2_S1 = enrichGO(tS1_glist, OrgDb = "org.Hs.eg.db", 
                           ont = "ALL", keyType = "SYMBOL",
                           pAdjustMethod = "bonferroni",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1)
    
    GO_MIC2_S1 = res_MIC2_S1@result
    GO_MIC2_S1$count = as.numeric(str_split_fixed(GO_MIC2_S1$GeneRatio, "/", n =2)[,1])
    GO_MIC2_S1$total = as.numeric(str_split_fixed(GO_MIC2_S1$GeneRatio, "/", n =2)[,2])
    GO_MIC2_S1$ratio = GO_MIC2_S1$count / GO_MIC2_S1$total
    GO_MIC2_S1 = GO_MIC2_S1[GO_MIC2_S1$p.adjust < 0.05,]
    
    val = min(20, nrow(GO_MIC2_S1))
    GO_MIC2_S1_plot = GO_MIC2_S1[1:val,]
    GO_MIC2_S1_plot = GO_MIC2_S1_plot[order(GO_MIC2_S1_plot$ratio, decreasing = T),]
    GO_MIC2_S1_plot$Description = factor(GO_MIC2_S1_plot$Description, levels = GO_MIC2_S1_plot$Description)
    plot(ggplot(GO_MIC2_S1_plot, aes(x = ratio, y = rev(Description), colour = p.adjust)) +
      geom_point(aes(size = count)) +
      ggtitle(paste0(mg_state, " - LOADS1")) +
      scale_size_continuous(range = c(2,8))+
      scale_color_gradient2(low = "#DA2A2A", mid = "grey80", high = "#4483B8", midpoint = mean(GO_MIC2_S1_plot$p.adjust)) +
      scale_y_discrete(labels = label_wrap(70)) +
      theme_bw() + xlab("Gene Ratio") + ylab("Pathway") +
      theme(axis.text.y = element_text(size = 10)))
  }
  if(mg_state %in% names(tbl_RvC)){
    tS2 = tS2[tS2$p_val < 0.05 & abs(tS2$logFC) > log2(1.2),]
    tS2 = tS2[order(tS2$p_val),]
    tS2_glist = sort(tS2$gene, decreasing = T)
    res_MIC2_S2 = clusterProfiler::enrichGO(tS2_glist, OrgDb = "org.Hs.eg.db", 
                                            ont = "ALL", keyType = "SYMBOL",
                                            pAdjustMethod = "bonferroni",
                                            pvalueCutoff  = 1,
                                            qvalueCutoff  = 1)
    GO_MIC2_S2 = res_MIC2_S2@result
    GO_MIC2_S2$count = as.numeric(str_split_fixed(GO_MIC2_S2$GeneRatio, "/", n =2)[,1])
    GO_MIC2_S2$total = as.numeric(str_split_fixed(GO_MIC2_S2$GeneRatio, "/", n =2)[,2])
    GO_MIC2_S2$ratio = GO_MIC2_S2$count / GO_MIC2_S2$total
    GO_MIC2_S2 = GO_MIC2_S2[GO_MIC2_S2$p.adjust < 0.05,]
    
    val = min(20, nrow(GO_MIC2_S2))
    GO_MIC2_S2_plot = GO_MIC2_S2[1:val,]
    GO_MIC2_S2_plot = GO_MIC2_S2_plot[order(GO_MIC2_S2_plot$ratio, decreasing = T),]
    GO_MIC2_S2_plot$Description = factor(GO_MIC2_S2_plot$Description, levels = GO_MIC2_S2_plot$Description)
    plot(ggplot(GO_MIC2_S2_plot, aes(x = ratio, y = rev(Description), colour = p.adjust)) +
      ggtitle(paste0(mg_state, " - LOADS2")) +
      geom_point(aes(size = count)) +
      scale_size_continuous(range = c(2,8))+
      scale_color_gradient2(low = "#DA2A2A", mid = "grey80", high = "#4483B8", midpoint = mean(GO_MIC2_S2_plot$p.adjust)) +
      scale_y_discrete(labels = label_wrap(70)) +
      theme_bw() + xlab("Gene Ratio") + ylab("Pathway") +
      theme(axis.text.y = element_text(size = 10)))
    
    
  }
}

sig_tbl_S1 = lapply(tbl_BvC, function(u) filter(u, abs(logFC) > log2(1.2)))
sig_tbl_S1 = lapply(sig_tbl_S1, function(u) filter(u, p_val < 0.05))
sig_tbl_S2 = lapply(tbl_RvC, function(u) filter(u, abs(logFC) > log2(1.2)))
sig_tbl_S2 = lapply(sig_tbl_S2, function(u) filter(u, p_val < 0.05))

list_gene_S1 = unique(unlist(lapply(sig_tbl_S1, function(u) c(u$gene))))
list_gene_S2 = unique(unlist(lapply(sig_tbl_S2, function(u) c(u$gene))))

intersect(list_gene_S1, list_gene_S2)
