##### scRNA #####
load("pheno_ROSMAP.Rdata")
set.seed(3)
# source("scripts/utils.R")
#Load packages
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)
library(Seurat)
library(SeuratObject)
#### Load data ####

scRNA_OLIG = readRDS("./excitatory_neurons_set2.rds")
dim(scRNA_OLIG)

scRNA_OLIG = UpdateSeuratObject(scRNA_OLIG)
scRNA_OLIG = subset(scRNA_OLIG, subset = projid %in% pheno_ROSMAP$projid)

dim(scRNA_OLIG)
gc()

scRNA_OLIG_sce = as.SingleCellExperiment(scRNA_OLIG)

remove(scRNA_OLIG)

dim(scRNA_OLIG_sce)
summary(scRNA_OLIG_sce)
scRNA_OLIG_sce

gc()

scRNA_OLIG_sce <- scRNA_OLIG_sce[rowSums(counts(scRNA_OLIG_sce) > 0) > 0, ]
dim(scRNA_OLIG_sce)

# Perform QC
qc = perCellQCMetrics(scRNA_OLIG_sce)

# Detect outlisers
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)

scRNA_OLIG_sce <- scRNA_OLIG_sce[, !ol]

scRNA_OLIG_sce = scRNA_OLIG_sce[rowSums(counts(scRNA_OLIG_sce) > 0) > 0,]

# Keep expressed genes

scRNA_OLIG_sce <- scRNA_OLIG_sce[rowSums(counts(scRNA_OLIG_sce) > 1) >= 10, ]

scRNA_OLIG_sce <- computeLibraryFactors(scRNA_OLIG_sce)
scRNA_OLIG_sce <- logNormCounts(scRNA_OLIG_sce)

gc()

#### Filter for group ####
pheno_ROSMAP[is.na(pheno_ROSMAP$age_death),]$age_death = 90
pheno_ROSMAP_RvC = pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red","Control"),]
pheno_ROSMAP_BvC = pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("blue","Control"),]

coldata = scRNA_OLIG_sce@colData
coldata_RvC = coldata[coldata$projid %in% pheno_ROSMAP_RvC$projid,]
coldata_BvC = coldata[coldata$projid %in% pheno_ROSMAP_BvC$projid,]

scRNA_sce_RvC = scRNA_OLIG_sce[,colnames(scRNA_OLIG_sce) %in% rownames(coldata_RvC)]
scRNA_sce_BvC = scRNA_OLIG_sce[,colnames(scRNA_OLIG_sce) %in% rownames(coldata_BvC)]

scRNA_sce_RvC@colData$cluster_id = scRNA_sce_RvC@colData$cell_type_high_resolution
scRNA_sce_RvC@colData$group_id = pheno_ROSMAP_RvC[match(scRNA_sce_RvC@colData$projid, pheno_ROSMAP_RvC$projid),]$subtype
scRNA_sce_RvC@colData$gender = pheno_ROSMAP_RvC[match(scRNA_sce_RvC@colData$projid, pheno_ROSMAP_RvC$projid),]$msex
scRNA_sce_RvC@colData$age_death = pheno_ROSMAP_RvC[match(scRNA_sce_RvC@colData$projid, pheno_ROSMAP_RvC$projid),]$age_death
scRNA_sce_RvC@colData$pmi = pheno_ROSMAP_RvC[match(scRNA_sce_RvC@colData$projid, pheno_ROSMAP_RvC$projid),]$pmi
scRNA_sce_RvC@colData$sample_id = paste(scRNA_sce_RvC@colData$group_id, scRNA_sce_RvC@colData$projid, sep = "_")
scRNA_sce_BvC@colData$cluster_id = scRNA_sce_BvC@colData$cell_type_high_resolution
scRNA_sce_BvC@colData$group_id = pheno_ROSMAP_BvC[match(scRNA_sce_BvC@colData$projid, pheno_ROSMAP_BvC$projid),]$subtype
scRNA_sce_BvC@colData$gender = pheno_ROSMAP_BvC[match(scRNA_sce_BvC@colData$projid, pheno_ROSMAP_BvC$projid),]$msex
scRNA_sce_BvC@colData$age_death = pheno_ROSMAP_BvC[match(scRNA_sce_BvC@colData$projid, pheno_ROSMAP_BvC$projid),]$age_death
scRNA_sce_BvC@colData$pmi = pheno_ROSMAP_BvC[match(scRNA_sce_BvC@colData$projid, pheno_ROSMAP_BvC$projid),]$pmi
scRNA_sce_BvC@colData$sample_id = paste(scRNA_sce_BvC@colData$group_id, scRNA_sce_BvC@colData$projid, sep = "_")

scRNA_sce_RvC = prepSCE(scRNA_sce_RvC)
t(table(scRNA_sce_RvC$cluster_id, scRNA_sce_RvC$sample_id))

#### Pseudobulk ####

pb_RvC <- aggregateData(scRNA_sce_RvC,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb_RvC)
t(head(assay(pb_RvC)))

#Filter by expr at the pseudobulk level
#Filter gene is by expression, filter samples is outlier detection
formula = ~gender+pmi+age_death
cd = as.data.frame(colData(pb_RvC))
cd2=cd[,c("gender","pmi","age_death")]
design=model.matrix(formula,cd2)

res_RvC <- pbDS(pb_RvC, design = design, min_cells = 10, filter = "both", verbose = TRUE)
tbl_RvC <- res_RvC$table[[1]]

# tbl_RvC = lapply(tbl_RvC, function(u) 
#   filter(u, abs(logFC) > log2(1.2)))
# tbl_RvC = lapply(tbl_RvC, function(u) 
#   filter(u, p_val < 0.05))
# tbl_RvC

gc()

scRNA_sce_BvC = prepSCE(scRNA_sce_BvC)
t(table(scRNA_sce_BvC$cluster_id, scRNA_sce_BvC$sample_id))



pb_BvC <- aggregateData(scRNA_sce_BvC,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
assayNames(pb_BvC)
t(head(assay(pb_BvC)))


#Filter gene is by expression, filter samples is outlier detection
formula = ~gender+pmi+age_death
cd = as.data.frame(colData(pb_BvC))
cd2=cd[,c("gender","pmi","age_death")]
design=model.matrix(formula,cd2)

cd[is.na(cd$pmi),]$pmi = mean(design[,3])
cd2=cd[,c("gender","pmi","age_death")]
design=model.matrix(formula,cd2)

res_BvC <- pbDS(pb_BvC, design = design, min_cells = 10, filter = "both", verbose = TRUE)
tbl_BvC <- res_BvC$table[[1]]

# tbl_BvC = lapply(tbl_BvC, function(u) 
#   filter(u, abs(logFC) > log2(1.2)))
# tbl_BvC = lapply(tbl_BvC, function(u) 
#   filter(u, p_val < 0.05))

# #### Order ####
# 
# for(i in seq(1:length(tbl_BvC))){
#   print(dim(tbl_RvC[[i]]))
#   tbl_RvC[[i]] = tbl_RvC[[i]][order(tbl_RvC[[i]]$p_val),]
#   
# }
# for(i in seq(1:length(tbl_BvC))){
#   print(dim(tbl_BvC[[i]]))
#   tbl_BvC[[i]] = tbl_BvC[[i]][order(tbl_BvC[[i]]$p_val),]
# }
# 
# # Check overlap between red and blue OLIGroglia state DEA results
# overlap_OLIG_states = c()
# for(i in assayNames(pb_BvC)){
#   if(i %in% names(tbl_BvC) & i %in% names(tbl_RvC)){
#     print(i)
#     print(intersect(tbl_BvC[[i]]$gene, tbl_RvC[[i]]$gene))
#   }
# }
###

scRNA_OLIG_sce@colData$cluster_id = scRNA_OLIG_sce@colData$cell_type_high_resolution
scRNA_OLIG_sce@colData$group_id = pheno_ROSMAP[match(scRNA_OLIG_sce@colData$projid, pheno_ROSMAP$projid),]$diag
scRNA_OLIG_sce@colData$gender = pheno_ROSMAP[match(scRNA_OLIG_sce@colData$projid, pheno_ROSMAP$projid),]$msex
scRNA_OLIG_sce@colData$age_death = pheno_ROSMAP[match(scRNA_OLIG_sce@colData$projid, pheno_ROSMAP$projid),]$age_death
scRNA_OLIG_sce@colData$pmi = pheno_ROSMAP[match(scRNA_OLIG_sce@colData$projid, pheno_ROSMAP$projid),]$pmi
scRNA_OLIG_sce@colData$sample_id = paste(scRNA_OLIG_sce@colData$group_id, scRNA_OLIG_sce@colData$projid, sep = "_")
t(table(scRNA_OLIG_sce$cluster_id, scRNA_OLIG_sce$sample_id))
scRNA_OLIG_sce = prepSCE(scRNA_OLIG_sce)

pb_ADvC <- aggregateData(scRNA_OLIG_sce,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb_ADvC)
t(head(assay(pb_ADvC)))

#Filter by expr at the pseudobulk level
#Filter gene is by expression, filter samples is outlier detection
formula = ~gender+pmi+age_death
cd = as.data.frame(colData(pb_ADvC))
cd2=cd[,c("gender","pmi","age_death")]
design=model.matrix(formula,cd2)

cd[is.na(cd$pmi),]$pmi = mean(design[,3])
cd2=cd[,c("gender","pmi","age_death")]
design=model.matrix(formula,cd2)

res_ADvC <- pbDS(pb_ADvC, design = design, min_cells = 10, filter = "both", verbose = TRUE)
tbl_ADvC <- res_ADvC$table[[1]]


# 
# intersect_R = tbl_RvC[[i]][tbl_RvC[[i]]$gene %in% intersect(tbl_BvC[[i]]$gene, tbl_RvC[[i]]$gene),]
# intersect_R = intersect_R[order(intersect_R$gene),]
# intersect_B = tbl_BvC[[i]][tbl_BvC[[i]]$gene %in% intersect(tbl_BvC[[i]]$gene, tbl_RvC[[i]]$gene),]
# intersect_B = intersect_B[order(intersect_B$gene),]
# 
# intersect_RB = data.frame(intersect_R$gene, intersect_R$logFC, intersect_B$logFC)
# intersect_RB$dir_R = NA
# intersect_RB$dir_B = NA
# 
# for(i in seq(1:nrow(intersect_RB))){
#   if(intersect_RB$intersect_R.logFC[i] > 0){
#     intersect_RB$dir_R[i] = "+"
#   } else {
#   intersect_RB$dir_R[i] = "-"
#   }
#   
#   if(intersect_RB$intersect_B.logFC[i] > 0){
#     intersect_RB$dir_B[i] = "+"
#   } else {
#     intersect_RB$dir_B[i] = "-"
#   }
# }


###########################################

load("snRNA_Ne_Ol/res_excitatory_set1.Rdata")
res_ex1 = list(tbl_ADvC, tbl_BvC, tbl_RvC, pb_ADvC, pb_BvC, pb_RvC)

load("snRNA_Ne_Ol/res_excitatory_set2.Rdata")
res_ex2 = list(tbl_ADvC, tbl_BvC, tbl_RvC, pb_ADvC, pb_BvC, pb_RvC)

load("snRNA_Ne_Ol/res_excitatory_set3.Rdata")
res_ex3 = list(tbl_ADvC, tbl_BvC, tbl_RvC, pb_ADvC, pb_BvC, pb_RvC)

load("snRNA_Ne_Ol/res_inhibitory.Rdata")
res_inh = list(tbl_ADvC, tbl_BvC, tbl_RvC, pb_ADvC, pb_BvC, pb_RvC)

load("snRNA_Ne_Ol/res_scRNA_oligos.Rdata")
res_oli = list(tbl_ADvC, tbl_BvC, tbl_RvC, pb_ADvC, pb_BvC, pb_RvC)

remove(tbl_RvC, tbl_BvC, tbl_ADvC, pb_RvC, pb_BvC, pb_ADvC)

res_ex_all = list()
tmp_table = data.frame(matrix(ncol = 9, nrow = 0))
for(i in names(res_ex1[[1]])){
  tmp_table = rbind(tmp_table, res_ex1[[1]][[i]])
}
for(i in names(res_ex2[[1]])){
  tmp_table = rbind(tmp_table, res_ex2[[1]][[i]])
}
for(i in names(res_ex3[[1]])){
  tmp_table = rbind(tmp_table, res_ex3[[1]][[i]])
}

res_ex_all[["ADvC"]] = tmp_table

tmp_table = data.frame(matrix(ncol = 9, nrow = 0))
for(i in names(res_ex1[[2]])){
  tmp_table = rbind(tmp_table, res_ex1[[2]][[i]])
}
for(i in names(res_ex2[[2]])){
  tmp_table = rbind(tmp_table, res_ex2[[2]][[i]])
}
for(i in names(res_ex3[[2]])){
  tmp_table = rbind(tmp_table, res_ex3[[2]][[i]])
}

res_ex_all[["BvC"]] = tmp_table

tmp_table = data.frame(matrix(ncol = 9, nrow = 0))
for(i in names(res_ex1[[3]])){
  tmp_table = rbind(tmp_table, res_ex1[[3]][[i]])
}
for(i in names(res_ex2[[3]])){
  tmp_table = rbind(tmp_table, res_ex2[[3]][[i]])
}
for(i in names(res_ex3[[3]])){
  tmp_table = rbind(tmp_table, res_ex3[[3]][[i]])
}

res_ex_all[["RvC"]] = tmp_table

###

res_inh_all = list()
tmp_table = data.frame(matrix(ncol = 9, nrow = 0))
for(i in names(res_inh[[1]])){
  tmp_table = rbind(tmp_table, res_inh[[1]][[i]])
}

res_inh_all[["ADvC"]] = tmp_table

tmp_table = data.frame(matrix(ncol = 9, nrow = 0))
for(i in names(res_inh[[2]])){
  tmp_table = rbind(tmp_table, res_inh[[2]][[i]])
}


res_inh_all[["BvC"]] = tmp_table

tmp_table = data.frame(matrix(ncol = 9, nrow = 0))
for(i in names(res_inh[[3]])){
  tmp_table = rbind(tmp_table, res_inh[[3]][[i]])
}


res_inh_all[["RvC"]] = tmp_table


#########################

load("snRNA_Ne_Ol/cell_counts_exc1.Rdata")
load("snRNA_Ne_Ol/cell_counts_exc2.Rdata")
load("snRNA_Ne_Ol/cell_counts_exc3.Rdata")
load("snRNA_Ne_Ol/cell_counts_olig.Rdata")
load("snRNA_Ne_Ol/cell_counts_inh.Rdata")

table_cells_exc1_ADvC = as.data.frame(table_cells_exc1_ADvC)
table_cells_exc2_ADvC = as.data.frame(table_cells_exc2_ADvC)
table_cells_exc3_ADvC = as.data.frame(table_cells_exc3_ADvC)
table_cells_exc_ADvC = rbind(table_cells_exc1_ADvC, table_cells_exc2_ADvC,table_cells_exc3_ADvC)
table_cells_exc_ADvC = dcast(data = table_cells_exc_ADvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

table_cells_exc1_RvC = as.data.frame(table_cells_exc1_RvC)
table_cells_exc2_RvC = as.data.frame(table_cells_exc2_RvC)
table_cells_exc3_RvC = as.data.frame(table_cells_exc3_RvC)
table_cells_exc_RvC = rbind(table_cells_exc1_RvC, table_cells_exc2_RvC,table_cells_exc3_RvC)
table_cells_exc_RvC = dcast(data = table_cells_exc_RvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")


table_cells_exc1_BvC = as.data.frame(table_cells_exc1_BvC)
table_cells_exc2_BvC = as.data.frame(table_cells_exc2_BvC)
table_cells_exc3_BvC = as.data.frame(table_cells_exc3_BvC)
table_cells_exc_BvC = rbind(table_cells_exc1_BvC, table_cells_exc2_BvC,table_cells_exc3_BvC)
table_cells_exc_BvC = dcast(data = table_cells_exc_BvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

table_cells_inh_ADvC = as.data.frame(table_cells_inh_ADvC)
table_cells_inh_ADvC = dcast(data = table_cells_inh_ADvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")
table_cells_inh_RvC = as.data.frame(table_cells_inh_RvC)
table_cells_inh_RvC = dcast(data = table_cells_inh_RvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")
table_cells_inh_BvC = as.data.frame(table_cells_inh_BvC)
table_cells_inh_BvC = dcast(data = table_cells_inh_BvC,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

table_cells_olig_ADvC = as.data.frame(table_cells_olig_ADvC)
table_cells_olig_RvC = as.data.frame(table_cells_olig_RvC)
table_cells_olig_BvC = as.data.frame(table_cells_olig_BvC)

remove(table_cells_exc1_ADvC, table_cells_exc2_ADvC, table_cells_exc3_ADvC,
       table_cells_exc1_BvC, table_cells_exc2_BvC, table_cells_exc3_BvC,
       table_cells_exc1_RvC, table_cells_exc2_RvC, table_cells_exc3_RvC)

table_cells_exc_ADvC = table_cells_exc_ADvC %>% separate(Var1, c('diag', 'id'))
table_cells_exc_RvC = table_cells_exc_RvC %>% separate(Var1, c('diag', 'id'))
table_cells_exc_BvC = table_cells_exc_BvC %>% separate(Var1, c('diag', 'id'))


table_cells_exc_ADvC$`Exc L2-3 CBLN2 LINC02306`
boxplot(`Exc L2-3 CBLN2 LINC02306` ~ diag, data = table_cells_exc_ADvC)

table_cells_exc_all = rbind(table_cells_exc_ADvC,
                            table_cells_exc_RvC[table_cells_exc_RvC$diag == "red",],
                            table_cells_exc_BvC[table_cells_exc_BvC$diag == "blue",])

table_cells_inh_ADvC = table_cells_inh_ADvC %>% separate(Var1, c('diag', 'id'))
table_cells_inh_RvC = table_cells_inh_RvC %>% separate(Var1, c('diag', 'id'))
table_cells_inh_BvC = table_cells_inh_BvC %>% separate(Var1, c('diag', 'id'))

table_cells_inh_all = rbind(table_cells_inh_ADvC,
                            table_cells_inh_RvC[table_cells_inh_RvC$diag == "red",],
                            table_cells_inh_BvC[table_cells_inh_BvC$diag == "blue",])

table_cells_olig_ADvC = table_cells_olig_ADvC %>% separate(Var1, c('diag', 'id'))
table_cells_olig_RvC = table_cells_olig_RvC %>% separate(Var1, c('diag', 'id'))
table_cells_olig_BvC = table_cells_olig_BvC %>% separate(Var1, c('diag', 'id'))

table_cells_olig_all = rbind(table_cells_olig_ADvC,
                            table_cells_olig_RvC[table_cells_olig_RvC$diag == "red",],
                            table_cells_olig_BvC[table_cells_olig_BvC$diag == "blue",])
boxplot(Freq ~ diag, data = table_cells_olig_all) 

save(table_cells_exc_all, table_cells_inh_all, table_cells_olig_all, file = "snRNA_Ne_Ol/cell_counts_all.Rdata")
