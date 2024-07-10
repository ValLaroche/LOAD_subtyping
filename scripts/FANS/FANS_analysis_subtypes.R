##### FACS celltype enrichment #####
# Load FACS EWAS results
load("FACS/IRF8_allcpg.rdata")
AllLME_IRF8 = AllLME
load("FACS/NeuN_allcpg.rdata")
AllLME_NeuN = AllLME
load("FACS/Sox10_allcpg.rdata")
AllLME_Sox10 = AllLME
load("FACS/Trip neg_allcpg.rdata")
AllLME_TripNeg = AllLME
remove(AllLME)

AllLME_IRF8 = t(AllLME_IRF8)
AllLME_NeuN = t(AllLME_NeuN)
AllLME_Sox10 = t(AllLME_Sox10)
AllLME_TripNeg = t(AllLME_TripNeg)

head(AllLME_IRF8)
head(AllLME_NeuN)
head(AllLME_Sox10)
head(AllLME_TripNeg)

identical(rownames(AllLME_IRF8), rownames(AllLME_Sox10))

# Concatenate EWAS results
All_CT = cbind(AllLME_IRF8[,c(1,5)], AllLME_Sox10[,c(1,5)], AllLME_NeuN[,c(1,5)], AllLME_TripNeg[,c(1,5)])
head(All_CT)
colnames(All_CT) = c("Mic_est", "Mic_pval",
                     "Olig_est", "Olig_pval",
                     "NeuN_est", "NeuN_pval",
                     "Ast_est", "Ast_pval")
All_CT = All_CT[rownames(All_CT) %in% rownames(betas_ROSMAP),]
All_CT = as.data.frame(All_CT)

All_CT$Mic_pval = p.adjust(All_CT$Mic_pval, method = "bonf")
range(All_CT$Mic_est)

All_CT$Ast_pval = p.adjust(All_CT$Ast_pval, method = "bonf")
range(All_CT$Ast_est)

All_CT$NeuN_pval = p.adjust(All_CT$NeuN_pval, method = "bonf")
range(All_CT$NeuN_est)

All_CT$Olig_pval = p.adjust(All_CT$Olig_pval, method = "bonf")
range(All_CT$Olig_est)

All_CT_blue = All_CT[rownames(All_CT) %in% rownames(sig_methyl_450k_blueC),]

All_CT_red = All_CT[rownames(All_CT) %in% rownames(sig_methyl_450k_redC),]

All_CT_grey = All_CT[rownames(All_CT) %in% rownames(sig_methyl_450k_greyC),]

# save(All_CT, All_CT_red, All_CT_blue, All_CT_AD, file = "FACS/All_LME_FACS450k.rdata")
load("FACS/All_LME_FACS450k.rdata")

# Filter for significance for each cell type
pvalue_threshold = 1e-10
effect_threshold = 0.5

AllLME_IRF8 = as.data.frame(AllLME_IRF8)
AllLME_IRF8 = AllLME_IRF8[AllLME_IRF8$"p-value" < pvalue_threshold & abs(AllLME_IRF8$Value) > effect_threshold,]
AllLME_IRF8 = AllLME_IRF8[rownames(AllLME_IRF8) %in% rownames(betas_ROSMAP),]
AllLME_IRF8 = AllLME_IRF8[order(AllLME_IRF8$`p-value`),]
AllLME_IRF8$p.bonf = p.adjust(AllLME_IRF8$`p-value`, method = "bonf")
write.csv(AllLME_IRF8, file = "SUBTYPING-PAPER-TABLEFIG/Sig_MIC_FANS.csv")

AllLME_NeuN = as.data.frame(AllLME_NeuN)
AllLME_NeuN = AllLME_NeuN[AllLME_NeuN$"p-value" < pvalue_threshold & abs(AllLME_NeuN$Value) > effect_threshold,]
AllLME_NeuN = AllLME_NeuN[rownames(AllLME_NeuN) %in% rownames(betas_ROSMAP),]
AllLME_NeuN = AllLME_NeuN[order(AllLME_NeuN$`p-value`),]
AllLME_NeuN$p.bonf = p.adjust(AllLME_NeuN$`p-value`, method = "bonf")
write.csv(AllLME_NeuN, file = "SUBTYPING-PAPER-TABLEFIG/Sig_NEUN_FANS.csv")

AllLME_Sox10 = as.data.frame(AllLME_Sox10)
AllLME_Sox10 = AllLME_Sox10[AllLME_Sox10$"p-value" < pvalue_threshold & abs(AllLME_Sox10$Value) > effect_threshold,]
AllLME_Sox10 = AllLME_Sox10[rownames(AllLME_Sox10) %in% rownames(betas_ROSMAP),]
AllLME_Sox10 = AllLME_Sox10[order(AllLME_Sox10$`p-value`),]
AllLME_Sox10$p.bonf = p.adjust(AllLME_Sox10$`p-value`, method = "bonf")
write.csv(AllLME_Sox10, file = "SUBTYPING-PAPER-TABLEFIG/Sig_OLIG_FANS.csv")

AllLME_TripNeg = as.data.frame(AllLME_TripNeg)
AllLME_TripNeg = AllLME_TripNeg[AllLME_TripNeg$"p-value" < pvalue_threshold & abs(AllLME_TripNeg$Value) > effect_threshold,]
AllLME_TripNeg = AllLME_TripNeg[rownames(AllLME_TripNeg) %in% rownames(betas_ROSMAP),]
AllLME_TripNeg = AllLME_TripNeg[order(AllLME_TripNeg$`p-value`),]
AllLME_TripNeg$p.bonf = p.adjust(AllLME_TripNeg$`p-value`, method = "bonf")
write.csv(AllLME_TripNeg, file = "SUBTYPING-PAPER-TABLEFIG/Sig_AST_FANS.csv")

# Perform over-representation tests
table_res_fisher = data.frame(matrix(ncol = 6, nrow = 4))
colnames(table_res_fisher) = c("Blue_pval","Blue_odds", "Red_pval","Red_odds", "Grey_pval", "Grey_odds")
rownames(table_res_fisher) = c("NeuN","Ast","Mic","Oligo")

####  blue ####
table_blue_Olig = matrix(c(table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                           table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[1] - table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2] - table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1] + table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2]),
                         nrow = 2)
fisher_blue_Olig = fisher.test(table_blue_Olig, alternative = "greater")
table_res_fisher[4,1] = fisher_blue_Olig$p.value
table_res_fisher[4,2] = fisher_blue_Olig$estimate

table_blue_NeuN = matrix(c(table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[1] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1] + table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2]),
                         nrow = 2)
fisher_blue_NeuN = fisher.test(table_blue_NeuN, alternative = "greater")
table_res_fisher[1,1] = fisher_blue_NeuN$p.value
table_res_fisher[1,2] = fisher_blue_NeuN$estimate

table_blue_Ast = matrix(c(table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[2],
                          table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[1] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2] - table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1] + table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2]),
                        nrow = 2)
fisher_blue_Ast = fisher.test(table_blue_Ast, alternative = "greater")
table_res_fisher[2,1] = fisher_blue_Ast$p.value
table_res_fisher[2,2] = fisher_blue_Ast$estimate

table_blue_Mic = matrix(c(table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                          table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[1] - table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2] - table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1]),
                        nrow = 2)
fisher_blue_Mic = fisher.test(table_blue_Mic, alternative = "greater")
table_res_fisher[3,1] = fisher_blue_Mic$p.value
table_res_fisher[3,2] = fisher_blue_Mic$estimate

#### red ####
table_red_Olig = matrix(c(table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                          table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[1] - table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                          table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2] - table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                          table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1]),
                        nrow = 2)
fisher_red_Olig = fisher.test(table_red_Olig, alternative = "greater")
table_res_fisher[4,3] = fisher_red_Olig$p.value
table_res_fisher[4,4] = fisher_red_Olig$estimate

table_red_NeuN = matrix(c(table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                          table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[1] - table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                          table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2] - table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                          table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1]),
                        nrow = 2)
fisher_red_NeuN = fisher.test(table_red_NeuN, alternative = "greater")
table_res_fisher[1,3] = fisher_red_NeuN$p.value
table_res_fisher[1,4] = fisher_red_NeuN$estimate

table_red_Ast = matrix(c(table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[2],
                         table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[1] - table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[2],
                         table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2] - table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[2],
                         table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1]),
                       nrow = 2)
# fisher_red_Ast = fisher.test(table_red_Ast, alternative = "greater")
table_res_fisher[2,3] = 0
table_res_fisher[2,4] = 0

table_red_Mic = matrix(c(table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                         table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[1] - table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                         table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2] - table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                         table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1]),
                       nrow = 2)
fisher_red_Mic = fisher.test(table_red_Mic, alternative = "greater")
table_res_fisher[3,3] = fisher_red_Mic$p.value
table_res_fisher[3,4] = fisher_red_Mic$estimate

#### grey ####
table_grey_Olig = matrix(c(table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2],
                           table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[1] - table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2] - table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2],
                           table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1]),
                         nrow = 2)
fisher_grey_Olig = fisher.test(table_grey_Olig, alternative = "greater")
table_res_fisher[4,5] = fisher_grey_Olig$p.value
table_res_fisher[4,6] = fisher_grey_Olig$estimate

table_grey_NeuN = matrix(c(table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                           table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[1] - table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2] - table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                           table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1]),
                         nrow = 2)
fisher_grey_NeuN = fisher.test(table_grey_NeuN, alternative = "greater")
table_res_fisher[1,5] = fisher_grey_NeuN$p.value
table_res_fisher[1,6] = fisher_grey_NeuN$estimate

table_grey_Ast = matrix(c(table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[2],
                          table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[1] - table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2] - table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[2],
                          table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1]),
                        nrow = 2)
# fisher_grey_Ast = fisher.test(table_grey_Ast, alternative = "greater")
table_res_fisher[2,5] = 0
table_res_fisher[2,6] = 0

table_grey_Mic = matrix(c(table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                          table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[1] - table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2] - table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                          table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1]),
                        nrow = 2)
fisher_grey_Mic = fisher.test(table_grey_Mic)
table_res_fisher[3,5] = fisher_grey_Mic$p.value
table_res_fisher[3,6] = fisher_grey_Mic$estimate



# Generate barplot of the results

to_barplot = data.frame(matrix(ncol = 8, nrow = 12))
colnames(to_barplot) = c("subtype","cell_type","total_non_sig","total_sig",
                         "celltype_sig","celltype_non_sig","pvalue","odds")
to_barplot$subtype = c(rep("blue", 4), rep("red", 4), rep("grey", 4))
to_barplot$cell_type = rep(c("NeuN", "Ast", "Mic", "Olig"), 3)

to_barplot$total_non_sig = c(table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1],
                             table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1],
                             table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1],
                             table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1],
                             table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1],
                             table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1],
                             table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1],
                             table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1],
                             table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[1],
                             table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[1],
                             table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[1],
                             table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[1])

to_barplot$total_sig = c(table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2],
                         table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2],
                         table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2],
                         table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2],
                         table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2],
                         table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2],
                         table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2],
                         table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2],
                         table(All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold)[2],
                         table(All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold)[2],
                         table(All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold)[2],
                         table(All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold)[2])

to_barplot$celltype_sig = c(table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[2],
                            table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[2],
                            table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[2],
                            table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[2],
                            table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[2],
                            0,
                            table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[2],
                            table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[2],
                            table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[2],
                            0,
                            table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[2],
                            table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[2])

to_barplot$celltype_non_sig = c(table(All_CT_blue$NeuN_pval < pvalue_threshold & abs(All_CT_blue$NeuN_est) > effect_threshold)[1],
                                table(All_CT_blue$Ast_pval < pvalue_threshold & abs(All_CT_blue$Ast_est) > effect_threshold)[1],
                                table(All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold)[1],
                                table(All_CT_blue$Olig_pval < pvalue_threshold & abs(All_CT_blue$Olig_est) > effect_threshold)[1],
                                table(All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold)[1],
                                table(All_CT_red$Ast_pval < pvalue_threshold & abs(All_CT_red$Ast_est) > effect_threshold)[1],
                                table(All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold)[1],
                                table(All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold)[1],
                                table(All_CT_grey$NeuN_pval < pvalue_threshold & abs(All_CT_grey$NeuN_est) > effect_threshold)[1],
                                table(All_CT_grey$Ast_pval < pvalue_threshold & abs(All_CT_grey$Ast_est) > effect_threshold)[1],
                                table(All_CT_grey$Mic_pval < pvalue_threshold & abs(All_CT_grey$Mic_est) > effect_threshold)[1],
                                table(All_CT_grey$Olig_pval < pvalue_threshold & abs(All_CT_grey$Olig_est) > effect_threshold)[1])

to_barplot$pvalue = c(table_res_fisher$Blue_pval, table_res_fisher$Red_pval, table_res_fisher$Grey_pval)
to_barplot$odds = c(table_res_fisher$Blue_odds, table_res_fisher$Red_odds, table_res_fisher$Grey_odds)
to_barplot$total_cpg = rowSums(to_barplot[,3:4])
to_barplot$total_subtype = rowSums(to_barplot[,5:6])

to_barplot$total_perc_sig = to_barplot$total_sig / to_barplot$total_cpg
to_barplot$total_perc_nonsig = to_barplot$total_non_sig / to_barplot$total_cpg
to_barplot$celltype_perc_sig = to_barplot$celltype_sig / to_barplot$total_subtype
to_barplot$celltype_perc_nonsig = to_barplot$celltype_non_sig / to_barplot$total_subtype

to_barplot_percs = to_barplot[,c(1,2,11,12,13,14)]

to_barplot_melt = melt(to_barplot_percs)
to_barplot_melt$cpg_level = str_split_fixed(as.character(to_barplot_melt$variable), "_",2)[,1]
to_barplot_melt$subtype = as.character(to_barplot_melt$subtype)
to_barplot_melt[to_barplot_melt$cpg_level == "celltype",]$cpg_level = to_barplot_melt[to_barplot_melt$cpg_level == "celltype",]$subtype

to_barplot_melt = to_barplot_melt[-(5:12),]

to_barplot_pvalue = to_barplot[,c(1,2,7)]
to_barplot_melt = merge(to_barplot_melt, to_barplot_pvalue)

to_barplot_melt[to_barplot_melt$variable %in% c("total_perc_sig","total_perc_nonsig","celltype_perc_nonsig"),]$pvalue = NA

to_barplot_melt$pvalue = as.numeric(format(to_barplot_melt$pvalue, digits = 3, scientific = TRUE))

to_barplot_melt$significant = factor(str_split_fixed(as.character(to_barplot_melt$variable), "_",3)[,3],
                                     levels = c("sig", "nonsig"))

to_barplot_melt$cell_type = factor(to_barplot_melt$cell_type, levels = c("NeuN","Olig","Mic","Ast"))

to_barplot_melt = to_barplot_melt[to_barplot_melt$significant == "sig",]
to_barplot_melt[to_barplot_melt$cpg_level == "total",]$cpg_level = "Complete set"
to_barplot_melt[to_barplot_melt$cpg_level == "red",]$cpg_level = "RED subtype"
to_barplot_melt[to_barplot_melt$cpg_level == "blue",]$cpg_level = "BLUE subtype"
to_barplot_melt[to_barplot_melt$cpg_level == "grey",]$cpg_level = "Unassigned AD"

to_barplot_melt$cpg_level = factor(to_barplot_melt$cpg_level, levels = c("Complete set", "BLUE subtype", "RED subtype", "Unassigned AD"))

to_barplot_melt[9,]$pvalue = NA
to_barplot_melt[13,]$pvalue = NA
to_barplot_melt[is.na(to_barplot_melt$pvalue),]$pvalue = 1

to_barplot_melt = to_barplot_melt[to_barplot_melt$cell_type != "Ast",]

to_barplot_melt$asterisks = NA

ggplot(to_barplot_melt, aes(x = cell_type, y = value, fill = cpg_level)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = scales::percent(value, 0.1)), position = position_dodge(width = 0.9), vjust = -0.5) +
  geom_text(aes(label = asterisks), position = position_dodge(width = 0.9), vjust = -0.5, size = 8) +
  # geom_text(aes(y = 0.15, label = scales::scientific(pvalue)), position = position_dodge(width = 0.9), vjust = -0.5) +
  # geom_text(aes(y = 0.16, label = if_else(cpg_level == "total", "Fisher's test P-value :", ""))) +
  geom_text(aes(y = -0.007, label =cpg_level ), position = position_dodge(width = 0.9), vjust = -0.5, angle = 45, size = 3) +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15)) + 
  scale_x_discrete(labels=c("NeuN" = "Neuronal", "Olig" = "Oligodendrocytes",
                            "Mic" = "Microglia","Ast" = "Astrocytes")) +
  xlab("Cell type") + ylab("Percentage of significance in the CpG set") +
  # geom_hline(yintercept=1.05) + geom_hline(yintercept=1.15) +
  scale_fill_discrete(type = c("black", "steelblue1", "brown2", "grey"), name = "CpG set") + 
  theme_bw() +
  theme(axis.text.x = element_text(size=10))


# Get significant CpGs for cell types
BG_Mic = All_CT[All_CT$Mic_pval < pvalue_threshold & abs(All_CT$Mic_est) > effect_threshold,]
BG_Neun = All_CT[All_CT$NeuN_pval < pvalue_threshold & abs(All_CT$NeuN_est) > effect_threshold,]
BG_Olig = All_CT[All_CT$Olig_pval < pvalue_threshold & abs(All_CT$Olig_est) > effect_threshold,]
BG_Ast = All_CT[All_CT$Ast_pval < pvalue_threshold & abs(All_CT$Ast_est) > effect_threshold,]
write.csv(BG_Mic, file ="SUBTYPING-PAPER-TABLEFIG/Sig_MIC_FANS.csv")
write.csv(BG_Neun, file ="SUBTYPING-PAPER-TABLEFIG/Sig_NEUN_FANS.csv")
write.csv(BG_Olig, file ="SUBTYPING-PAPER-TABLEFIG/Sig_OLIG_FANS.csv")
write.csv(BG_Ast, file ="SUBTYPING-PAPER-TABLEFIG/Sig_AST_FANS.csv")

# Perform subtype-specific downstream for each cell type based on overrepresentation
## MIC
All_CT_red_mic = All_CT_red[All_CT_red$Mic_pval < pvalue_threshold & abs(All_CT_red$Mic_est) > effect_threshold,]
All_CT_blue_mic = All_CT_blue[All_CT_blue$Mic_pval < pvalue_threshold & abs(All_CT_blue$Mic_est) > effect_threshold,]
intersect(rownames(All_CT_blue_mic), rownames(All_CT_red_mic))
FANS_Gometh_Red = gometh(sig.cpg = rownames(All_CT_red_mic),
                         all.cpg = rownames(BG_Mic$X),
                         collection = c("GO", "KEGG"),
                         array.type = c("450K", "EPIC"),
                         sig.genes = FALSE)
FANS_Gometh_Red = FANS_Gometh_Red[FANS_Gometh_Red$P.DE < 0.05,]
FANS_Gometh_Red = FANS_Gometh_Red[order(FANS_Gometh_Red$P.DE),]
write.csv(FANS_Gometh_Red, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Red_MIC.csv")

FANS_Gometh_Blue = gometh(sig.cpg = rownames(All_CT_blue_mic),
                          all.cpg = rownames(BG_Mic$X),
                          collection = c("GO", "KEGG"),
                          array.type = c("450K", "EPIC"),
                          sig.genes = FALSE)
FANS_Gometh_Blue = FANS_Gometh_Blue[FANS_Gometh_Blue$P.DE < 0.05,]
FANS_Gometh_Blue = FANS_Gometh_Blue[order(FANS_Gometh_Blue$P.DE),]
write.csv(FANS_Gometh_Blue, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Blue_MIC.csv")

## NeuN RED
All_CT_red_neun = All_CT_red[All_CT_red$NeuN_pval < pvalue_threshold & abs(All_CT_red$NeuN_est) > effect_threshold,]
FANS_Gometh_Red = gometh(sig.cpg = rownames(All_CT_red_neun),
                         all.cpg = rownames(BG_Neun$X),
                         collection = c("GO", "KEGG"),
                         array.type = c("450K", "EPIC"),
                         sig.genes = FALSE)
FANS_Gometh_Red = FANS_Gometh_Red[FANS_Gometh_Red$P.DE < 0.05,]
FANS_Gometh_Red = FANS_Gometh_Red[order(FANS_Gometh_Red$P.DE),]
write.csv(FANS_Gometh_Red, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Red_NEUN.csv")



## Oligo RED
All_CT_red_olig = All_CT_red[All_CT_red$Olig_pval < pvalue_threshold & abs(All_CT_red$Olig_est) > effect_threshold,]
FANS_Gometh_Red = gometh(sig.cpg = rownames(All_CT_red_olig),
                         all.cpg = rownames(BG_Olig$X),
                         collection = c("GO", "KEGG"),
                         array.type = c("450K", "EPIC"),
                         sig.genes = FALSE)
FANS_Gometh_Red = FANS_Gometh_Red[FANS_Gometh_Red$P.DE < 0.05,]
FANS_Gometh_Red = FANS_Gometh_Red[order(FANS_Gometh_Red$P.DE),]
write.csv(FANS_Gometh_Red, file ="SUBTYPING-PAPER-TABLEFIG/FANS_Gometh_Red_OLIG.csv")

