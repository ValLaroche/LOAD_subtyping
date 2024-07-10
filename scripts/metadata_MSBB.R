MSBB <- function(){##### Detailing metadata #####

setwd("D:/valentin/main")
MSBB_all = read.csv("./MSBB/metadata/MSBB_assay_MethylationArray_metadata.csv") #709 samples
pheno_MSBB = read.csv("./MSBB/metadata/MSBB_individual_metadata.csv") #1221 samples

MSBB_filter = MSBB_all[MSBB_all %in% pheno_MSBB,]#604 samples

MSBB_metadata_filter = pheno_PFC[pheno_MSBB$BBNId %in% MSBB_filter$BBNId,]#1170 samples

MSBB_filter = MSBB_filter[order(MSBB_filter$BBNId),]
MSBB_metadata_filter = MSBB_metadata_filter[order(MSBB_metadata_filter$BBNId),]

MSBB_metadata_filter = merge(MSBB_metadata_filter, MSBB_filter, by = "BBNId", all.y = TRUE)
MSBB_metadata_filter = MSBB_metadata_filter[MSBB_metadata_filter$BR == "Prefrontal",] #584 samples

####
pheno_MSBB[pheno_MSBB$ageDeath == "90+",]$ageDeath = 90 #377
pheno_MSBB$ageDeath = as.numeric(pheno_MSBB$ageDeath)
table(pheno_MSBB$ageDeath)
pheno_MSBB = pheno_MSBB[pheno_MSBB$ageDeath > 65,] #361
pheno_MSBB$age_conc = NA ### Making age into simpler intervals
pheno_MSBB[pheno_MSBB$ageDeath < summary(pheno_MSBB$ageDeath)[2],]$age_conc = 1
pheno_MSBB[pheno_MSBB$ageDeath >= summary(pheno_MSBB$ageDeath)[2] & 
                      pheno_MSBB$ageDeath < summary(pheno_MSBB$ageDeath)[3],]$age_conc = 2
pheno_MSBB[pheno_MSBB$ageDeath >= summary(pheno_MSBB$ageDeath)[3] & 
                      pheno_MSBB$ageDeath < summary(pheno_MSBB$ageDeath)[5],]$age_conc = 3
pheno_MSBB[pheno_MSBB$ageDeath >= summary(pheno_MSBB$ageDeath)[5],]$age_conc = 4

pheno_MSBB$dementia = NA
pheno_MSBB$diag = NA
#CDR filter
table(pheno_MSBB$CDR)

#Dementia assessment
pheno_MSBB$dementia[pheno_MSBB$MMSE < 24] = 1
pheno_MSBB$dementia[pheno_MSBB$CDR >= 1] = 1
pheno_MSBB$dementia[pheno_MSBB$CDR < 1] = 0
pheno_MSBB$dementia[is.na(pheno_MSBB$dementia)] = 2# Anything we can't conclude
table(pheno_MSBB$dementia)
#2 should be removed later on (table)
#0   1 
#91 270
#Diagnosis
pheno_MSBB$diag[pheno_MSBB$CERAD <= 1 & pheno_MSBB$Braak <= 2 & pheno_MSBB$dementia == 0] = "Control" 
pheno_MSBB$diag[pheno_MSBB$CERAD == 0 & pheno_MSBB$Braak <= 3 & pheno_MSBB$dementia == 0] = "Control" 
pheno_MSBB$diag[pheno_MSBB$CERAD >= 1 & pheno_MSBB$Braak >= 3 & pheno_MSBB$dementia == 0] = "AsymAD"
pheno_MSBB$diag[pheno_MSBB$CERAD >= 2 & pheno_MSBB$Braak >= 3 & pheno_MSBB$dementia == 1] = "AD" 
pheno_MSBB$diag[is.na(pheno_MSBB$diag)] = "Unassigned"

#AD     AsymAD    Control Unassigned 
#216         37         50         58

##### Curating dataset #####
MSBB_ok = MSBB_metadata_filter[MSBB_metadata_filter$diag != "Unassigned",]

# Clinical AD
#We searched for clinical diagnoses of AD made by doctors directly on unassigned samples to get more samples
MSBB_unassigned = MSBB_metadata_filter[MSBB_metadata_filter$diag == "Unassigned",]
MSBB_unassigned[grep("Alzheimer", MSBB_unassigned$clinical),]$clinical
MSBB_unassigned[grep("Alzheimer", MSBB_unassigned$clinical),]$diag = "AD"
MSBB_ok = rbind(MSBB_ok, MSBB_unassigned[MSBB_unassigned$diag != "Unassigned",])
MSBB_unassigned = MSBB_unassigned[MSBB_unassigned$diag == "Unassigned",]
# We move from 217 AD to 314 AD

# Clinical dementia
MSBB_unassigned[grep("dementia", MSBB_unassigned$clinical),]$clinical
MSBB_unassigned[grep("dementia", MSBB_unassigned$clinical),]$dementia = 1

MSBB_unassigned$diag[MSBB_unassigned$Cerad <= 1 & MSBB_unassigned$Braak_tangle <= 2 & MSBB_unassigned$dementia == 0] = "Control" 
MSBB_unassigned$diag[MSBB_unassigned$Cerad == 0 & MSBB_unassigned$Braak_tangle <= 3 & MSBB_unassigned$dementia == 0] = "Control" 
MSBB_unassigned$diag[MSBB_unassigned$Cerad >= 1 & MSBB_unassigned$Braak_tangle >= 3 & MSBB_unassigned$dementia == 0] = "AsymAD"
MSBB_unassigned$diag[MSBB_unassigned$Cerad >= 2 & MSBB_unassigned$Braak_tangle >= 3 & MSBB_unassigned$dementia == 1] = "AD" 
table(MSBB_unassigned$diag)
MSBB_ok = rbind(MSBB_ok, MSBB_unassigned[MSBB_unassigned$diag != "Unassigned",])
MSBB_unassigned = MSBB_unassigned[MSBB_unassigned$diag == "Unassigned",]
#313 AD to 323 AD

# Braak tangle and Cerad NA's can't be fixed - no data
# No abnormalities = controls
MSBB_unassigned[grep("No abnormality detected", MSBB_unassigned$clinical),]$clinical
MSBB_unassigned[grep("No abnormality detected", MSBB_unassigned$clinical),]$diag = "Control"
MSBB_ok = rbind(MSBB_ok, MSBB_unassigned[MSBB_unassigned$diag != "Unassigned",])
MSBB_unassigned = MSBB_unassigned[MSBB_unassigned$diag == "Unassigned",]
#104 Controls to 144 controls



# Rest
MSBB_unassigned$clinical
MSBB_unassigned$clinical_plus
MSBB_unassigned$pathology
MSBB_unassigned$Braak_LB

# Final = We exclude MSBB_unassigned 83 remaining samples because not enough information OR other types of Dementia
# Final = We exclude the 35 AsymAD samples and age <65
#313 AD and 140 Controls

MSBB_ok_age = MSBB_ok[MSBB_ok$Age >= 65,]
table(MSBB_ok_age$diag)
MSBB_ok_ADCtrl = MSBB_ok_age[MSBB_ok_age$diag != "AsymAD",]

table(MSBB_ok_ADCtrl$diag)
table(MSBB_ok_ADCtrl$dementia)

table(MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "AD",]$MMSE)
table(MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "AD",]$CDR)
table(MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "AD",]$Braak_tangle)
table(MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "AD",]$Cerad)

table(MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "AD" & MSBB_ok_ADCtrl$CDR < 1,]$MMSE)



##### Curating Controls #####
MSBB_ok_ADCtrl = MSBB_ok_ADCtrl[order(MSBB_ok_ADCtrl$BBNId),]
MSBB_ok_ADCtrl$BBNId

MSBB_ok_ADCtrl$Thal_number = NA
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == "Thal phase 0 "] = 0
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == "Thal phase 1 "] = 1
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == "Thal phase 2 "] = 2
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == "Thal phase 3 "] = 3
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == "Thal phase 4 "] = 4
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == "Thal phase 5 "] = 5
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == "Thal phase not assessed "] = NA
MSBB_ok_ADCtrl$Thal_number[MSBB_ok_ADCtrl$Thal.stage == ""] = NA
table(MSBB_ok_ADCtrl$Thal.stage)
table(MSBB_ok_ADCtrl$Thal_number)

MSBB_Ctrl = MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "Control",]
MSBB_Ctrl$new_diag = NA
MSBB_Ctrl$Thal.stage


#Braak stage 0-2
table(MSBB_Ctrl$BraakTangle_numeric)
MSBB_Ctrl_BTnonOK = MSBB_Ctrl[which(MSBB_Ctrl$Braak_tangle >= 4),]

# Argyrophilic
grep("+Ndgen.NonADTau.AGD", MSBB_Ctrl$pathology, fixed = TRUE)
MSBB_Ctrl_Argnonok = MSBB_Ctrl[grep("+Ndgen.NonADTau.AGD", MSBB_Ctrl$pathology, fixed = TRUE),]
MSBB_Ctrl_Argnonok$pathology  

# Thal and Cerad 0
table(MSBB_Ctrl$Cerad, MSBB_Ctrl$Thal_number)
MSBB_Ctrl_amynonok = MSBB_Ctrl[which(MSBB_Ctrl$Cerad > 0 | MSBB_Ctrl$Thal_number > 0),]

# TDP43
grep("+Ndgen.TDP.TDPaccomp", MSBB_Ctrl$pathology, fixed = TRUE)
MSBB_Ctrl_TDP43 = MSBB_Ctrl[grep("+Ndgen.TDP.TDPaccomp", MSBB_Ctrl$pathology, fixed = TRUE),]
MSBB_Ctrl_TDP43$pathology

# Alpha synuclein - no info ??
grep("a-", MSBB_Ctrl$pathology, fixed = TRUE)
MSBB_Ctrl_alpha = MSBB_Ctrl[grep("a-", MSBB_Ctrl$pathology, fixed = TRUE),]
MSBB_Ctrl_alpha$pathology

# Cerebrovascular

grep("vascular", MSBB_Ctrl$clinical)
grep("vascular", MSBB_Ctrl$clinical_plus)
grep("vascular", MSBB_Ctrl$pathology)

final_cereb = c(grep("vascular", MSBB_Ctrl$clinical, fixed = TRUE),
                grep("vascular", MSBB_Ctrl$clinical_plus, fixed = TRUE),
                grep("+CVD.Vasc.Athero.Sevathero", MSBB_Ctrl$pathology, fixed = TRUE))
MSBB_Ctrl_cereb = MSBB_Ctrl[unique(final_cereb),]
MSBB_Ctrl_cereb$clinical
MSBB_Ctrl_cereb$clinical_plus
MSBB_Ctrl_cereb$pathology

# Leukemia - None ?
grep("LEUKEMIA", MSBB_Ctrl$Cause_of_Death)

# Neuropathologic correlate - None ?
grep("psychia", MSBB_Ctrl$clinical_plus)


# Final assessment
final_ctrl_nonok = c(grep("vascular", MSBB_Ctrl$clinical, fixed = TRUE),
                     grep("vascular", MSBB_Ctrl$clinical_plus, fixed = TRUE),
                     grep("+CVD.Vasc.Athero.Sevathero", MSBB_Ctrl$pathology, fixed = TRUE),
                     grep("+Ndgen.TDP.TDPaccomp", MSBB_Ctrl$pathology, fixed = TRUE),
                     which(MSBB_Ctrl$Cerad > 1 | MSBB_Ctrl$Thal_number > 2),
                     grep("+Ndgen.NonADTau.AGD", MSBB_Ctrl$pathology, fixed = TRUE),
                     which(MSBB_Ctrl$Braak_tangle >= 3))
final_ctrl_nonok = unique(final_ctrl_nonok)


MSBB_remove = MSBB_Ctrl[-final_ctrl_nonok,]
MSBB_remove$Basename

##### Train / Test models #####
#### Normal sampling #### 
# training_ID = splitstackshape::stratified(MSBB_metadata_AD, c("age_conc","Gender"),size = 0.6)
# training_ID = training_ID$Basename
# test_ID = MSBB_metadata_AD[!MSBB_metadata_AD$Basename %in% training_ID,]$Basename
#### Kstone sampling ####

#### Split AD Ctrl ####

MSBB_ok_ADCtrl = MSBB_ok_ADCtrl[!MSBB_ok_ADCtrl$BBNId %in% MSBB_remove$BBNId,]
table(MSBB_ok_ADCtrl$diag)

save(MSBB_ok_ADCtrl, file = "C:/Users/p70077107/Desktop/MSBB_samples_final_verified_25-8-23.rData")


MSBB_ok_ADCtrl = MSBB_ok_ADCtrl[order(MSBB_ok_ADCtrl$Basename),]
mSet_betas = mSet_betas[,order(colnames(mSet_betas))]
head(colnames(mSet_betas))
head(MSBB_ok_ADCtrl$Basename)
table(MSBB_ok_ADCtrl$diag)

MSBB_ok_ADCtrl$age_conc = NA ### Making age into simpler intervals
MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$Age < summary(MSBB_ok_ADCtrl$Age)[2],]$age_conc = 1
MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$Age >= summary(MSBB_ok_ADCtrl$Age)[2] & 
                MSBB_ok_ADCtrl$Age < summary(MSBB_ok_ADCtrl$Age)[3],]$age_conc = 2
MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$Age >= summary(MSBB_ok_ADCtrl$Age)[3] & 
                MSBB_ok_ADCtrl$Age < summary(MSBB_ok_ADCtrl$Age)[5],]$age_conc = 3
MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$Age >= summary(MSBB_ok_ADCtrl$Age)[5],]$age_conc = 4

MSBB_metadata_AD = MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "AD",]
MSBB_metadata_Ctrl = MSBB_ok_ADCtrl[MSBB_ok_ADCtrl$diag == "Control",]
mSet_betas_AD = mSet_betas[,colnames(mSet_betas) %in% MSBB_metadata_AD$Basename]
mSet_betas_Ctrl = mSet_betas[,colnames(mSet_betas) %in% MSBB_metadata_Ctrl$Basename]
}
