block <- function(){
##### Detailing metadata #####

setwd("D:/valentin/main")
BDR_all = read.csv("raw_data/BDR_FULL_PHENOTYPING.csv", sep = ";") #709 samples
pheno_PFC = read.csv("./BDR_pheno.csv", sep = "\t") #1221 samples

BDR_filter = BDR_all[BDR_all$BBNId %in% pheno_PFC$BBNId,]#604 samples

BDR_metadata_filter = pheno_PFC[pheno_PFC$BBNId %in% BDR_filter$BBNId,]#1170 samples

BDR_filter = BDR_filter[order(BDR_filter$BBNId),]
BDR_metadata_filter = BDR_metadata_filter[order(BDR_metadata_filter$BBNId),]

BDR_metadata_filter = merge(BDR_metadata_filter, BDR_filter, by = "BBNId", all.y = TRUE)
BDR_metadata_filter = BDR_metadata_filter[BDR_metadata_filter$BR == "Prefrontal",] #584 samples
#To 584 because of either missing metadata OR no overlapping PFC/Occipital sample

BDR_metadata_filter$Clinical.Diagnoses = str_conv(BDR_metadata_filter$Clinical.Diagnoses, encoding = "UTF-8")
BDR_metadata_filter$Pathology = str_conv(BDR_metadata_filter$Pathology, encoding = "UTF-8")
BDR_metadata_filter$Clin.Diagnosis.Details = str_conv(BDR_metadata_filter$Clin.Diagnosis.Details, encoding = "UTF-8")
# BDR_metadata_filter$Clinical.Diagnoses = make.names(BDR_metadata_filter$Clinical.Diagnoses)
# BDR_metadata_filter$Clinical.Diagnoses = gsub("."," ", BDR_metadata_filter$Clinical.Diagnoses)

BDR_metadata_filter$Braak_tangle = BDR_metadata_filter$BraakTangle_numeric
BDR_metadata_filter$Braak_LB = BDR_metadata_filter$Braak.LB.stage
BDR_metadata_filter$Cerad = BDR_metadata_filter$Cerad
BDR_metadata_filter$CDR = BDR_metadata_filter$Final.CDR
BDR_metadata_filter$Gender = BDR_metadata_filter$Gender.x
BDR_metadata_filter$Age = BDR_metadata_filter$Age.x
BDR_metadata_filter$Cause_of_Death = BDR_metadata_filter$Cause.Of.Death
BDR_metadata_filter$MMSE = BDR_metadata_filter$Final.MMSE
BDR_metadata_filter$APOE = BDR_metadata_filter$Polymorphisms
BDR_metadata_filter$clinical = BDR_metadata_filter$Clinical.Diagnoses
BDR_metadata_filter$clinical_plus = BDR_metadata_filter$Clin.Diagnosis.Details
BDR_metadata_filter$pathology = BDR_metadata_filter$Pathology


BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD no plaques "] = 0
BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD sparse neuritic plaques "] = 1
BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD moderate density of neuritic plaques "] = 2
BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == "CERAD high density of neuritic plaques "] = 3
BDR_metadata_filter$Cerad[BDR_metadata_filter$Cerad == ""] = NA

BDR_metadata_filter$age_conc = NA ### Making age into simpler intervals
BDR_metadata_filter[BDR_metadata_filter$Age < summary(BDR_metadata_filter$Age)[2],]$age_conc = 1
BDR_metadata_filter[BDR_metadata_filter$Age >= summary(BDR_metadata_filter$Age)[2] & 
                      BDR_metadata_filter$Age < summary(BDR_metadata_filter$Age)[3],]$age_conc = 2
BDR_metadata_filter[BDR_metadata_filter$Age >= summary(BDR_metadata_filter$Age)[3] & 
                      BDR_metadata_filter$Age < summary(BDR_metadata_filter$Age)[5],]$age_conc = 3
BDR_metadata_filter[BDR_metadata_filter$Age >= summary(BDR_metadata_filter$Age)[5],]$age_conc = 4

pheno$dementia = NA
pheno$diag = NA
#CDR filter
BDR_metadata_filter$CDR[BDR_metadata_filter$CDR == -9] = NA
BDR_metadata_filter$CDR[BDR_metadata_filter$CDR == 18] = NA
BDR_metadata_filter$CDR[BDR_metadata_filter$CDR == 6] = NA
#Dementia assessment
pheno$dementia = NA
pheno$diag = NA

pheno$dementia[pheno$MMSE < 24] = 1
pheno$dementia[pheno$CDR >= 1] = 1
pheno$dementia[pheno$CDR < 1] = 0
pheno$dementia[is.na(pheno$dementia)] = 2 # Anything we can't conclude
#2 should be removed later on (table)
#0   1   2 
#157 340  87 
#Diagnosis
pheno$diag[pheno$CERAD <= 1 & pheno$Braak <= 2 & pheno$dementia == 0] = "Control" 
pheno$diag[pheno$CERAD == 0 & pheno$Braak <= 3 & pheno$dementia == 0] = "Control" 
pheno$diag[pheno$CERAD >= 1 & pheno$Braak >= 3 & pheno$dementia == 0] = "AsymAD"
pheno$diag[pheno$CERAD >= 2 & pheno$Braak >= 3 & pheno$dementia == 1] = "AD" 
pheno$diag[is.na(pheno$diag)] = "Unassigned"

#AD     AsymAD    Control Unassigned 
#217         35        104        228

##### Curating dataset #####
BDR_ok = BDR_metadata_filter[BDR_metadata_filter$diag != "Unassigned",]

# Clinical AD
#We searched for clinical diagnoses of AD made by doctors directly on unassigned samples to get more samples
BDR_unassigned = BDR_metadata_filter[BDR_metadata_filter$diag == "Unassigned",]
BDR_unassigned[grep("Alzheimer", BDR_unassigned$clinical),]$clinical
BDR_unassigned[grep("Alzheimer", BDR_unassigned$clinical),]$diag = "AD"
BDR_ok = rbind(BDR_ok, BDR_unassigned[BDR_unassigned$diag != "Unassigned",])
BDR_unassigned = BDR_unassigned[BDR_unassigned$diag == "Unassigned",]
# We move from 217 AD to 314 AD

# Clinical dementia
BDR_unassigned[grep("dementia", BDR_unassigned$clinical),]$clinical
BDR_unassigned[grep("dementia", BDR_unassigned$clinical),]$dementia = 1

BDR_unassigned$diag[BDR_unassigned$Cerad <= 1 & BDR_unassigned$Braak_tangle <= 2 & BDR_unassigned$dementia == 0] = "Control" 
BDR_unassigned$diag[BDR_unassigned$Cerad == 0 & BDR_unassigned$Braak_tangle <= 3 & BDR_unassigned$dementia == 0] = "Control" 
BDR_unassigned$diag[BDR_unassigned$Cerad >= 1 & BDR_unassigned$Braak_tangle >= 3 & BDR_unassigned$dementia == 0] = "AsymAD"
BDR_unassigned$diag[BDR_unassigned$Cerad >= 2 & BDR_unassigned$Braak_tangle >= 3 & BDR_unassigned$dementia == 1] = "AD" 
table(BDR_unassigned$diag)
BDR_ok = rbind(BDR_ok, BDR_unassigned[BDR_unassigned$diag != "Unassigned",])
BDR_unassigned = BDR_unassigned[BDR_unassigned$diag == "Unassigned",]
#313 AD to 323 AD

# Braak tangle and Cerad NA's can't be fixed - no data
# No abnormalities = controls
BDR_unassigned[grep("No abnormality detected", BDR_unassigned$clinical),]$clinical
BDR_unassigned[grep("No abnormality detected", BDR_unassigned$clinical),]$diag = "Control"
BDR_ok = rbind(BDR_ok, BDR_unassigned[BDR_unassigned$diag != "Unassigned",])
BDR_unassigned = BDR_unassigned[BDR_unassigned$diag == "Unassigned",]
#104 Controls to 144 controls



# Rest
BDR_unassigned$clinical
BDR_unassigned$clinical_plus
BDR_unassigned$pathology
BDR_unassigned$Braak_LB

# Final = We exclude BDR_unassigned 83 remaining samples because not enough information OR other types of Dementia
# Final = We exclude the 35 AsymAD samples and age <65
#313 AD and 140 Controls

BDR_ok_age = BDR_ok[BDR_ok$Age >= 65,]
table(BDR_ok_age$diag)
BDR_ok_ADCtrl = BDR_ok_age[BDR_ok_age$diag != "AsymAD",]

table(BDR_ok_ADCtrl$diag)
table(BDR_ok_ADCtrl$dementia)

table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$MMSE)
table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$CDR)
table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$Braak_tangle)
table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD",]$Cerad)

table(BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "AD" & BDR_ok_ADCtrl$CDR < 1,]$MMSE)



##### Curating Controls #####
BDR_ok_ADCtrl = BDR_ok_ADCtrl[order(BDR_ok_ADCtrl$BBNId),]
BDR_ok_ADCtrl$BBNId

BDR_ok_ADCtrl$Thal_number = NA
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == "Thal phase 0 "] = 0
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == "Thal phase 1 "] = 1
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == "Thal phase 2 "] = 2
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == "Thal phase 3 "] = 3
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == "Thal phase 4 "] = 4
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == "Thal phase 5 "] = 5
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == "Thal phase not assessed "] = NA
BDR_ok_ADCtrl$Thal_number[BDR_ok_ADCtrl$Thal.stage == ""] = NA
table(BDR_ok_ADCtrl$Thal.stage)
table(BDR_ok_ADCtrl$Thal_number)

BDR_Ctrl = BDR_ok_ADCtrl[BDR_ok_ADCtrl$diag == "Control",]
BDR_Ctrl$new_diag = NA
BDR_Ctrl$Thal.stage


#Braak stage 0-2
table(BDR_Ctrl$BraakTangle_numeric)
BDR_Ctrl_BTnonOK = BDR_Ctrl[which(BDR_Ctrl$Braak_tangle >= 4),]

# Argyrophilic
grep("+Ndgen.NonADTau.AGD", BDR_Ctrl$pathology, fixed = TRUE)
BDR_Ctrl_Argnonok = BDR_Ctrl[grep("+Ndgen.NonADTau.AGD", BDR_Ctrl$pathology, fixed = TRUE),]
BDR_Ctrl_Argnonok$pathology  

# Thal and Cerad 0
table(BDR_Ctrl$Cerad, BDR_Ctrl$Thal_number)
BDR_Ctrl_amynonok = BDR_Ctrl[which(BDR_Ctrl$Cerad > 0 | BDR_Ctrl$Thal_number > 0),]

# TDP43
grep("+Ndgen.TDP.TDPaccomp", BDR_Ctrl$pathology, fixed = TRUE)
BDR_Ctrl_TDP43 = BDR_Ctrl[grep("+Ndgen.TDP.TDPaccomp", BDR_Ctrl$pathology, fixed = TRUE),]
BDR_Ctrl_TDP43$pathology

# Alpha synuclein - no info ??
grep("a-", BDR_Ctrl$pathology, fixed = TRUE)
BDR_Ctrl_alpha = BDR_Ctrl[grep("a-", BDR_Ctrl$pathology, fixed = TRUE),]
BDR_Ctrl_alpha$pathology

# Cerebrovascular

grep("vascular", BDR_Ctrl$clinical)
grep("vascular", BDR_Ctrl$clinical_plus)
grep("vascular", BDR_Ctrl$pathology)

final_cereb = c(grep("vascular", BDR_Ctrl$clinical, fixed = TRUE),
                grep("vascular", BDR_Ctrl$clinical_plus, fixed = TRUE),
                grep("+CVD.Vasc.Athero.Sevathero", BDR_Ctrl$pathology, fixed = TRUE))
BDR_Ctrl_cereb = BDR_Ctrl[unique(final_cereb),]
BDR_Ctrl_cereb$clinical
BDR_Ctrl_cereb$clinical_plus
BDR_Ctrl_cereb$pathology

# Leukemia - None ?
grep("LEUKEMIA", BDR_Ctrl$Cause_of_Death)

# Neuropathologic correlate - None ?
grep("psychia", BDR_Ctrl$clinical_plus)


# Final assessment
final_ctrl_nonok = c(grep("vascular", BDR_Ctrl$clinical, fixed = TRUE),
                     grep("vascular", BDR_Ctrl$clinical_plus, fixed = TRUE),
                     grep("+CVD.Vasc.Athero.Sevathero", BDR_Ctrl$pathology, fixed = TRUE),
                     grep("+Ndgen.TDP.TDPaccomp", BDR_Ctrl$pathology, fixed = TRUE),
                     which(BDR_Ctrl$Cerad > 1 | BDR_Ctrl$Thal_number > 2),
                     grep("+Ndgen.NonADTau.AGD", BDR_Ctrl$pathology, fixed = TRUE),
                     which(BDR_Ctrl$Braak_tangle >= 3))
final_ctrl_nonok = unique(final_ctrl_nonok)


BDR_remove = BDR_Ctrl[-final_ctrl_nonok,]
BDR_remove$Basename

}
