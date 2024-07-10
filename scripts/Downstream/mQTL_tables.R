mQTL_tables <- function(){
  # Extract mQTLs for the red/blue significant sets to test for genetic colocalization
  table_mqtls_EPICredC = data.frame()
  table_mqtls_450kredC = data.frame()
  table_mqtls_EPICblueC = data.frame()
  table_mqtls_450kblueC = data.frame()
  
  chr = dir("./mQTLdb/",full.names=TRUE)
  
  # Get mQTLs per chromosome
  for(chrfile in chr){
    tmp_table = read.csv(chrfile)
    
    tmp_table_EPICredC = tmp_table[tmp_table$CpG %in% rownames(sig_methyl_EPIC_redC),]
    tmp_table_EPICredC = tmp_table_EPICredC[tmp_table_EPICredC$p < 1e-5,]
    
    table_mqtls_EPICredC = rbind(table_mqtls_EPICredC, tmp_table_EPICredC)
    
    tmp_table_450kredC = tmp_table[tmp_table$CpG %in% rownames(sig_methyl_450k_redC),]
    tmp_table_450kredC = tmp_table_450kredC[tmp_table_450kredC$p < 1e-5,]
    
    table_mqtls_450kredC = rbind(table_mqtls_450kredC, tmp_table_450kredC)
    
    tmp_table_EPICblueC = tmp_table[tmp_table$CpG %in% rownames(sig_methyl_EPIC_blueC),]
    tmp_table_EPICblueC = tmp_table_EPICblueC[tmp_table_EPICblueC$p < 1e-5,]
    
    table_mqtls_EPICblueC = rbind(table_mqtls_EPICblueC, tmp_table_EPICblueC)
    
    tmp_table_450kblueC = tmp_table[tmp_table$CpG %in% rownames(sig_methyl_450k_blueC),]
    tmp_table_450kblueC = tmp_table_450kblueC[tmp_table_450kblueC$p < 1e-5,]
    
    table_mqtls_450kblueC = rbind(table_mqtls_450kblueC, tmp_table_450kblueC)

    
  }
  table_mqtls_EPICredC$N = nrow(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red", "Control"),]) +
    nrow(pheno_PITT[pheno_PITT$subtype %in% c("red", "Control"),])
  table_mqtls_450kredC$N = nrow(pheno_UKBBN[pheno_UKBBN$subtype %in% c("red", "Control"),]) +
    nrow(pheno_PITT[pheno_PITT$subtype %in% c("red", "Control"),]) +
    nrow(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("red", "Control"),])
  table_mqtls_EPICblueC$N = nrow(phenoUKBBN[pheno_UKBBN$subtype %in% c("blue", "Control"),]) +
    nrow(pheno_PITT[pheno_PITT$subtype %in% c("blue", "Control"),])
  table_mqtls_450kblueC$N = nrow(pheno_UKBBN[pheno_UKBBN$subtype %in% c("blue", "Control"),]) +
    nrow(pheno_PITT[pheno_PITT$subtype %in% c("blue", "Control"),]) +
    nrow(pheno_ROSMAP[pheno_ROSMAP$subtype %in% c("blue", "Control"),])
  
  table_MAF = data.frame()
  
  # Add MAF values per mQTL
  MAF = dir("./maf/",full.names=TRUE)
  
  table_mqtls_EPICredC$MAF = NA
  table_mqtls_450kredC$MAF = NA
  table_mqtls_EPICblueC$MAF = NA
  table_mqtls_450kblueC$MAF = NA
  
  for(MAFfile in MAF){
    tmp_MAF = read.table(MAFfile)
    table_MAF = rbind(table_MAF,tmp_MAF)
  }
  
  colnames(table_MAF) = c("RSID","MAF")
  
  #Add MAF values and save tables for both EPIC/450K
  colnames(table_mqtls_EPICredC) = toupper(colnames(table_mqtls_EPICredC))
  MAF_table = table_MAF[table_MAF$RSID %in% table_mqtls_EPICredC$RSID,]
  table_mqtls_EPICredC = table_mqtls_EPICredC[table_mqtls_EPICredC$RSID %in% MAF_table$RSID,]
  table_mqtls_EPICredC = table_mqtls_EPICredC[order(table_mqtls_EPICredC$RSID),]
  MAF_table = MAF_table[order(MAF_table$RSID),]
  
  table_mqtls_EPICredC$MAF = MAF_table[match(table_mqtls_EPICredC$RSID, MAF_table$RSID),]$MAF
  write.csv(table_mqtls_EPICredC, file = "./table_mqtls_EPICredC.csv")
  
  colnames(table_mqtls_450kredC) = toupper(colnames(table_mqtls_450kredC))
  MAF_table = table_MAF[table_MAF$RSID %in% table_mqtls_450kredC$RSID,]
  table_mqtls_450kredC = table_mqtls_450kredC[table_mqtls_450kredC$RSID %in% MAF_table$RSID,]
  table_mqtls_450kredC = table_mqtls_450kredC[order(table_mqtls_450kredC$RSID),]
  MAF_table = MAF_table[order(MAF_table$RSID),]
  
  table_mqtls_450kredC$MAF = MAF_table[match(table_mqtls_450kredC$RSID, MAF_table$RSID),]$MAF
  write.csv(table_mqtls_450kredC, file = "./table_mqtls_450kredC.csv")
  
  
  colnames(table_mqtls_EPICblueC) = toupper(colnames(table_mqtls_EPICblueC))
  MAF_table = table_MAF[table_MAF$RSID %in% table_mqtls_EPICblueC$RSID,]
  table_mqtls_EPICblueC = table_mqtls_EPICblueC[table_mqtls_EPICblueC$RSID %in% MAF_table$RSID,]
  table_mqtls_EPICblueC = table_mqtls_EPICblueC[order(table_mqtls_EPICblueC$RSID),]
  MAF_table = MAF_table[order(MAF_table$RSID),]
  
  table_mqtls_EPICblueC$MAF = MAF_table[match(table_mqtls_EPICblueC$RSID, MAF_table$RSID),]$MAF
  write.csv(table_mqtls_EPICblueC, file = "./table_mqtls_EPICblueC.csv")
  
  
  colnames(table_mqtls_450kblueC) = toupper(colnames(table_mqtls_450kblueC))
  MAF_table = table_MAF[table_MAF$RSID %in% table_mqtls_450kblueC$RSID,]
  table_mqtls_450kblueC = table_mqtls_450kblueC[table_mqtls_450kblueC$RSID %in% MAF_table$RSID,]
  table_mqtls_450kblueC = table_mqtls_450kblueC[order(table_mqtls_450kblueC$RSID),]
  MAF_table = MAF_table[order(MAF_table$RSID),]
  
  table_mqtls_450kblueC$MAF = MAF_table[match(table_mqtls_450kblueC$RSID, MAF_table$RSID),]$MAF
  write.csv(table_mqtls_450kblueC, file = "./table_mqtls_450kblueC.csv")
  
  
}

