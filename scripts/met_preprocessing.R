met_preprocessing <- function(){
  
  library(wateRmelon)
  load("./mSet454.BDR.Rdata")
 
  outliers <- outlyx(mSet, plot=FALSE)
  
  print(outliers)
  
  mSet_outliers = mSet[,!outliers$out]
  
  bsc <- bscon(mSet_outliers)
  
  mSet_pfilter <- pfilter(mSet_outliers)
  
  estimateSex(betas(mSet_pfilter), do_plot=FALSE)
  
  estimateCellCounts.wmln(mSet_pfilter, referencePlatform = "IlluminaHumanMethylationEPIC",
                          compositeCellType = "DLPFC", cellTypes = c("NeuN_neg","NeuN_pos"))

  mSet_BMIQ <- BMIQ(mSet_pfilter)
  
  qu <- qual(betas(melon), betas(das))
  
  dmrse_row(mSet_pfilter)
  
  genki(mSet_pfilter)
  
  seabi(mSet_pfilter, sex=pData(mSet_pfilter)$sex, X=fData(mSet_pfilter)$CHR=='X')
  
  mSet_betas <- betas(mSet_pfilter)
  mSet_Mvalues <- beta2m(mSet_betas)
  pwod_bet <- pwod(bet)
  
}