unmelt <- function() {
  setwd("C:/Users/P70077107/Desktop/PhD")
  resmq = read.delim("mqres/txt/evidence.txt")
  unique(resmq$Raw.file)
  
  resmq_filt = data.frame(resmq$Proteins,
                          resmq$Raw.file,
                          resmq$MS.MS.count)
  
  resmq_table = dcast(resmq_filt, resmq.Proteins ~ resmq.Raw.file,value.var = "resmq.MS.MS.count", fun.aggregate = sum)
}