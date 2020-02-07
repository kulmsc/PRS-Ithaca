
ss <- read.table("preSummStat",stringsAsFactors=F)
colnames(ss) <- c("CHR","MarkerName","BP","Allele1","Allele2","FRQ","INFO","BETA","SE","P","SAMPSIZE","FRQ2")
ss$Z <- qnorm(ss$P/2)
write.table(ss,"forImpute",col.names=T,row.names=F,quote=F)
