train=read.table("train.snpinfo",stringsAsFactors=F,header=T)
ss=read.table("preSummStat",stringsAsFactors=F,header=T)

newFreq=rep(0,nrow(ss))
newFreq[ss$ref == train$alleleA]=ss$reffreq[ss$ref == train$alleleA]
newFreq[ss$ref != train$alleleA]=1-ss$reffreq[ss$ref != train$alleleA]
ss$ref=train$alleleA
ss$alt=train$alleleB
ss$reffreq=newFreq
write.table(ss,"newSummStat",quote=F,sep='\t',col.names=T,row.names=F)
