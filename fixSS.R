valid=read.table("testDosInfo",stringsAsFactors=F,header=T)
train=read.table("trainDosInfo",stringsAsFactors=F,header=T)
ss=read.table("preSummStat",stringsAsFactors=F,header=T)
ss=ss[order(ss[,2]),]

newFreq=rep(0,nrow(ss))
newFreq[ss$ref == train$alleleA]=ss$reffreq[ss$ref == train$alleleA]
newFreq[ss$ref != train$alleleA]=1-ss$reffreq[ss$ref != train$alleleA]
ss$ref=train$alleleA
ss$alt=train$alleleB
ss$reffreq=newFreq
write.table(ss,"newSummStat",quote=F,sep='\t',col.names=T,row.names=F)
