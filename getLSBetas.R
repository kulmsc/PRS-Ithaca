pseudo=readRDS("lassosum.pseudovalidate.rds")
splitval=readRDS("lassosum.splitvalidate.rds")
val=readRDS("lassosum.validate.rds")
pipe=readRDS("lassosum.lassosum.pipeline.rds")
ss=read.table("preSummStat",stringsAsFactors=F,header=T)

ssout=ss[pipe$sumstats$order,]
ssout$beta=val$best.beta
write.table(ssout,"ss1",row.names=F,col.names=F,quote=F,sep='\t')

ssout$beta=splitval$best.beta[,2]
write.table(ssout,"ss2",row.names=F,col.names=F,quote=F,sep='\t')

ssout$beta=pseudo$best.beta
write.table(ssout,"ss3",row.names=F,col.names=F,quote=F,sep='\t')
