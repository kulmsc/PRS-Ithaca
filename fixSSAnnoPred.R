ss=read.table("preSummStat",stringsAsFactors=F)
newSS=ss[,c(1,2,4,5,3,8,10)]
newSS[,1]=paste0("chr",newSS[,1])
newSS[,6]=exp(newSS[,6])
colnames(newSS)=c("hg19chrc","snpid","a1","a2","bp","or","p")
write.table(newSS,"annoSummStat",row.names=F,col.names=T,quote=F,sep=' ')
