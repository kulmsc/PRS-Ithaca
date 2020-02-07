library("FNN")
args <- commandArgs(trailingOnly=TRUE)

ss <- read.table("win.ss",stringsAsFactors=F,header=F)

#pTOs <- read.table("pToStat",stringsAsFactors=T,header=T)
#pTOs <- pTOs[-which(is.na(pTOs[,2])),]
#newP <- as.matrix(ss[,10])
#oldP <- as.matrix(pTOs[,1])
#oldStat <- as.matrix(pTOs[,2])
#colnames(newP)<-"pval"
#colnames(oldP)<-"pval"
#colnames(oldStat)<-"stat"
#newStat <- knn.reg(oldP,newP,oldStat,3)

#ss$se=abs(ss[,8])/abs(newStat$pred)
ss$se=ss[,9]
out <- data.frame(ss[,c(2,4,5,6,8,11,10)])
colnames(out)=c("SNP","A1","A2","freq","b","se","p")
out$n=rep(args[1],nrow(out))
write.table(out,"ss.ma",col.names=T,row.names=F,quote=F,sep=" ")

