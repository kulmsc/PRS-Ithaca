library("FNN")

args=commandArgs(trailingOnly=TRUE)
herit=as.numeric(args[1])


totalAnnot=readRDS("totalAnnot.RDS")
preSummStat=read.table("preSummStat",stringsAsFactors=F)
coefs=read.table("coefs",stringsAsFactors=F,header=T) #given in the github instructions

preSummStat=preSummStat[preSummStat[,4]!=preSummStat[,5],]

#make sure the snps match between annotation and presummstat
totalAnnot=totalAnnot[totalAnnot$SNP %in% preSummStat[,2],]
preSummStat=preSummStat[preSummStat[,2] %in% totalAnnot$SNP,]
oldAnnot=data.frame(totalAnnot)

#order to match them up
totalAnnot=totalAnnot[order(totalAnnot$SNP),]
preSummStat=preSummStat[order(preSummStat[,2]),]

#subset the columns correctly
snps=totalAnnot[,3]
totalAnnot=totalAnnot[,5:ncol(totalAnnot)]

#form the tau vector and then multiply by annot to get heritability estimates for every snp
coefs=coefs[1:ncol(totalAnnot),]
coefs=as.matrix(coefs[,8])
coefs=coefs/herit
totalAnnot=as.matrix(totalAnnot)
ans=totalAnnot %*% coefs
ans=cbind(snps,ans)
write.table(ans,"funcEnrich",row.names=F,col.names=F,quote=F,sep=' ')

#now fortmat the summstat as they want it
ss=preSummStat[,c(1:5,10,8,10)]
ss[,8]=abs(qnorm(ss[,8]/2))*sign(ss[,7])
colnames(ss)=c("CHR","SNP","BP","A1","A2","P","BETA","Z")
write.table(ss,"functSummStat",col.names=T,row.names=F,quote=F,sep='\t')

#format again but making chi square statistics
#pTOs <- read.table("pToStat",stringsAsFactors=T,header=T)
#pTOs <- pTOs[-which(is.na(pTOs[,2])),]
#ss <- read.table("preSummStat",stringsAsFactors=F,header=F)
#ss=ss[ss[,2] %in% oldAnnot$SNP,]
#ss=ss[order(ss[,2]),]
#newP <- as.matrix(ss[,10])
#oldP <- as.matrix(pTOs[,1])
#oldStat <- as.matrix(pTOs[,2])
#colnames(newP)<-"pval"
#colnames(oldP)<-"pval"
#colnames(oldStat)<-"stat"
#newStat <- knn.reg(oldP,newP,oldStat,3)
#ss=preSummStat[,c(1:5,10,8,10)]
#ss[,8]=newStat$pred
#colnames(ss)=c("CHR","SNP","BP","A1","A2","P","BETA","CHISQ")
#write.table(ss,"summStat_2",col.names=T,row.names=F,quote=F,sep='\t')
