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
coefs=as.matrix(coefs[,8])
coefs=coefs/herit
totalAnnot=as.matrix(totalAnnot)
ans=totalAnnot %*% coefs
ans=cbind(rep(preSummStat[1,1],nrow(ans)),snps,ans)
write.table(ans,"funcEnrich",row.names=F,col.names=F,quote=F,sep=' ')
