library("GraBLD",lib.loc="/home/sdk2004/bin")

args = commandArgs(TRUE)

summStat=read.table("preSummStat",stringsAsFactors=F)
#gwasRes=read.table("gwasRes.assoc.logistic",stringsAsFactors=F,header=T)

rawFile=load_geno("forGwas.raw")
normGeno=full_normal_geno(rawFile)
#forGwasFam=read.table("forGwas.fam",stringsAsFactors=F)
#gwasFam=read.table("gwas.fam",stringsAsFactors=F)
#gwasFam=gwasFam[gwasFam$V1 %in% forGwasFam$V1,]
gwasFam=read.table("our.fam",stringsAsFactors=F)

gwasBeta=apply(normGeno,2,function(val) coef(summary(glm(factor(gwasFam$V6) ~ val, family="binomial")))[2,1])
justBeta=cbind(rep(args[1],length(gwasBeta)),gwasBeta)

#gwasRes=gwasRes[complete.cases(gwasRes),]
#gwasRes=gwasRes[!(duplicated(gwasRes$SNP)),]

#summStat=summStat[summStat$V2 %in% gwasRes$SNP,]
fullAnno=summStat[,c(2,8)]
fullAnno[,1]=as.factor(fullAnno[,1])

#justBeta=gwasRes[,c(1,7)]
write.table(justBeta,"temp_univariate_beta.txt",sep='\t',row.names=F,col.names=F,quote=F)
ldout=LDadj(geno_raw=normGeno)

doGrab <- function(iter,inSummStat=summStat,depth=5,bag=0.5,shri=0.001,sigVal=0.00001,reallyDo=TRUE){
	if(reallyDo){
		totalBetas=c()
		for(i in 1:5){
			#sigVal=0.05
			#ans="try-error"
			#while(ans[1]=="try-error"){
				ans=GraB(betas=justBeta,annotations=fullAnno,pos=2,pos_sign=2,abs_effect=2,which.var=2,trait_name="temp",
					sig=sigVal,steps=i,interact_depth=depth,bag_frac=bag,shrink=shri)
			#	sigVal=sigVal+0.05
			#}
			totalBetas=c(totalBetas,ans)
		}

		#ldAdj=readRDS(paste0("../grabldLD_OUT.",args[1],".rds"))
		#rsids=read.table(paste0("../grabldAllRsids.",args[1]),stringsAsFactor=F,header=F)

		#ldAdj=ldAdj[!duplicated(rsids[,1]),]
		#rsids=rsids[!duplicated(rsids[,1]),]

		#ldAdj=ldAdj[rsids %in% summStat[,2]]

		inSummStat[,8]=totalBetas/ldout[,1]
		write.table(inSummStat,paste0("summStat",iter),row.names=F,col.names=F,sep='\t',quote=F)
		inSummStat[,8]=totalBetas
		iter=iter+1
		write.table(inSummStat,paste0("summStat",iter),row.names=F,col.names=F,sep='\t',quote=F)
	} else {
		inSummStat[,8]=inSummStat[,8]/ldout[,1]
		write.table(inSummStat,paste0("summStat",iter),row.names=F,col.names=F,sep='\t',quote=F)
	}
}

#doGrab(1,reallyDo=FALSE)
#try(doGrab(2))

#print("complete 1")
#try(doGrab(3,bag=0.25))

#print("complete 2")
#try(doGrab(4,bag=0.75))

#print("complete 2.5")
#try(doGrab(5,sigVal=0.1))

print("complete 3")
try(doGrab(1,sigVal=0.001))

print("complete 4")
try(doGrab(2,sigVal=0.00000001))

print("complete 5")
try(doGrab(3,depth=3))

print("complete 6")
try(doGrab(4,depth=10))


