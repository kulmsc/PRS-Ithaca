library(lassosum)
source("lassosum.splitVal.R")

args = commandArgs(trailingOnly=TRUE) #sampleSize Chr file

print(args[1])
ss = read.table(args[4],stringsAsFactors=F)
cor = p2cor(p=ss[,10],n=as.numeric(args[1]),sign=ss[,8])
print("done cor")

#3 Fold CV ##################################################################################
decoder=read.table("~/prsDatabase/nameNumberDecoder",stringsAsFactors=F,header=F)
files=read.table("~/prsDatabase/fileDecoder",stringsAsFactors=F,header=T)
allScores=read.table("allScoresCompiled.train",header=T,stringsAsFactors=F)
famFile=read.table(paste0("../new.",args[2],".fam"))

#This chunk of code will make the pheno or score file
fileName=strsplit(args[3],".",fixed=T)[[1]][1]
code=files[files$author==fileName,3]
allDisease=c()
for(code in strsplit(code,",")[[1]]){
        tempDecode=decoder[decoder[,2]==code,1]
        finalDecode=gsub("/",".",tempDecode)
        allDisease=c(allDisease,finalDecode)
}
manyScore=allScores[,colnames(allScores) %in% allDisease]
if(length(allDisease)>1){
        score=apply(manyScore,1,sum)
        score[score>1]=1
} else {
        score=manyScore
}

groupCuts=seq(1,nrow(famFile),length.out=4)
group1=famFile[groupCuts[1]:groupCuts[2],1]
group2=famFile[groupCuts[2]:groupCuts[3],1]
group3=famFile[groupCuts[3]:groupCuts[4],1]
testGroups=list(group1,group2,group3)
allLambda=vector(length=3)
allS=vector(length=3)
for(i in 1:3){
	testPheno=score[famFile[,1] %in% testGroups[[i]]]
	print("begin pipeline")
	out = lassosum.pipeline(cor=cor,LDblocks="EUR.hg19",
		chr=ss[,1],pos=ss[,3],A1=ss[,4],A2=ss[,5],max.ref.bfile.n=1000000,
		ref.bfile=paste0("../new.",args[2]),test.bfile=paste0("../new.",args[2]),
		keep.ref=!(famFile[,1] %in% testGroups[[i]]),keep.test=famFile[,1] %in% testGroups[[i]])
	print("done pipeline")
	v <- validate(out, pheno=testPheno)
	allLambda[i]=v$best.lambda
	allS[i]=v$best.s
}
out = lassosum.pipeline(cor=cor,LDblocks="EUR.hg19",max.ref.bfile.n=1000000,
                chr=ss[,1],pos=ss[,3],A1=ss[,4],A2=ss[,5],
                ref.bfile=paste0("../new.",args[2]),
                lambda=mean(allLambda),s=mean(allS))
thisSS=ss[ss[,3] %in% out$sumstats$pos,]
thisSS[,8]=out$beta[[1]][,1]
write.table(thisSS,"summStat.temp.1",quote=F,sep='\t',col.names=F,row.names=F)


#split validate #########################################################################################
out.reuse = lassosum.pipeline(cor=cor,LDblocks="EUR.hg19",max.ref.bfile.n=1000000,
                chr=ss[,1],pos=ss[,3],A1=ss[,4],A2=ss[,5],
                ref.bfile=paste0("../new.",args[2]),test.bfile=paste0("../new.",args[2]))
save.image("test.RData")
sv <- splitvalidate.mine(out.reuse,pheno=score)
out = lassosum.pipeline(cor=cor,LDblocks="EUR.hg19",max.ref.bfile.n=1000000,
                chr=ss[,1],pos=ss[,3],A1=ss[,4],A2=ss[,5],
                ref.bfile=paste0("../new.",args[2]),lambda=sv$best.lambda,s=sv$best.s)
thisSS=ss[ss[,3] %in% out$sumstats$pos,]
thisSS[,8]=out$beta[[1]][,1]
write.table(thisSS,"summStat.temp.2",quote=F,sep='\t',col.names=F,row.names=F)

#pseudo validate ##########################################################################################
v <- pseudovalidate(out.reuse)
out = lassosum.pipeline(cor=cor,LDblocks="EUR.hg19",
                chr=ss[,1],pos=ss[,3],A1=ss[,4],A2=ss[,5],max.ref.bfile.n=1000000,
                ref.bfile=paste0("../new.",args[2]),lambda=v$best.lambda,s=v$best.s)
thisSS=ss[ss[,3] %in% out$sumstats$pos,]
thisSS[,8]=out$beta[[1]][,1]
write.table(thisSS,"summStat.temp.3",quote=F,sep='\t',col.names=F,row.names=F)
