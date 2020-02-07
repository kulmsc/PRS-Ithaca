
args = commandArgs(TRUE)

decoder=read.table("nameNumberDecoder",stringsAsFactors=F,header=F)
files=read.table("fileDecoder",stringsAsFactors=F,header=T)
fileName=strsplit(args[2],".",fixed=T)[[1]][1]
code=files[files$author==fileName,3]
allDisease=c()
for(code in strsplit(code,",")[[1]]){
	tempDecode=decoder[decoder[,2]==code,1]
	finalDecode=gsub("/",".",tempDecode)
	allDisease=c(allDisease,finalDecode)
}


allScores=read.table("allScoresCompiled.train",header=T,stringsAsFactors=F)
covars=read.table("covars",header=T,row.names=1,stringsAsFactors=F)
fam=read.table(paste0("../new.",args[1],".fam"),stringsAsFactors=F)
fam[covars$sex=="M",5]=1
fam[covars$sex=="F",5]=2

manyScore=allScores[,colnames(allScores) %in% allDisease]
if(length(allDisease)>1){
	score=apply(manyScore,1,sum)
	score[score>1]=1
} else {
	score=manyScore
}
fam[,6]=score+1

fullCovar=cbind(fam[,1:2],covars[,c(2,3,4,5,7)])

write.table(fam,"gwas.fam",row.names=F,col.names=F,quote=F,sep='\t')
write.table(fullCovar,"fullCovar.cov",row.names=F,col.names=F,quote=F,sep='\t')


#########################################################################

#geno=read.table("test.raw",header=T,stringsAsFactors=F)

#geno=geno[,7:ncol(geno)]
#score=allScores[,colnames(allScores)==disease]
#score=factor(score)
#df=cbind(score,covars)

#pval=vector(length=ncol(geno))
#beta=vector(length=ncol(geno))

#gwas <- function(snp,df){
#	res=glm(score ~ . + geno[,1],data=df,family="binomial")
#	beta=coef(summary(res))[9,1]
#	pval=coef(summary(res))[9,4]
#	return(list(beta,pval))
#}

#ans=apply(geno,2,gwas,df)

#beta[1]=coef(summary(res))[9,1]
#pval[1]=coef(summary(res))[9,4]
