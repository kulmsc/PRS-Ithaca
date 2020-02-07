args = commandArgs(TRUE)

decoder=read.table("nameNumberDecoder",stringsAsFactors=F,header=F)
files=read.table("fileDecoder",stringsAsFactors=F,header=T)
fileName=strsplit(args[1],".",fixed=T)[[1]][1]
code=files[files$author==fileName,3]
allDisease=c()
for(code in strsplit(code,",")[[1]]){
	tempDecode=decoder[decoder[,2]==code,1]
	finalDecode=gsub("/",".",tempDecode)
	allDisease=c(allDisease,finalDecode)
}


allScores=read.table("allScoresCompiled.train",header=T,stringsAsFactors=F)

manyScore=allScores[,colnames(allScores) %in% allDisease]
if(length(allDisease)>1){
	score=apply(manyScore,1,sum)
	score[score>1]=1
} else {
	score=manyScore
}

write.table(score,"diseaseStatus",row.names=F,col.names=F,quote=F,sep='\t')
