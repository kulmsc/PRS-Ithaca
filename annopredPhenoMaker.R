args = commandArgs(TRUE)

authorName=strsplit(args[2],".",fixed=T)[[1]][1]
allResponses <- read.table("~/prsDatabase/allResponses",stringsAsFactors=F,header=T)
score <- allResponses[,colnames(allResponses)==authorName]

#decoder=read.table("nameNumberDecoder",stringsAsFactors=F,header=F)
#files=read.table("fileDecoder",stringsAsFactors=F,header=T)
#fileName=strsplit(args[2],".",fixed=T)[[1]][1]
#code=files[files$author==fileName,3]
#allDisease=c()
#for(code in strsplit(code,",")[[1]]){
#        tempDecode=decoder[decoder[,2]==code,1]
#        finalDecode=gsub("/",".",tempDecode)
#        allDisease=c(allDisease,finalDecode)
#}

#allScores=read.table("allScoresCompiled.train",header=T,stringsAsFactors=F)
#manyScore=allScores[,colnames(allScores) %in% allDisease]
#if(length(allDisease)>1){
#        score=apply(manyScore,1,sum)
#        score[score>1]=1
#} else {
#        score=manyScore
#}

print(length(score))
famFile=read.table(paste0("../new.",args[1],".fam"),stringsAsFactors=F)
famFile[,6]=score+1
write.table(famFile,"our.fam",sep=' ',row.names=F,col.names=F,quote=F)
