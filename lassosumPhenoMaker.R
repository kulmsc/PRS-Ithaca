args = commandArgs(TRUE)


authorName=strsplit(args[2],".",fixed=T)[[1]][1]
allResponses <- read.table("~/prsDatabase/allResponses",stringsAsFactors=F,header=T)
score <- allResponses[,colnames(allResponses)==authorName]

print(length(score))
famFile=read.table(paste0("../new.",args[1],".fam"),stringsAsFactors=F)
famFile[,6]=score+1
if(sum(famFile[,6]==2)<2000){
	print("small")
	keepFam=c(sample(famFile[famFile[,6]==1,1],2000),famFile[famFile[,6]==2,1])
} else {
	print("big")
	keepFam=c(sample(famFile[famFile[,6]==1,1],2000), sample(famFile[famFile[,6]==2,1],2000))
}

famFile=famFile[famFile[,1] %in% keepFam,]

write.table(famFile,"our.fam",sep=' ',row.names=F,col.names=F,quote=F)
