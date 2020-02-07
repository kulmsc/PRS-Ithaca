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
manyScore=allScores[,colnames(allScores) %in% allDisease]
if(length(allDisease)>1){
        score=apply(manyScore,1,sum)
        score[score>1]=1
} else {
        score=manyScore
}

samp=read.table("phase.eid",stringsAsFactors=F)
allCase=samp[score==1,]
allControl=samp[score==0,]
trainCase=allCase[1:round(length(allCase)/2)]
validCase=allCase[(round(length(allCase)/2)+1):length(allCase)]
trainControl=allControl[1:round(length(allCase)/2)]
validControl=allControl[(round(length(allCase)/2)+1):length(allCase)]

comboFam=c(trainCase,trainControl)
trainFam=cbind(comboFam,comboFam,rep(0,length(comboFam)),rep(0,length(comboFam)),rep(0,length(comboFam)),rep(-9,length(comboFam)))
comboFam=c(validCase,validControl)
validFam=cbind(comboFam,comboFam,rep(0,length(comboFam)),rep(0,length(comboFam)),rep(0,length(comboFam)),rep(-9,length(comboFam)))

comboFam=c(trainCase,trainControl)
trainPheno=cbind(comboFam,comboFam,c(rep(1,length(trainCase)),rep(0,length(trainControl))))
comboFam=c(validCase,validControl)
validPheno=cbind(comboFam,comboFam,c(rep(1,length(validCase)),rep(0,length(validControl))))

colnames(trainPheno)=c("FID","IID","Outcome")
colnames(validPheno)=c("FID","IID","Outcome")

trainFam=rbind(trainFam,validFam)
trainPheno=rbind(trainPheno,validPheno)

write.table(trainFam,"train.fam",quote=F,sep=' ',col.names=F,row.names=F)
#write.table(validFam,"valid.fam",quote=F,sep=' ',col.names=F,row.names=F)
write.table(trainPheno,"train.phen",quote=F,sep='\t',col.names=T,row.names=F)
#write.table(validPheno,"valid.phen",quote=F,sep='\t',col.names=T,row.names=F)
