allNames=read.table("allNames2",stringsAsFactors=F)
#print(allNames)

for(fName in allNames[,1]){
	print(fName)
	allNormal=list()
	for(i in 1:22){
		print(paste0("../chr",as.character(i),"/",fName,".ss.",as.character(i),".gz"))
		gFile=read.table(paste0("../chr",i,"/",fName,".ss.",i,".gz"),stringsAsFactors=F)
		allNormal[[i]]=gFile
	}
	print("reading in each ldpred")
	for(i in 1:3){
		system(paste("zcat",paste0(fName,".",i,".final.set.ldpred.gz"),"| grep chrom > readFrom"))
		print(i)
		fixFile=read.table("readFrom",stringsAsFactors=F,header=F)
		blankFix=c()
		for(j in 1:22){
			gFile=allNormal[[j]]
			subFixFile=fixFile[fixFile[,1]==paste0("chrom_",j),]
			subGFile=gFile[gFile[,2] %in% subFixFile[,3],]
			
			subGFile=subGFile[order(subGFile[,2]),]
			subFixFile=subFixFile[order(subFixFile[,3]),]
			subGFile[,8]=subFixFile[,7]
			blankFix=rbind(blankFix,subGFile)
		}
		write.table(blankFix,paste0(fName,".",i,".new.stuff"),row.names=F,col.names=F,sep='\t',quote=F)
	}
}
