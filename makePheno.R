
allCodes=c(1081,1222,1289,1408,1471)
pheno=read.table("shortPheno.csv",sep=',',stringsAsFactors=F,header=T,row.names=1)
pheno=pheno[!(is.na(pheno$X22009.0.1)),8:36]

output=cbind(rownames(pheno),rownames(pheno))
for(code in allCodes){
	response=apply(pheno,1,function(x) ifelse(code %in% x,1,0))
	response[response==1]=2
	response[response==0]=1
	output=cbind(output,response)
}

write.table(output,"phenoCoding",row.names=F,col.names=F,quote=F,sep='\t')
