updbeta=read.table("updBeta",stringsAsFactors=F)
binbeta=read.table("binBeta",stringsAsFactors=F)

for(bin in 1:101){
	adj=binbeta[bin,2]
	if(sum(updbeta[,3]==bin) > 0){
		updbeta[updbeta[,3]==bin,2]=updbeta[updbeta[,3]==bin,2]*adj
	}
}

write.table(updbeta,"finalBeta",quote=F,sep='\t',row.names=F,col.names=F)
