library("GraBLD",lib.loc="/home/sdk2004/bin")

makeMean <- function(x){
	val=mean(na.omit(x))
	x[is.na(x)]=val
	return(x)
}

args = commandArgs(TRUE)

#read in main data
geno_data=load_geno(paste0("forLD.",args[1],".raw"))

#impute means and remove those with af=1
vecGenoData=apply(geno_data,2,makeMean)
geno_data=matrix(vecGenoData,nrow=nrow(geno_data),ncol=ncol(geno_data))
badCol=apply(geno_data,2,function(x) length(unique(x)))
if(sum(badCol==1)>0){
	geno_data=geno_data[,-which(badCol==1)]
} 


#write down which of the rsids did not make it
con <- file(paste0("forLD.",args[1],".raw"),"r")
firstLine <- readLines(con,n=1)
close(con)
rsids=strsplit(firstLine,' ')[[1]][-(1:6)]
rsids=unlist(lapply(rsids,strsplit,'_'))
rsids=rsids[seq(1,length(rsids),2)]
if(sum(badCol==1)>0){
	badRsids=rsids[which(badCol==1)]
	writeRsids=rsids[-which(badCol==1)]
} else {
	badRsids="rs0"
	writeRsids=rsids
}
write.table(badRsids,paste0("grabldBadRsids.",args[1]),quote=F,row.names=F,col.names=F)
write.table(writeRsids,paste0("grabldAllRsids.",args[1]),quote=F,row.names=F,col.names=F)

#normalize geno_data and write back for GWAS
geno_norm=full_normal_geno(geno_data)
LD_OUT=LDadj(geno_raw=geno_norm)
saveRDS(LD_OUT,file=paste0("grabldLD_OUT.",args[1],".rds"))
