#plink --bfile new --exclude dupIDs --r2 --ld-window-kb 500 --ld-window 1000 --ld-window-r2 0.5 --ld-snp-list totalRsids

#NOTE!!!! THIS IS NOT PERFECT, WILL ONLY GIVE ABOUT 99% OF THE SNPS THAT THE FULL CLUMP WOULD

ld=read.table("plink.ld",header=T,stringsAsFactors=F)
totRsids=read.table("totalRsids",stringsAsFactors=F)
report=read.table("headTotalSummStat",stringsAsFactors=F)

ldSpec=ld[(ld$SNP_A %in% totRsids$V1 & ld$SNP_B %in% totRsids$V1),]
ldSpec=ldSpec[ldSpec$SNP_A != ldSpec$SNP_B,]

specRsids=unique(c(ldSpec$SNP_A,ldSpec$SNP_B))
specPval=report[report$V2 %in% specRsids,10]
specRsids=report[report$V2 %in% specRsids,2]
specRsids=specRsids[order(specPval)]
specPval=specPval[order(specPval)]


removeSnp=c()
i=1
while(i<length(specRsids)){
	lowID=specRsids[i]
	nowRemove=unique(c(ldSpec[ldSpec$SNP_A == lowID,6],ldSpec[ldSpec$SNP_B == lowID,3]))
	specRsids=specRsids[!(specRsids %in% nowRemove)]
	ldSpec=ldSpec[ldSpec$SNP_A != lowID | ldSpec$SNP_B != lowID,]
	removeSnp=c(removeSnp,nowRemove)
	i=i+1
}

doneRsids=totRsids[!(totRsids %in% removeSnp)]
print("donersids length is",length(doneRsids))

write.table(doneRsids,"doneRsids",quote=F,row.names=F,col.names=F)

