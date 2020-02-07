
args = commandArgs(trailingOnly=TRUE)
frq <- read.table("our.frq",stringsAsFactors=F)
ss <- read.table(args[1],stringsAsFactors=F)

frq <- frq[order(frq[,2]),]
ss <- ss[order(ss[,2]),]

ss[ss[,4]==frq[,3] & ss[,5]==frq[,4],6] <- frq[ss[,4]==frq[,3] & ss[,5]==frq[,4],5]
ss[ss[,4]==frq[,4] & ss[,5]==frq[,3],6] <- frq[ss[,4]==frq[,4] & ss[,5]==frq[,3],5]

write.table(ss,args[1],quote=F,sep='\t',row.names=F,col.names=F)
