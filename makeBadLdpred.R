x <- read.table("preSummStat", stringsAsFactors=F)
y <- data.frame(x)
y <- y[order(y[,3]),]

#sorted by rsid - sorted by pos
#x[,8] <- y[,8]
#write.table(x, "preSummStat", row.names=F, col.names=F, quote=F, sep='\t')

#sorted by character position
x <- x[order(as.character(x[,3])), ]
x[,8] <- y[,8]
write.table(x, "preSummStat", row.names=F, col.names=F, quote=F, sep='\t')
