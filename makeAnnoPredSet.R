args <- commandArgs(trailingOnly = TRUE)

x <- read.table("preSummStat", stringsAsFactors=F)
y <- read.table(args[1], stringsAsFactors=F)

x <- x[complete.cases(x),]
y <- y[complete.cases(y),]

x <- x[x[,2] %in% y[,3],]
y <- y[y[,3] %in% x[,2],]

x <- x[order(x[,2]),]
y <- y[order(y[,3]),]

x[,8] <- y[,7]
write.table(x,"outSummStat", row.names=F, col.names=F, quote=F, sep='\t')
