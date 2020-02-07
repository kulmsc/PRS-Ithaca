args <- commandArgs(trailingOnly = TRUE)

x <- read.table("preSummStat", stringsAsFactors = FALSE)
y <- read.table(args[1], stringsAsFactors = FALSE, header = T)

x <- x[complete.cases(x),]
y <- y[complete.cases(y),]

x <- x[x$V2 %in% y$sid,]
y <- y[y$sid %in% x$V2,]

x <- x[order(x$V2),]
y <- y[order(y$sid),]

x$V8 <- y$ldpred_beta

write.table(x, "newSummStat", row.names = F, col.names = F, quote = F, sep = '\t')
