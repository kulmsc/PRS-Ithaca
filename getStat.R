ss=read.table("comboSummStat",header=T,stringsAsFactors=F)
tVal=qt(ss$P,df=ss$n)
se=ss$Direction/abs(tVal)
stat=(ss$Direction*se)^2
newSS=cbind(ss[,1:4],stat,ss[,6])
colnames(newSS)=c("Predictor","A1","A2","Direction","Stat","n")
write.table(newSS,"statSummStat",col.names=T,row.names=F,sep=' ',quote=F)
