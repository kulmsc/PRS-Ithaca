library(bigsnpr)

args = commandArgs(trailingOnly=TRUE)
author <- args[1]
chr <- args[2]
allResponses <- read.table("~/prsDatabase/allResponses",stringsAsFactors=F,header=T)

#snp_readBed("new.22.filter.bed")
obj.bigSNP <- snp_attach(paste0("../new.",chr,".filter.rds"))
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y <- allResponses[,colnames(allResponses)==author]

sumstats <- read.table("preSummStat",stringsAsFactors=F)

#should this number be changed as a parameter?
colnames(sumstats) <- c("chr","rsid","pos","a0","a1","freq","info","beta","se","p","n","freq2")

map <- obj.bigSNP$map[,-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map)
beta <- info_snp$beta
lpval <- -log10(info_snp$p)

pullInds <- sort(c(sample(which(y==0),sum(y==1)*5),which(y==1)))
y2 <- y[pullInds]
G2 <- big_copy(X=G,ind.row=pullInds,ind.col=match(info_snp$rsid,obj.bigSNP$map$marker.ID),backingfile="./subset")
CHR2 <- rep(obj.bigSNP$map$chromosome[1],ncol(G2))
POS2 <- obj.bigSNP$map$physical.pos[obj.bigSNP$map$marker.ID %in% info_snp$rsid]

replace_snp <- data.frame(info_snp)
replace_snp <- replace_snp[,c(1,5,2,3,4,6,7,8,9,10,11,12)]


i=1
for(trainFrac in c(0.25,0.5,0.75)){
  trainSize <- round(nrow(G2)*trainFrac)
  ind.train <- sample(nrow(G2), trainSize)
  ind.test <- setdiff(rows_along(G2), ind.train)

  all_keep <- snp_grid_clumping(G2, CHR2, POS2, ind.row = ind.train,
                              lpS = lpval, ncores = 4)
  multi_PRS <- snp_grid_PRS(G2, all_keep, beta, lpval, ind.row = ind.train,
                          n_thr_lpS = 50, ncores = 4)
  final_mod <- snp_grid_stacking(multi_PRS, y2[ind.train], ncores = 4, K = 4)
  replace_snp$beta <- final_mod$beta.G
  write.table(replace_snp,paste0("summStat",i),col.names=T,row.names=F,quote=F,sep='\t')
  i=i+1
}
