#this goes in the total work part
chr=22
cat new.fam | cut -f1 | sort -R | head -10000 > fam500
plink --bfile new --keep-fam fam500 --recode A --out forLD
Rscript ldGrabBLD.R $chr

#above is only for the total
cat summStat | cut -f2 > specificRsids
plink --bfile new --extract specificRsids --exclude badRsids --recode A --out forGrabld
plink --bfile new --extract specificRsids --exclude badRsids --make-bed --out forGrabld

Rscript prepGwas.R
mv gwas.fam forGrabld.fam
plink --bfile forGrabld --logistic hide-covar beta --covar fullCovar.cov --out gwasRes

Rscript finalGrab.R
