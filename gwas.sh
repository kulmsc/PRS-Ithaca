num=$1
chr=$2

cd dir$num

mv ../rsid.$num .
bgenix -g ../chr${chr}.bgen -incl-rsids rsid.$num > new.bgen
plink2 --bgen new.bgen --sample ../chr${chr}.sample --make-bed --out new
rm new.bgen
cat new.bim | cut -f2 | uniq -d > dupIDs

plink --bfile new --exclude dupIDs --keep-fam ../keepEID --pheno ../phenoCoding --all-pheno --allow-no-sex --assoc --out assoc.${chr}.${num}
rm new.*


cd ..
