fileName=$1
chrNum=$2
indexNum=$3
baseName="${fileName::-3}"


gunzip $fileName
cat $baseName | cut -f2 | tail -n +2 > rsids.$indexNum
bgenix -g ../chr${chrNum}.bgen -incl-rsids rsids.$indexNum > new.${indexNum}.bgen 2> bgenLog
plink2 --bgen new.${indexNum}.bgen --sample ../chr${chrNum}.sample --make-bed --out new.${indexNum} &> plink2Log
rm new.${indexNum}.bgen

cat new.${indexNum}.bim | cut -f2 | uniq -d > dupIDs.${indexNum}
cat fileHeader $baseName > temp.$indexNum ; mv temp.$indexNum ${baseName}.head
plink --bfile new.${indexNum} --exclude dupIDs.${indexNum} --remove-fam comboRemoveEIDS --clump ${baseName}.head --clump-r2 0.2 --out ${baseName} &> plinkLog
rm new.${indexNum}.bed
rm new.${indexNum}.bim
rm new.${indexNum}.fam


sed "s/^[ \t]*//" ${baseName}.clumped | sed 's/ \+ /\t/g' | sed '/^\s*$/d' | cut -f3 > temp.$indexNum
grep -f temp.$indexNum $baseName > temp2.$indexNum ; mv temp2.$indexNum $baseName
gzip $baseName
