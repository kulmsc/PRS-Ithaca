
dir=`pwd|cut -f4 -d'/' | cut -f2 -d'r'`
shortFile=$1
firstFile=`echo $shortFile | cut -f1-3 -d'.'`
file=${firstFile}.gz
chr=`echo $shortFile | cut -f3 -d'.'`

infoLim=0.3
zcat /home/cbsumezey/sdk2004/chr${chr}/$file | awk -v var="$infoLim" '$7 > var {print $0}' | \
        awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | \
        awk 'length($4) == 1 {print $0}' | awk 'length($5) == 1 {print $0}' | \
        awk '$10 > 0 {print $0}' | fgrep -v NA | sort -k2 | rev | uniq -f8 -u | rev > preSummStat
fgrep -w -f ../goodRsids.${chr} preSummStat > temp; mv temp preSummStat

cp /home/sdk2004/prsDatabase/makeLDPredStandard.py .
cp /home/sdk2004/prsDatabase/metaData .
python makeLDPredStandard.py
authorName=`echo $file | cut -f1 -d'.'`
sampleSize=`cat metaData | grep $authorName | cut -f2`
numberSnps=`cat metaData | grep $authorName | cut -f3`
let ldr=numberSnps/4500
let otherDir=$dir+30
taskset -c $dir,$otherDir python ~/bin/LDpred.py coord --gf=/home/cbsumezey/sdk2004/1000genomes/eur.chr$chr --ssf=summStat --N=$sampleSize --out=madeCoord --ssf-format=STANDARD
mv madeCoord goCoord.hd5

h5dump -d "/cord_data/chrom_${chr}/sids" goCoord.hd5 | fgrep rs | cut -f2 -d'"' | cut -f1 -d'\' > compRsid

paste /workdir/sdk2004/sets/chr${chr}/$shortFile compRsid | sort -k2 > comboLdpred
cat preSummStat | fgrep -w -f compRsid | sort -k2 > forMakeSummStat
cat forMakeSummStat | cut -f1-7 > firstPart
cat forMakeSummStat | cut -f9-10 > secondPart
cat comboLdpred | sed 's/ \+/\t/g' | cut -f2 > newBetas
paste firstPart newBetas secondPart > newSet
mv newSet /workdir/sdk2004/sets/chr${chr}/$shortFile
