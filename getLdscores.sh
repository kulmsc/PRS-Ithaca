file=Liu_1.sumher.22.1.gz
startName=`echo $file | cut -f1,2 -d'.'`
endName=`echo $file | cut -f4,5 -d'.'`

echo $startName
echo $endName

for chr in {1..22}; do
	echo ${startName}.${chr}.${endName}
	zcat /home/cbsumezey/sdk2004/chr${chr}/${startName}.${chr}.${endName} | awk '$7 > 0.9 {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | awk '$10 > 0 {print $0}' >> ${startName}.full
done

length=`cat ${startName}.full | wc -l`
authorName=`echo $file | cut -f1 -d'.'`
sampSize=`cat metaData | grep $authorName | cut -f2`
yes $sampSize | head -$length > sampSize

cat ${startName}.full | cut -f2,4,5 > part1
cat ${startName}.full | cut -f10,8 > part2
paste part1 sampSize part2 > ${startName}.done
