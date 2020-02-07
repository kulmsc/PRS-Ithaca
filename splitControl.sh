score=$1

if [ $# -eq 0 ]; then
        echo "no argument"
else


#maxVal=0
#while read line; do
#	val=$(echo $line | cut -d'.' -f5)
#	if [ $val > $maxVal ]; then
#		maxVal=$val
#	fi
#done <<< "$(ls /home/cbsumezey/sdk2004/chr1/*split*)"
#echo MAXVAL
#echo $maxVal

#rm [0-9][0-9]
#rm [0-9]

rm possMax
ls /home/cbsumezey/sdk2004/chr1/*split* | while read line; do
	echo $line | cut -d'.' -f5 >> possMax
done
maxVal=`cat possMax | sort -n | tail -1`

echo MAXVAL
echo $maxVal

for num in {22..22}; do
	echo chromosome $num

	cp /bscb/data/500k_ukb/imputed/ukb_imp_chr${num}_v2.bgen chr${num}.bgen 
	cp /bscb/data/500k_ukb/imputed/ukb_bgi_chr${num}_v2.bgi chr${num}.bgen.bgi 
	cp /bscb/bscb09/500k_ukb/imputed/ukb1994_imp_chr${num}_v2_s487406.sample chr${num}.sample	

	#for i in {1..$maxVal}; do
	for (( i=1; i<=$maxVal; i++ )); do
		cp /home/cbsumezey/sdk2004/chr${num}/*split.${i}.gz dir$i
       	done


	#for i in {1..$maxVal}; do
	for (( i=1; i<=$maxVal; i++ )); do
		./splitScore.sh $i $num $score &
		sleep 10
	done
	wait

	rm chr${num}.bgen
	rm chr${num}.bgen.bgi
	rm chr${num}.sample
	#for i in {1..3};  do
	for (( i=1; i<=$maxVal; i++ )); do
		rm dir${i}/*
	done
done

fi

date
