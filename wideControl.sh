score=$1
if [ $# -eq 0 ]; then
        echo "no argument"
else

for num in {10..22}; do
	echo chromosome $num

	cp /bscb/data/500k_ukb/imputed/ukb_imp_chr${num}_v2.bgen chr${num}.bgen 
	cp /bscb/data/500k_ukb/imputed/ukb_bgi_chr${num}_v2.bgi chr${num}.bgen.bgi 
	cp /bscb/bscb09/500k_ukb/imputed/ukb1994_imp_chr${num}_v2_s487406.sample chr${num}.sample	

	for i in {1..48}; do
		cp /home/cbsumezey/sdk2004/chr${num}/*split.${i}.gz dir$i
		cp /home/cbsumezey/sdk2004/chr${num}/*wide.${i}.gz dir$i
       	done


	for i in {1..48}; do
		./wideScore.sh $i $num $score &
		sleep 10
	done
	wait

	rm chr${num}.bgen
	rm chr${num}.bgen.bgi
	rm chr${num}.sample
	for i in {1..32};  do
		rm dir${i}/*
	done
done

fi
