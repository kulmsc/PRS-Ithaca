for num in {22..22}; do
	echo chromosome $num
	ls /home/cbsumezey/sdk2004/chr22 | grep -v split | while read ssFile; do
		cp /home/cbsumezey/sdk2004/chr${num}/$ssFile prep
	done

	#cp /bscb/data/500k_ukb/imputed/ukb_imp_chr${num}_v2.bgen chr${num}.bgen 
	#cp /bscb/data/500k_ukb/imputed/ukb_bgi_chr${num}_v2.bgi chr${num}.bgen.bgi 
	#cp /bscb/bscb09/500k_ukb/imputed/ukb1994_imp_chr${num}_v2_s487406.sample chr${num}.sample	

	cd prep
	cp ~/prsDatabase/comboRemoveEIDS .
	cp ~/prsDatabase/fileHeader .
	i=1
	ls *gz | while read line; do
		echo i $i
		../clump.sh $line $num $i &
		sleep 60
		let i++
		joblist=($(jobs -p))
		while (( ${#joblist[*]} >= 32 )); do
			sleep 2
			joblist=($(jobs -p))
		done	
	done
	rm `ls | grep -v gz`
	python ../bestSplit.py
	for i in {1..32}; do
		mv *.split.${i}.* ../dir$i
	done
	cd ..


	#for i in {1..32}; do
	#	./score.sh $i $num &
	#	sleep 60
	#done
	#wait

	#rm chr${num}.bgen
	#rm chr${num}.bgen.bgi
	#rm chr${num}.sample
	#for i in {1..32};  do
	#	rm dir${i}/*
	#done
done
