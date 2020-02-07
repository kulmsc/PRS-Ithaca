cp ~/prsDatabase/keepEID .
cp ~/prsDatabase/shortPheno.csv .
Rscript makePheno.R


for num in {1..22}; do
        echo chromosome $num

        cp /bscb/data/500k_ukb/imputed/ukb_imp_chr${num}_v2.bgen chr${num}.bgen
        cp /bscb/data/500k_ukb/imputed/ukb_bgi_chr${num}_v2.bgi chr${num}.bgen.bgi
        cp /bscb/bscb09/500k_ukb/imputed/ukb1994_imp_chr${num}_v2_s487406.sample chr${num}.sample
	cat /home/sdk2004/ukData/bgen/listOnly/chr$num | cut -f2 | grep rs > rsids
	
	numLines=$(cat rsids | wc -l)
	let "splitLines = $numLines / 32" 
	split -l $splitLines -d -a 1 rsids rsid.
	split -l $splitLines -d rsids rsid.
	rm rsid.0*


	for i in {1..32}; do
		./gwas.sh $i $num &
		sleep 30
	done
	wait

	rm *bgen
	rm *bgi
	rm *sample

done
