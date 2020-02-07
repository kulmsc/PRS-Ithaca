echo "begin scoring"
chr=$1
file=$2
score=$3
infoLim=0.3


echo "THE SCORE IS"
echo $score

#DETERMINE THE DIRECTORY #######################
toOpen=`cat possAccess`
keepGoing=True
while [ $keepGoing == "True" ]; do
	if [ $toOpen == "True" ]; then
		echo False > possAccess
		dir=`head -1 possDirs`
		grep -v -w $dir possDirs > tempPoss; mv tempPoss  possDirs
		echo True > possAccess
		keepGoing=False
	else
		toOpen=`cat possAccess`
		sleep 1
	fi
done
echo we are in directory $dir
##################################################



#DETERMINE IF WHOLE NEW.BED MUST BE MADE ####################
#files are ready
if [ -f  ready.$chr ]; then
	echo all prepared
#in the process of making files
elif [ -f willBe.${chr} ]; then
	keepGoing=True
	while [ $keepGoing == "True" ]; do
		if [ -f ready.$chr ]; then
			keepGoing=False
		else
			sleep $(( ( RANDOM % 10 )  + 1 ))
		fi
	done
#will have to make files
else
	echo yep > willBe.${chr}
	plim=`cat ~/prsDatabase/subsetCLUMP | cut -f1 -d'-' | sort -n | tail -1`
	if [ $score == "nps" ]; then
                numBgens=`ls chr*bgen | wc -l`
        else
                numBgens=`ls *bgen | wc -l`
        fi
	checkGoing=True
	while [ $checkGoing == "True" ];do
		if [ $numBgens -ge 1 ]; then
			echo waiting for no bgens
			numBgens=`ls *bgen | wc -l`
			sleep $(( ( RANDOM % 10 )  + 1 ))
		else
			checkGoing=False
		fi
	done

	echo must do the prep work
	cp /bscb/data/500k_ukb/imputed/ukb_imp_chr${chr}_v2.bgen chr${chr}.bgen
        cp /bscb/data/500k_ukb/imputed/ukb_bgi_chr${chr}_v2.bgi chr${chr}.bgen.bgi
        cp /bscb/bscb09/500k_ukb/imputed/ukb1994_imp_chr17_v2_s487406.sample chr${chr}.sample

	echo done copying bgens

	if [ $score == "clump" ] || [ $score == "winnersCurse-2d" ]; then
		ls /home/cbsumezey/sdk2004/chr${chr} | grep gz$ | grep -v split | while read rsLine; do
			zcat /home/cbsumezey/sdk2004/chr${chr}/$rsLine | awk -v var="$infoLim" '$7 > var {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | awk -v var="$plim" '$10 < var {print $0}' | cut -f2 >> allRsid.$chr
		done
	elif [ $score == "sumher" ]; then
		ls /home/cbsumezey/sdk2004/chr${chr} | grep gz$ | grep -v split | while read rsLine; do
			zcat /home/cbsumezey/sdk2004/chr${chr}/$rsLine | awk -v var="$infoLim" '$7 > var {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | cut -f2 >> allRsid.$chr
		done
		if [ ! -f ref.bed ]; then
			cp /home/cbsumezey/sdk2004/1000genomes/ref* .
			cp ~/prsDatabase/ldak5.linux .
		fi
	elif [ $score == "nps" ]; then
		ls /home/cbsumezey/sdk2004/chr${chr} | grep gz$ | grep -v split | while read rsLine; do
                        zcat /home/cbsumezey/sdk2004/chr${chr}/$rsLine | awk '$7 > 0.4 {print $0}' | awk '$6 > 0.05 {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | cut -f2 >> allRsid.$chr
                done

	else
		ls /home/cbsumezey/sdk2004/chr${chr} | grep gz$ | grep -v split | while read rsLine; do
			zcat /home/cbsumezey/sdk2004/chr${chr}/$rsLine | awk -v var="$infoLim" '$7 > var {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | cut -f2 >> allRsid.$chr
		done
	fi

	echo done with allRsid
	if [ $score == "nps" ]; then
		echo "went nps"
		cat allRsid.$chr | sort | uniq > temp.$chr ; mv temp.$chr allRsid.$chr
                bgenix -g chr${chr}.bgen -incl-rsids allRsid.$chr > new.${chr}.bgen
                rm chr${chr}.bgen chr${chr}.bgen.bgi
		bgenix -g new.${chr}.bgen -index

		cp /bscb/data/500k_ukb/imputed/ukb_mfi_chr${chr}_v2.txt mfi${chr}
		cat mfi${chr} | cut -f1 > part1$chr
		cat mfi${chr} | cut -f1,2,3,4,5 > part2$chr
		cat mfi${chr} | cut -f4,6 > part3$chr
		paste part1$chr part2$chr part3$chr > new.${chr}.mfi
		#bgenix -g new.${chr}.bgen -list | tail -n +3 | head -n -1 > snps${chr}
		#cat snps${chr} | cut -f2 > snpsOnly${chr}
		#fgrep -w -f snpsOnly${chr} mfi${chr} > matchMfi${chr}
		
		if [ ! -d "nps" ]; then
			cp -r ~/prsDatabase/nps .
		fi
	else
		echo "went else"
		cat allRsid.$chr | sort | uniq > temp.$chr ; mv temp.$chr allRsid.$chr
		bgenix -g chr${chr}.bgen -incl-rsids allRsid.$chr > new.${chr}.bgen
		rm chr${chr}.bgen chr${chr}.bgen.bgi
		plink2 --bgen new.${chr}.bgen --sample chr${chr}.sample --make-bed --out new.temp.${chr}
		rm chr${chr}.sample allRsid.$chr
		cat new.temp.${chr}.bim | cut -f2 | sort | uniq -d > dupIDs.$chr
		plink --bfile new.temp.${chr} --keep-fam phase.eid --exclude dupIDs.$chr --make-bed --out new.${chr}
		rm new.${chr}.bgen new.temp.${chr}.*
		if [ $score == "ldpred" ] || [ $score == "ldFunct" ] || [ $score == "annoPred" ]; then
			cp /home/cbsumezey/sdk2004/1000genomes/eur.chr${chr}.*  .
		elif [ $score == "grabld" ]; then
			cat new.${chr}.fam | cut -f1 | sort -R | head -1000 > fam500.${chr}
			plink --bfile new.$chr --keep-fam fam500.${chr} --recode A --out forLD.${chr}
			cp ~/prsDatabase/ldGrabBLD.R ldGrabBLD.${chr}.R
			Rscript ldGrabBLD.${chr}.R $chr
			rm forLD.${chr}.raw ldGrabBLD.${chr}.R
		fi
	fi
	if [ $score == "annoPred" ]; then
		cp ~/prsDatabase/phenoFamFile new.${chr}.fam
	fi
	echo goodToGo > ready.$chr
fi
#####################################################################



#THE SCORING #####################################
cd dir$dir
echo onto scoring within dir$dir

if [ $score == "clump" ]; then
	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > summStat
	cp ~/prsDatabase/fileHeader .
	cp ~/prsDatabase/subsetCLUMP .
	cat fileHeader summStat > temp ; mv temp summStat
	baseName="${file::-3}"

	i=1
	cat subsetCLUMP | while read sub; do
		plim=`echo $sub | cut -f1 -d'-'`
		r2lim=`echo $sub | cut -f2 -d'-'`
		plink --bfile ../new.${chr}  --clump summStat --clump-p1 $plim --clump-r2 $r2lim
		if [ -f plink.clumped ]; then
			sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
			fgrep -f doneRsids summStat > doneSummStat
			plink --bfile ../new.${chr} --score doneSummStat 2 4 8 no-sum
			mv doneSummStat ../sets/chr${chr}/${baseName}.${i}
			if [ -f plink.profile ]; then
				sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${i}
				gzip ../store${dir}/${baseName}.${i}
				rm plink.profile
			fi
			rm plink.clumped
		fi
		let i=i+1
	done


elif [ $score == "ldpred" ]; then
	cp /home/sdk2004/prsDatabase/makeLDPredStandard.py .
        cp /home/sdk2004/prsDatabase/metaData .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > preSummStat
	python makeLDPredStandard.py
	authorName=`echo $file | cut -f1 -d'.'`
	sampleSize=`cat metaData | grep $authorName | cut -f2`
	numberSnps=`cat metaData | grep $authorName | cut -f3`
	let ldr=numberSnps/3000
	baseName="${file::-3}"
	echo before action
	ls

	coord --gf=../eur.chr$chr --ssf=summStat --N=$sampleSize --out=madeCoord
	ldpred --coord=madeCoord --ld_radius=$ldr --PS=0.3,0.2,0.1 --N=$sampleSize --out=donePred

	i=1
	ls *_p[0-9].* | while read pred; do
		plink --bfile ../new.${chr} --score $pred 3 4 7 no-sum header
		if [ -f plink.profile ]; then
			sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.${i}
			gzip ../store${dir}/${baseName}.${chr}.${i}
		fi
		mv $pred ../sets/chr${chr}/${baseName}.${i}
		let i=i+1
	done

elif [ $score == "report" ]; then
	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > summStat
	baseName="${file::-3}"
	plink --bfile ../new.${chr}  --score summStat 2 4 8 no-sum
	if [ -f plink.profile ]; then
        	sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.1
        	gzip ../store${dir}/${baseName}.${chr}.1
	fi

#note this should only be applied to summary stats that have already been clumped
elif [ $score == "winnersCurse-lasso" ]; then
	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' | grep -v NA > preSummStat
	baseName="${file::-3}"
	cp ~/prsDatabase/winnersCurse.py .

	i=1
	for lambda in 0.001 0.01 0.1; do
		python winnersCurse.py lasso $lambda
	        plink --bfile ../new.${chr} --score summStat 2 4 8 no-sum
	        if [ -f plink.profile ]; then
	                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.${i}
	                gzip ../store${dir}/${baseName}.${chr}.${i}
	        fi
		mv summStat ../sets/chr${chr}/${baseName}.${i}
		let i=i+1
	done


#note this should only be applied to summary stats that have already been clumped
elif [ $score == "winnersCurse-like" ]; then
        zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' | grep -v NA > preSummStat

        cp ~/prsDatabase/winnersCurse.py .
        python winnersCurse.py likelihood 0
        baseName="${file::-3}"
        plink --bfile ../new.${chr}  --score summStat 2 4 8 no-sum
        if [ -f plink.profile ]; then
                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.1
                gzip ../store${dir}/${baseName}.${chr}.1
        fi
	mv summStat ../sets/chr${chr}/${baseName}.1

elif [ $score == "grabld" ]; then
	cp ~/prsDatabase/fileDecoder .
	cp ~/prsDatabase/nameNumberDecoder .
	cp ~/prsDatabase/allScoresCompiled.train .
	cp ~/prsDatabase/covars .
	cp ~/prsDatabase/prepGwas.R .
	cp ~/prsDatabase/finalGrab.R .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > preSummStat
	cat preSummStat | cut -f2 > specificRsids
	fgrep -w -f specificRsids ../allRsids.$chr > finalRsids
	lenFinal=`cat finalRsids | wc -l`
	if [ $lenFinal -gt 5 ]; then
		Rscript prepGwas.R $chr $file
		plink --memory 3000 --bfile ../new.$chr --fam gwas.fam --extract finalRsids --logistic hide-covar beta --covar fullCovar.cov --out gwasRes
		Rscript finalGrab.R $chr

		baseName="${file::-3}"
		for i in {1..4}; do
			plink --bfile ../new.${chr}  --score summStat$i 2 4 8 no-sum
			if [ -f plink.profile ]; then
				sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.${i}
				gzip ../store${dir}/${baseName}.${chr}.${i}
	       	 	fi
			mv summStat$i ../sets/chr${chr}/${baseName}.${i}
		done
	fi

elif [ $score == "lassosum" ]; then
	cp ~/prsDatabase/fileDecoder .
        cp ~/prsDatabase/nameNumberDecoder .
        cp ~/prsDatabase/allScoresCompiled.train .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > preSummStat
	cat preSummStat | cut -f2 > specificRsids
	name=`echo $file | cut -f1 -d'.'`
	snps=`cat ~/prsDatabase/metaData | grep $name | cut -f3`	

	cp ~/prsDatabase/lassosum.R .
	cp ~/prsDatabase/lassosum.splitVal.R .
	split -l 80000 --additional-suffix .split preSummStat
	ls *.split | while read line; do
		Rscript lassosum.R $snps $chr $file $line
		cat summStat.temp.1 >> summStat1 ; cat summStat.temp.2 >> summStat2 ; cat summStat.temp.3 >> summStat3
	done


	baseName="${file::-3}"
        for i in {1..3}; do
        	plink --bfile ../new.${chr}  --score summStat$i 2 4 8 no-sum
                if [ -f plink.profile ]; then
                	sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.${i}
                        gzip ../store${dir}/${baseName}.${chr}.${i}
                fi
        	mv summStat$i ../sets/chr${chr}/${baseName}.${i}
	done

elif [ $score == "prsCS" ]; then
	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > fullSummStat
	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' | cut -f2,4,5,8,10 > preSummStat
	cat preSummStat | cut -f1 > specificRsids
	plink --bfile ../new.${chr} --extract specificRsids --make-just-bim --out ss

	name=`echo $file | cut -f1 -d'.'`
        sampSize=`cat ~/prsDatabase/metaData | grep $name | cut -f2`

	cp ~/PRScs/* .
	cat preSummStat >> summStatHeader; mv summStatHeader preSummStat
	i=1
	for phi in 0.000001 0.0001 0.01 1; do
		python PRScs.py --ref_dir=/home/sdk2004/ldblk_1kg --bim_prefix=ss --sst_file=preSummStat --n_gwas=$sampSize --out_dir=. --chrom=$chr --phi=$phi

		pstFile=`ls pst*`
		cat $pstFile | cut -f2 > prsRsids
		fgrep -f prsRsids fullSummStat > preSummStat
		cat preSummStat | cut -f1-7 > preSS1
		cat preSummStat | cut -f9-10 > preSS2
		cat $pstFile | cut -f6 > newBeta
		paste preSS1 newBeta preSS2 > summStat

		baseName="${file::-3}"
	        plink --bfile ../new.${chr}  --score summStat 2 4 8 no-sum
	        if [ -f plink.profile ]; then
	                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.1
	                gzip ../store${dir}/${baseName}.${chr}.1
	        fi
	        mv summStat ../sets/chr${chr}/${baseName}.${i}
		let i=i+1
	done

elif [ $score == "tweedy" ]; then
	cp ~/prsDatabase/tweedy.R .
	cp ~/prsDatabase/metaData .
        cp ~/prsDatabase/fileHeader .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' | awk '$6 > 0 {print $0}' > preSummStat
	cat fileHeader preSummStat > temp ; mv temp preSummStat

	plink --bfile ../new.${chr}  --clump preSummStat --clump-p1 0.05 --clump-r2 0.25
	sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
	fgrep -f doneRsids preSummStat > summStat


	Rscript tweedy.R $file
	baseName="${file::-3}"	

	for i in {1..3}; do
		plink --bfile ../new.${chr}  --score summStat${i} 2 4 8 no-sum
        	if [ -f plink.profile ]; then
	                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.${i}
	                gzip ../store${dir}/${baseName}.${chr}.${i}
	        fi
	        mv summStat${i} ../sets/chr${chr}/${baseName}.${i}
	done

elif [ $score == "winnersCurse-2d" ]; then
	echo "starting winnersCurse-2d"
	cp /home/cbsumezey/sdk2004/genomeAnnotations/ptSnps .
	cp /home/cbsumezey/sdk2004/genomeAnnotations/conservedSnps .
	cp ~/prsDatabase/fileHeader .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > summStat
        cat summStat | fgrep -f ptSnps > ptSummStat
	cat summStat | fgrep -f conservedSnps > conSummStat
	cat summStat | fgrep -v -f ptSnps | fgrep -v -f conservedSnps > temp ; mv temp summStat

	cat fileHeader ptSummStat > temp ; mv temp ptSummStat
	cat fileHeader conSummStat > temp ; mv temp conSummStat
	cat fileHeader summStat > temp; mv temp summStat
        baseName="${file::-3}"
	lengthFam=`cat ../new.${chr}.fam | wc -l`

	for plim in 0.00000001 0.01 0.1; do
		plink --bfile ../new.${chr}  --clump summStat --clump-p1 $plim --clump-r2 0.1
                if [ -f plink.clumped ]; then
                        sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
                        fgrep -f doneRsids summStat > doneSummStat
			rm plink.clumped
                        plink --bfile ../new.${chr} --score doneSummStat 2 4 8 no-sum
			mv doneSummStat stat.${plim}
                        if [ -f plink.profile ]; then
				sed -e 's/ [ ]*/\t/g' plink.profile | sed '/^\s*$/d' | cut -f7 | tail -n +2 > ss.score.$plim
			else
				yes 0 | head -n $lengthFam > ss.score.$plim
                        fi
		else
			yes 0 | head -n $lengthFam > ss.score.$plim
                fi
	done

	i=1
	for fileType in pt con; do
		for plim in 0.0001 0.1 0.5; do
	                plink --bfile ../new.${chr}  --clump ${fileType}SummStat --clump-p1 $plim --clump-r2 0.1
	                if [ -f plink.clumped ]; then
        	                sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
                	        fgrep -f doneRsids ${fileType}SummStat > doneSummStat
				rm plink.clumped
                        	plink --bfile ../new.${chr} --score doneSummStat 2 4 8 no-sum
                        	if [ -f plink.profile ]; then
					for plim2 in 0.00000001 0.01 0.1; do
						sed -e 's/ [ ]*/\t/g' plink.profile | sed '/^\s*$/d' | cut -f7 | tail -n +2 > newScore
						paste ss.score.$plim2 newScore | awk '{print $1+$2}' > ../store${dir}/${baseName}.${i}
						cat doneSummStat stat.$plim2 > ../sets/chr${chr}/${baseName}.${i}
						gzip ../store${dir}/${baseName}.${i}
						let i=i+1
					done
					rm plink.profile
                        	else
                                	yes 0 | head -n $lengthFam > ${fileType}.score.$plim
					for plim2 in 0.00000001 0.01 0.1; do
                                                paste ss.score.$plim2 plink.profile | awk '{print $1+$2}' > ../store${dir}/${baseName}.${i}
						cat doneSummStat stat.$plim2 > ../sets/chr${chr}/${baseName}.${i}
						gzip ../store${dir}/${baseName}.${i}
						let i=i+1
                                        done
                        	fi
                	fi
        	done
	done

elif [ $score == "nps" ]; then
	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | awk '$7 > 0.4 {print $0}' | awk '$6 > 0.05 {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' > preSummStat
	baseName="${file::-3}"
	cp ~/prsDatabase/npsHeader .
        cp ~/prsDatabase/fileDecoder .
        cp ~/prsDatabase/nameNumberDecoder .
        cp ~/prsDatabase/allScoresCompiled.train .
	cp ~/prsDatabase/npsPhenoMaker.R .
	cp ~/prsDatabase/fixSS.R .
	cp ~/prsDatabase/minHeader .
	Rscript npsPhenoMaker.R $chr $file #need to get train and valid files

	#common SNPs
	cat preSummStat | cut -f2 > inclRsids
	bgenix -g ../new.${chr}.bgen -incl-rsids inclRsids > small.bgen
	
	#filter samples
	cat valid.fam | cut -f1 -d' ' > inclValid
	qctool_v2.1-dev -g small.bgen -s ../chr${chr}.sample -incl-samples inclValid -og chrom${chr}.valid.QC1.bgen -os chrom${chr}.valid.QC1.sample
	qctool_v2.1-dev -g chrom${chr}.valid.QC1.bgen -s chrom${chr}.valid.QC1.sample -ofiletype dosage -og chrom${chr}.valid.QC1.dosage.gz

	#harmonize
	cat preSummStat | cut -f1,3,4,5 > part1; cat preSummStat | cut -f4,10 > part2; cat preSummStat | cut -f6,8 > part3
	paste part1 part2 part3 > minSumStat	
	cat minHeader minSumStat > temp; mv temp minSumStat
	Rscript ../nps/ukbb_support/harmonize_summstats.R minSumStat . valid $chr

	#filter variants
	zgrep -v -w -F -f valid.UKBB_rejected_SNPIDs chrom${chr}.valid.QC1.dosage.gz | gzip -c > temp; mv temp chrom${chr}.valid.QC1.dosage.gz
	grep -v -w -F -f valid.UKBB_rejected_SNPIDs chrom${chr}.valid.QC1.dosage.gz ../mfi$chr > chrom${chr}.ukb_mfi.QC2.txt

	#make fam file
	zcat chrom${chr}.valid.QC1.dosage.gz | head -n 1 | tr ' ' "\n" | tail -n +7 | awk '{ print( $1 " " $1 " 0 0 0 -9" ) }' > valid.QC1.fam

	#run the analysis
	cp ~/prsDatabase/stdgt .
	prefix=train
	N=`zcat chrom${chr}.valid.QC1.dosage.gz | head -n 1 | tr " " "\n" | tail -n +7 | wc -l`
	M=`zcat chrom${chr}.valid.QC1.dosage.gz | tail -n +2 | wc -l | cut -d' ' -f1`
	zcat chrom${chr}.valid.QC1.dosage.gz | ./stdgt $N $M chrom${chr}.valid.QC1
	gzip -f chrom${chr}.valid.QC1.stdgt	#run this program, normalizes the genotypes

	Rscript ../nps/npsR/nps_init.R summStat . train.fam train.phen valid.QC1 80 . $chr
	
	cp train.stdgt.gz chrom${chr}.train.stdgt.gz
	i=1
	for winShift in 0 20 40 60; do
		Rscript ../nps/npsR/nps_decor.R . $chr $winShift

		Rscript ../nps/npsR/nps_prune.R . $chr $winShift

		Rscript ../nps/npsR/nps_gwassig.R . $chr $winShift

		Rscript ../nps/npsR/nps_prep_part.R . $winShift 10 10
	
		Rscript ../nps/npsR/nps_part.R . $chr $winShift

		Rscript ../nps/npsR/nps_weight.R . $winShift $chr
	
		Rscript ../nps/npsR/nps_back2snpeff.R . $chr $winShift

		paste keep1 train.win_${winShift}.adjbetahat.chrom${chr}.txt keep2 > totalSummStat
		cp totalSummStat ../sets/chr${chr}/${baseName}.${i}
		let i=i+1
	done
	rm -r log

elif [ $score == "annoPred" ]; then
	baseName="${file::-3}"
	cp ~/prsDatabase/fixSSAnnoPred.R .
	cp ~/prsDatabase/fileDecoder .
	cp ~/prsDatabase/nameNumberDecoder .
	cp ~/prsDatabase/metaData .
	cp ~/prsDatabase/annopredPhenoMaker.R .
	cp ~/prsDatabase/allScoresCompiled.train .
	#Rscript annopredPhenoMaker.R $chr $file

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > preSummStat	
	Rscript fixSSAnnoPred.R

	author=`echo $file | cut -f1 -d'.'`
	sampSize=`cat metaData | grep $author | cut -f2`
	ln -s /home/cbsumezey/sdk2004/AnnoPred/ref ref
	split -l 25000 summStat newSS
	
	rm compBeta*
	ls newSS* | while read newSS; do
		rm test* sumstats* tier*
		python /home/cbsumezey/sdk2004/AnnoPred/AnnoPred.py --sumstats=${newSS} --ref_gt=../eur.chr${chr} --val_gt=../new.${chr} --coord_out=testCOORD --N_sample=${sampSize} --annotation_flag="tier3" --P=0.1 --local_ld_prefix=testLD --out=test --temp_dir=.
		cat test_h2_non_inf_betas* | tail -n +2 >> compBeta1
		cat test_h2_inf_betas* | tail -n +2 >> compBeta2
		cat test_pT_non_inf_betas* | tail -n +2 >> compBeta3
		cat test_pT_inf_betas* | tail -n +2 >> compBeta4	
	done

	cat preSummStat | cut -f1-7 > part1
        cat preSummStat | cut -f9-10 > part2
        i=1
        ls compBeta* | while read line; do
		sed 's/ \+ /\t/g' $line | cut -f7 > justBeta
                paste part1 justBeta part2 > summStat
                plink --bfile ../new.${chr}  --score summStat 2 4 8 no-sum
                if [ -f plink.profile ]; then
                        sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.${i}
                        gzip ../store${dir}/${baseName}.${chr}.${i}
                fi
                mv summStat ../sets/chr${chr}/${baseName}.${i}
		let i=i+1
        done

	


fi

rm *

cd ..
###################################################


echo $dir >> possDirs
echo done >> dones.$chr
numDone=`cat dones.$chr | wc -l`
toBeDone=`cat /home/cbsumezey/sdk2004/chr${chr} | grep gz$ | grep -v split | wc -l`
if [ $numDone == $toBeDone ]; then
	rm new.${chr}.*
fi
