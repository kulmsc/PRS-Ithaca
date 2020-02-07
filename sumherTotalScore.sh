echo "begin scoring"
chr=$1
file=$2
score=$3
infoLim=0.3

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
if [ -f  ready.$chr ]; then
	echo all prepared
elif [ -f willBe.${chr} ]; then
	keepGoing=True
	while [ $keepGoing == "True" ]; do
		if [ -f ready.$chr ]; then
			keepGoing=False
		else
			sleep $(( ( RANDOM % 10 )  + 1 ))
		fi
	done
else
	echo yep > willBe.${chr}
	plim=`cat ~/prsDatabase/subsetCLUMP | cut -f1 -d'-' | sort -n | tail -1`
	numBgens=`ls *bgen | wc -l`
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

	if [ $score == "clump" | $score == "winnersCurse-2d" ]; then
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
	else
		ls /home/cbsumezey/sdk2004/chr${chr} | grep gz$ | grep -v split | while read rsLine; do
			zcat /home/cbsumezey/sdk2004/chr${chr}/$rsLine | awk -v var="$infoLim" '$7 > var {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | cut -f2 >> allRsid.$chr
		done
	fi

	echo this is allRsid
	head allRsid.$chr
	echo this is ls
	lst
	echo done with allRsid
	cat allRsid.$chr | sort | uniq > temp.$chr ; mv temp.$chr allRsid.$chr
	bgenix -g chr${chr}.bgen -incl-rsids allRsid.$chr > new.${chr}.bgen
	rm chr${chr}.bgen chr${chr}.bgen.bgi
	plink2 --bgen new.${chr}.bgen --sample chr${chr}.sample --make-bed --out new.temp.${chr}
	rm chr${chr}.sample allRsid.$chr
	cat new.temp.${chr}.bim | cut -f2 | sort | uniq -d > dupIDs.$chr
	plink --bfile new.temp.${chr} --keep-fam phase.eid --exclude dupIDs.$chr --make-bed --out new.${chr}
	rm new.${chr}.bgen new.temp.${chr}.*
	if [ $score == "ldpred" ]; then
		cp /home/cbsumezey/sdk2004/1000genomes/eur.chr${chr}.*  .
	elif [ $score == "grabld" ]; then
		cat new.${chr}.fam | cut -f1 | sort -R | head -1000 > fam500.${chr}
		plink --bfile new.$chr --keep-fam fam500.${chr} --recode A --out forLD.${chr}
		cp ~/prsDatabase/ldGrabBLD.R ldGrabBLD.${chr}.R
		Rscript ldGrabBLD.${chr}.R $chr
		rm forLD.${chr}.raw ldGrabBLD.${chr}.R
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
			grep -f doneRsids summStat > doneSummStat
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
	python PRScs.py --ref_dir=/home/sdk2004/ldblk_1kg --bim_prefix=ss --sst_file=preSummStat --n_gwas=$sampSize --out_dir=. --chrom=$chr

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
        mv summStat ../sets/chr${chr}/${baseName}.1

elif [ $score == "tweedy" ]; then
	cp ~/prsDatabase/tweedy.R .
	cp ~/prsDatabase/metaData .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > preSummStat
	Rscript tweedy.R $file
	baseName="${file::-3}"	

	for i in {1..3}; do
		plink --bfile ../new.${chr}  --score summStat${i} 2 4 8 no-sum
        	if [ -f plink.profile ]; then
	                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${chr}.${i}
	                gzip ../store${dir}/${baseName}.${chr}.${i}
	        fi
	        mv summStat ../sets/chr${chr}/${baseName}.${i}
	done

elif [ $score == "winnersCurse-2d" ]; then
	cp /home/cbsumezey/sdk2004/genomeAnnotations/ptSnps .
	cp /home/cbsumezey/sdk2004/genomeAnnotations/conservedSnps .
	cp ~/prsDatabase/fileHeader .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | sort -k2 | rev | uniq -f8 -u | rev | awk '$10 > 0 {print $0}' > summStat
        cat summStat | fgrep -f ptSnps > ptSummStat
	cat summStat | fgrep -f conservedSnps > conSummStat
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
				mv plink.profile ss.score.$plim
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
						paste ss.score.$plim2 plink.profile | awk '{print $1+$2}' > ../store${dir}/${baseName}.${i}
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

elif [ $score == "sumher" ]; then
	cp ~/prsDatabase/metaData .
	cp ~/prsDatabase/sumherHeader .
	cp ~/prsDatabase/mhc.snps .
	cp getStat.R .

	zcat /home/cbsumezey/sdk2004/chr${chr}/$file | awk -v var="$infoLim" '$7 > var {print $0}' | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' > preSummStat
	cat preSummStat | cut -f2,4,5,8,10 > cutSummStat
	authorName=`echo $file | cut -f1 -d'.'`
	numControls=`cat metaData | grep $authorName | cut -f3`
	numCases=`cat metaData | grep $authorName | cut -f4`
	her=`cat metaData | grep $authorName | cut -f5`
	let totalPop=$numControls+$numCases	
	lengthSummStat=`cat preSummStat | wc -l`
	yes $totalPop | head -n $lengthSummStat > popCol
	paste cutSummStat popCol > comboSummStat
	cat comboSummStat | tr '\t' ' ' > newSummStat
	cat sumherHeader newSummStat > comboSummStat


	#Part 1 - may have to combine chromos to get better heritability estimations
	cat comboSummStat | cut -f1 -d' ' | sort | uniq -d > toRemove
	cat comboSummStat | fgrep -v -f toRemove > temp; mv temp comboSummStat
	Rscript getStat.R
	awk < statSummStat '(NR>1 && $5>$6/99){print $1}' > statBig
	if [ `cat statBig | wc -l` -gt 0 ];
		../ldak5.linux --remove-tags toGo --bfile ../ref --top-preds statBig --window-kb 1200 --min-cor .1
		cat mhc.snps toGo.out > excl
	else
		mv mhc.snps excl
	fi

	../ldak5.linux --cut-weights work_chr$chr --bfile ../ref --extract statSummStat --chr $chr
	../ldak5.linux --calc-weights-all work_chr$chr --bfile ../ref --extract statSummStat --chr $chr
	mv work_chr${chr}/weights.short work.weights
	../ldak5.linux --calc-tagging work_ldak --bfile ../ref --extract statSummStat --weights work.weights --power -.25 --window-kb 1200
	../ldak5.linux --sum-hers work_ldak --tagfile work_ldak.tagging --summary statSummStat --exclude excl

	#Part 2
	../ldak5.linux --calc-tagging work_ldak.ann --bfile ../ref --extract statSummStat --weights work.weights --power -.25 --window-kb 1200 --annotation-number 24 --annotation-prefix ../annSnps/ann_snps.
	../ldak5.linux --sum-hers work_ldak.ann --tagfile work_ldak.ann.tagging --summary statSummStat --exclude excl
	for i in {1..24}; do 
		grep "Enrich_A$i " work_ldak.ann.enrich | awk -v i=$i '{a+=$5/$6/$6;b+=1/$6/$6}END {print "Annotation:", i, "Av. Enr.:", a/b, "SD:", 1/sqrt(b), "Num. Traits:", NR}'
	done
	rm props.ldak
	for i in {1..24}; do 
		grep "Share_A$i " work_ldak.ann.share | awk -v i=$i '{a+=$2/$3/$3;b+=1/$3/$3}END {print a/b, 1/sqrt(b)}' >> props.ldak
	done
	grep Share_Base work_ldak.ann.share | awk -v i=$i '{a+=$2/$3/$3;b+=1/$3/$3}END {print a/b, 1/sqrt(b)}' >> props.ldak


	#Part 3
	#all been done!
	#../ldak5.linux --cut-weights ss_chr$chr --bfile ../ref --extract comboSummStat --chr $chr
	#../ldak5.linux --calc-weights-all ss_chr$chr --bfile ../ref --extract comboSummStat --chr $chr
	#mv ss_chr${chr}/weights.short work.weights
	../ldak5.linux --calc-tagging work_ldak.ann --bfile ../ref --extract comboSummStat --weights work.weights --power -.25 --window-kb 1200 --annotation-number 24 --annotation-prefix ../annSnps/ann_snps.
	../ldak5.linux --sum-hers work_ldak.ann --tagfile work_ldak.ann.tagging --summary comboSummStat --exclude excl
	#calculate the tagging
	../ldak5.linux --calc-tagging work_ldak.non --bfile ../ref --extract comboSummStat --weights work.weights --power -.25 --window-kb 1200 --reduce NO
	../ldak5.linux --calc-tagging work_ldak.ann.non --bfile ../ref --extract comboSummStat --weights work.weights --power -.25 --window-kb 1200 --annotation-number 24 --annotation-prefix ../annSnps/ann_snps. --reduce NO
	her=`grep Her_ALL work_ldak.hers | awk '{print $2}'`
	../ldak5.linux --calc-exps work_ldak.non --tagfile work_ldak.non.tagging --her $her
	#calculate single snp heritability
	#../ldak5.linux --calc-exps work_ldak.non --tagfile work_ldak.non.tagging --her $her
	#../ldak5.linux --calc-exps work_ldak.ann.non --tagfile work_ldak.ann.non.tagging --her $her


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
