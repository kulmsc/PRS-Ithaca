echo "begin scoring"
chr=$1
file=$2
score=$3
infoLim=0.3
mafLim=0.01
baseName="${file::-3}"
export OMP_NUM_THREADS=4

echo "THE SCORE IS"
echo score $score

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
			sleep $(( ( RANDOM % 30 )  + 1 ))
		fi
	done
#will have to make files
else
	echo yep > willBe.${chr}
	plim=`cat ~/prsDatabase/subsetCLUMP | cut -f1 -d'-' | sort -n | tail -1`
        numBgens=`ls *bgen | wc -l`
	checkGoing=True
	while [ $checkGoing == "True" ];do
		if [ $numBgens -ge 1 ]; then
			echo waiting for no bgens
			numBgens=`ls *bgen | wc -l`
			sleep $(( ( RANDOM % 30 )  + 1 ))
		else
			checkGoing=False
		fi
	done

	echo copying bgens
	#using v2
	#cp /bscb/data/500k_ukb/imputed/ukb_imp_chr${chr}_v2.bgen chr${chr}.bgen
        #cp /bscb/data/500k_ukb/imputed/ukb_bgi_chr${chr}_v2.bgi chr${chr}.bgen.bgi
	cp /bscb/bscb09/500k_ukb/imputed/ukb1994_imp_chr17_v2_s487406.sample chr${chr}.sample
	#using v3
	cp /bscb/ak735_0001/data/UK_Biobank/ukb_data/genetic/imputed/bgen_files/ukb_imp_chr${chr}_v3.bgen chr${chr}.bgen
	cp /home/cbsumezey/sdk2004/alexBgis/ukb_imp_chr${chr}_v3.bgen.bgi chr${chr}.bgen.bgi


	echo making list of rsids
	ls /home/cbsumezey/sdk2004/chr${chr} | grep gz$ | fgrep -v split | while read rsLine; do
		zcat /home/cbsumezey/sdk2004/chr${chr}/$rsLine | awk -v var="$infoLim" '$7 > var {print $0}' | \
		awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | \
		awk 'length($4) == 1 {print $0}' | awk 'length($5) == 1 {print $0}' | \
		awk '$10 > 0 {print $0}' | fgrep -v NA | fgrep -v nf | cut -f2 | sort | uniq -u >> allRsid.$chr
	done

	echo going through the plink processing steps
	bgenix -g chr${chr}.bgen -incl-rsids allRsid.$chr > new.${chr}.bgen
	rm chr${chr}.bgen chr${chr}.bgen.bgi
	plink2 --bgen new.${chr}.bgen --sample chr${chr}.sample --keep-fam phase.eid --make-bed --out new.temp.${chr}
	plink --bfile new.temp.${chr} --freq --out new.temp.${chr}
	sed 's/  */\t/g' new.temp.${chr}.frq  | tail -n +2 | awk '$5 < 0.01 {print $2}' | sort | uniq -u > badRsids.${chr}
	cat new.temp.${chr}.bim | cut -f2 | sort | uniq -d >> badRsids.${chr}
	plink --bfile new.temp.${chr} --exclude badRsids.${chr} --make-bed --out new.${chr}
	cat new.${chr}.bim | cut -f2 > goodRsids.${chr}
	rm new.${chr}.bgen new.temp.${chr}.* chr${chr}.sample

	echo goodToGo > ready.$chr
fi

#for ldpred LDPredFunct annoPred
if [ $score == "lassosum" ] || [ $score == "ldpred" ] || [ $score == "sbayesr" ]; then
	if [ ! -f "ukbb.ref.${chr}.bim" ]; then
		cp /home/cbsumezey/sdk2004/ukbbRef/ukbb.ref.${chr}.*  .
	fi
fi

if [ ! -d "LDpred-funct" ]; then
	cp -r ~/prsDatabase/LDpred-funct .
fi

#for grabld
if [ $score == "grabld" ]; then
	if [ ! -f "fam500.${chr}" ]; then
		cat new.${chr}.fam | cut -f1 | sort -R | head -1000 > fam500.${chr}
		plink --bfile new.$chr --keep-fam fam500.${chr} --recode A --out forLD.${chr}
		cp ~/prsDatabase/ldGrabBLD.R ldGrabBLD.${chr}.R
		Rscript ldGrabBLD.${chr}.R $chr
		rm forLD.${chr}.raw ldGrabBLD.${chr}.R
	fi
fi
#for stackCT
if [ $score == "stackCT" ]; then
	if [ ! -f new.${chr}.filter.rds ]; then
		if [ ! -f willBeStack.$chr ]; then
			echo willBeStack > willBeStack.$chr
			rm new.${chr}.filter*
			plink --bfile new.$chr --geno 0.01 --fill-missing-a2 --make-bed --out new.${chr}.filter
			Rscript ~/prsDatabase/makeStackRds.R $chr
			rm new.${chr}.filter.bed
			rm new.${chr}.filter.bim
			rm new.${chr}.filter.fam
		else
			while [ ! -f new.${chr}.filter.rds ];do
				sleep 30
			done
		fi
	fi
fi
#####################################################################



#THE SCORING #####################################
cd dir$dir
echo onto scoring within dir$dir
#to write down info
echo $dir > info
echo $chr >> info
echo $file >> info
echo $score >> info


#files shared by all methods
ls ../store* | fgrep ss > storeRecord
zcat /home/cbsumezey/sdk2004/chr${chr}/$file | awk -v var="$infoLim" '$7 > var {print $0}' | \
	awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' | \
	awk 'length($4) == 1 {print $0}' | awk 'length($5) == 1 {print $0}' | \
	awk '$10 > 0 {print $0}' | fgrep -v NA | fgrep -v nf | sort -k2 | rev | uniq -f8 -u | rev > preSummStat
fgrep -w -f ../goodRsids.${chr} preSummStat > temp; mv temp preSummStat



if [ $score == "clump" ]; then
	cp ~/prsDatabase/fileHeader .

	cat fileHeader preSummStat > temp ; mv temp preSummStat

	i=1
	for plim in 0.00000005 0.00005 0.05 0.5; do
		for r2lim in 0.25 0.5 0.75; do
			numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
                        if [ $numLine -eq 0 ]; then
				plink --threads 1 --bfile ../new.${chr} --clump preSummStat --clump-p1 $plim --clump-r2 $r2lim
				if [ -f plink.clumped ]; then
					sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
					fgrep -f doneRsids preSummStat > summStat
					plink --threads 1 --bfile ../new.${chr} --score summStat 2 4 8 no-sum
					mv summStat ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
					if [ -f plink.profile ]; then
						sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
						gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
						rm plink.profile
					fi
					rm plink.clumped
				fi
			fi
			let i=i+1
		done
	done


elif [ $score == "ldpred" ]; then
        echo head ldpred
        ls
	numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr} | wc -l`
        if [ $numLine -lt 5 ]; then

	cp /home/sdk2004/prsDatabase/makeLDPredStandard.py .
        cp /home/sdk2004/prsDatabase/metaData .

	python makeLDPredStandard.py
	authorName=`echo $file | cut -f1 -d'.'`
	sampleSize=`cat metaData | grep $authorName | cut -f2`
	numberSnps=`cat metaData | grep $authorName | cut -f3`
	let ldr=numberSnps/4500
	let otherDir=$dir+30

	python makeLDPredStandard.py --N=$sampleSize
	rm madeCoord*
	i=1

	taskset -c $dir,$otherDir ldpred coord --gf=/home/cbsumezey/sdk2004/1000genomes/eur.chr${chr} --ssf=summStat --N=$sampleSize --out=madeCoord --ssf-format=STANDARD
        #taskset -c $dir,$otherDir ldpred coord --gf=/home/cbsumezey/sdk2004/ukbbRef/ukbb.ref.${chr} --ssf=summStat --N=$sampleSize --out=madeCoord --ssf-format=STANDARD

	for f in 0.5 0.3 0.1 0.05 0.01; do
                echo f $f
		numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
        	if [ $numLine -eq 0 ]; then
                	taskset -c $dir,$otherDir ldpred gibbs --cf=madeCoord --ldr=$ldr --f=$f --N=$sampleSize --out=donePred --ldf=None
			ls donePred*_p* | fgrep -v inf | while read pred; do
                                echo $pred
				mv $pred forScore.$i
			done
		fi
		let i=i+1
        done

        echo MOVING ON
	ls forScore* | while read pred; do
                echo pred $pred
		i=`echo $pred | cut -f2 -d'.'`
		numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
                if [ $numLine -eq 0 ]; then
                        cp ~/prsDatabase/makeLDPredSet.R .
                        Rscript makeLDPredSet.R $pred
			plink --threads 1 --bfile ../new.${chr} --score newSummStat 2 4 8 sum
			if [ -f plink.profile ]; then
				sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
				gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
				mv newSummStat ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
			fi
		fi
		let i=i+1
	done
	fi

elif [ $score == "report" ]; then
	numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.1.gz | wc -l`
        if [ $numLine -eq 0 ]; then
	plink --threads 1 --bfile ../new.${chr}  --score preSummStat 2 4 8 sum
		if [ -f plink.profile ]; then
	        	sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.1
                        #sed 's/ \+/\t/g' plink.profile | cut -f5 | tail -n +2 > temp
                        #numSnp=`awk '{ total += $1 } END { print total/NR }' temp`
                        #echo $baseName $numSnp >> ../recAllSnps
	        	gzip ../store${dir}/${baseName}.${score}.${chr}.1
		fi
	fi


elif [ "$score" = "oldLdpred" ] ; then
        cp /home/sdk2004/prsDatabase/makeLDPredStandard.py .
        cp /home/sdk2004/prsDatabase/metaData .
        
	python makeLDPredStandard.py preSummStat #$wideName

	authorName=`echo $file | cut -f1 -d'.'`
        sampleSize=`cat metaData | grep $authorName | cut -f2`
        let otherDir=$dir+30

        snps=`cat summStat | wc -l`
        let ldr=snps/3000
	coord --gf=/home/cbsumezey/sdk2004/1000genomes/eur.chr${chr} --ssf=summStat --N=$sampleSize --out=madeCoord

	i=1
        for rho in  0.5 0.3 0.1 0.05 0.01; do
                ldpred --coord=madeCoord --ld_radius=$ldr --PS=$rho --N=$sampleSize --out summStat.adj
                cp ~/prsDatabase/makeLDPredSet.R .
		pred=`ls summStat.adj* | fgrep -v inf.txt`
                Rscript makeLDPredSet.R $pred

		plink --threads 1 --bfile ../new.${chr} --score newSummStat 2 4 8 sum
                if [ -f plink.profile ]; then
                        sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
                        gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
                        mv newSummStat ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
                fi
		rm summStat.adj*
		let i=i+1
        done





elif [ $score == "winnersCurseLasso" ]; then
	cp ~/prsDatabase/winnersCurse.py .
	cp ~/prsDatabase/fileHeader .
        cat fileHeader preSummStat > temp ; mv temp preSummStat
	plink --threads 1 --bfile ../new.${chr} --clump preSummStat --clump-p1 0.01 --clump-r2 0.1
	sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
        fgrep -f doneRsids preSummStat > temp; mv temp preSummStat

	i=1
	for lambda in 0.001 0.01 0.1; do
		numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
        	if [ $numLine -eq 0 ]; then
			python winnersCurse.py lasso $lambda
		        plink --threads 1 --bfile ../new.${chr} --score summStat 2 4 8 no-sum
		        if [ -f plink.profile ]; then
		                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
		                gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
		        fi
			mv summStat ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
		fi
		let i=i+1
	done


elif [ $score == "winnersCurseLike" ]; then
	numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.1.gz | wc -l`
        if [ $numLine -eq 0 ]; then
	        cp ~/prsDatabase/winnersCurse.py .
		cp ~/prsDatabase/fileHeader .
        	cat fileHeader preSummStat > temp ; mv temp preSummStat
		plink --threads 1 --bfile ../new.${chr} --clump preSummStat --clump-p1 0.01 --clump-r2 0.1
	        sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
        	fgrep -f doneRsids preSummStat > temp; mv temp preSummStat

	        python winnersCurse.py likelihood 0
	        plink --threads 1 --bfile ../new.${chr}  --score summStat 2 4 8 no-sum
	        if [ -f plink.profile ]; then
	                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.1
	                gzip ../store${dir}/${baseName}.${score}.${chr}.1
	        fi
		mv summStat ../sets/chr${chr}/${baseName}.${score}.${chr}.1
	fi


elif [ $score == "grabld" ]; then
	cp ~/prsDatabase/finalGrab.R .
	cp ~/prsDatabase/prepGwas.R .
	cp ~/prsDatabase/lassosumPhenoMaker.R .
	
	Rscript lassosumPhenoMaker.R $chr $file
	cat preSummStat | cut -f2 > specificRsids
	lenFinal=`cat specificRsids | wc -l`
	numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr} | wc -l`
        if [ $numLine -lt 4 ]; then

	if [ $lenFinal -gt 5 ]; then
		plink --bfile ../new.$chr --keep-fam our.fam --extract specificRsids --recode A --out forGwas

		Rscript finalGrab.R $chr

		for i in {1..4}; do
			numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
		        if [ $numLine -eq 0 ]; then
				if [ -f summStat$i ]; then
					plink --bfile ../new.${chr}  --score summStat$i 2 4 8 no-sum
					if [ -f plink.profile ]; then
						sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
						gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
			       	 	fi
					mv summStat$i ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
				fi
			fi
		done
	fi
	fi

elif [ $score == "lassosum" ]; then
	cp ~/prsDatabase/lassosumPhenoMaker.R .
	cp ~/prsDatabase/lsHeader .
	cp ~/prsDatabase/getLSBetas.R .

	Rscript lassosumPhenoMaker.R $chr $file
	cat preSummStat | cut -f2 > specificRsids
	plink --bfile ../new.$chr --keep-fam our.fam --extract specificRsids --make-bed --out test
	mv test.bed our.bed; mv test.bim our.bim

	plink --bfile ../ukbb.ref.$chr --extract specificRsids --make-bed --out eur
	cat lsHeader preSummStat > temp; mv temp preSummStat
	author=`echo $file | cut -f1 -d'.'`
	sampSize=`cat ~/prsDatabase/metaData | grep $author | cut -f2`	

	let otherDir=dir+30
	counter=1
	typeLasso=2
	numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr} | wc -l`
        if [ $numLine -lt 3 ]; then
		if [ $typeLasso -eq 1 ]; then
			taskset -c $dir,$otherDir lassosum --data preSummStat --chr chr --pos pos --A1 A1 --A2 A2 --pval pval --beta beta --n $sampSize \
				 --nthreads 1 --LDblocks EUR.hg19 --ref.bfile our 
		elif [ $typeLasso -eq 2 ]; then
			taskset -c $dir,$otherDir lassosum --data preSummStat --chr chr --pos pos --A1 A1 --A2 A2 --pval pval --beta beta --n $sampSize \
				 --nthreads 1 --LDblocks EUR.hg19 --ref.bfile our --test.bfile our
		elif [ $typeLasso -eq 3 ]; then
			taskset -c $dir,$otherDir lassosum --data preSummStat --chr chr --pos pos --A1 A1 --A2 A2 --pval pval --beta beta --n $sampSize \
				 --nthreads 1 --LDblocks EUR.hg19 --ref.bfile eur --test.bfile our
		fi

		lassosum --lassosum.pipeline lassosum.lassosum.pipeline.rds --pseudovalidate
		lassosum --lassosum.pipeline lassosum.lassosum.pipeline.rds --validate
		lassosum --lassosum.pipeline lassosum.lassosum.pipeline.rds --splitvalidate

		Rscript getLSBetas.R

        	for i in {1..3}; do
			numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${counter}.gz | wc -l`
		        if [ $numLine -eq 0 ]; then
        		plink --threads 1 --bfile ../new.${chr}  --score ss$i 2 4 8 no-sum
                	if [ -f plink.profile ]; then
                		sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${counter}
                        	gzip ../store${dir}/${baseName}.${score}.${chr}.${counter}
                	fi
        		mv ss$i ../sets/chr${chr}/${baseName}.${score}.${chr}.${counter}
			fi
			let counter=counter+1
		done
	fi

elif [ $score == "prsCS" ]; then
        cp ~/prsDatabase/makePRSCSSet.R .
	cat preSummStat | cut -f2,4,5,8,10 > goSummStat
	cat goSummStat | cut -f1 > specificRsids
	plink --threads 1 --bfile ../new.${chr} --extract specificRsids --make-bed --out ss

	name=`echo $file | cut -f1 -d'.'`
        sampSize=`cat ~/prsDatabase/metaData | grep $name | cut -f2`

	cp ~/PRScs/* .
	i=1
	b=0.5
	a=1
	for phi in 0.000001 0.01 1; do
		numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
		if [ $numLine -eq 0 ]; then
			let otherDir=dir+30
			cat summStatHeader goSummStat > smallSummStat
			taskset -c $dir,$otherDir python PRScs.py --ref_dir=/home/sdk2004/ldblk_1kg --bim_prefix=ss --sst_file=smallSummStat --n_gwas=$sampSize --out_dir=. \
				 --chrom=$chr --phi=$phi --a=$a --b=$b --n_burnin=100 --n_iter=300
			cat pst* >> fullPst
			rm pst*
		
			pstFile=fullPst
                        Rscript makePRSCSSet.R $pstFile

		        plink --bfile ../new.${chr}  --score outSummStat 2 4 8 no-sum
		        if [ -f plink.profile ]; then
		                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
		                gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
		        fi
		        mv outSummStat ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
			rm fullPst
		fi
		let i=i+1
	done

elif [ $score == "tweedy" ]; then
	cp ~/prsDatabase/tweedy.R .
        cp ~/prsDatabase/fileHeader .
	cp ~/prsDatabase/metaData .

	cat fileHeader preSummStat > temp ; mv temp preSummStat

	plink --threads 1 --bfile ../new.${chr}  --clump preSummStat --clump-p1 0.05 --clump-r2 0.25
	sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
	fgrep -f doneRsids preSummStat > summStat

	Rscript tweedy.R $file

	for i in {1..3}; do
		numLine=`cat ../storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
                if [ $numLine -eq 0 ]; then
			plink --bfile ../new.${chr}  --score summStat${i} 2 4 8 no-sum
	        	if [ -f plink.profile ]; then
		                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
		                gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
		        fi
		        mv summStat${i} ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
		fi
	done

elif [ $score == "winnersCurse-2d" ]; then
	cp /home/cbsumezey/sdk2004/genomeAnnotations/ptSnps .
	cp /home/cbsumezey/sdk2004/genomeAnnotations/conservedSnps .
	cp ~/prsDatabase/fileHeader .

        cat preSummStat | fgrep -f ptSnps > ptSummStat
	cat preSummStat | fgrep -f conservedSnps > conSummStat
	cat preSummStat | fgrep -v -f ptSnps | fgrep -v -f conservedSnps > temp ; mv temp preSummStat

	cat fileHeader ptSummStat > temp ; mv temp ptSummStat
	cat fileHeader conSummStat > temp ; mv temp conSummStat
	cat fileHeader preSummStat > temp; mv temp preSummStat
	lengthFam=`cat ../new.${chr}.fam | wc -l`

	for plim in 0.00000001 0.01 0.1; do
		plink --threads 1 --bfile ../new.${chr}  --clump preSummStat --clump-p1 $plim --clump-r2 0.1
                if [ -f plink.clumped ]; then
                        sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
                        fgrep -f doneRsids preSummStat > summStat
			rm plink.clumped
                        plink --threads 1 --bfile ../new.${chr} --score summStat 2 4 8 no-sum
			mv summStat stat.${plim}
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
	                plink --threads 1 --bfile ../new.${chr}  --clump ${fileType}SummStat --clump-p1 $plim --clump-r2 0.1
	                if [ -f plink.clumped ]; then
        	                sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
                	        fgrep -f doneRsids ${fileType}SummStat > summStat
				rm plink.clumped
                        	plink --threads 1 --bfile ../new.${chr} --score summStat 2 4 8 no-sum
                        	if [ -f plink.profile ]; then
					for plim2 in 0.00000001 0.01 0.1; do
						numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
				                if [ $numLine -eq 0 ]; then
							sed -e 's/ [ ]*/\t/g' plink.profile | sed '/^\s*$/d' | cut -f7 | tail -n +2 > newScore
							paste ss.score.$plim2 newScore | awk '{print $1+$2}' > ../store${dir}/${baseName}.${score}.${chr}.${i}
							cat summStat stat.$plim2 > ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
							gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
						fi
						let i=i+1
					done
					rm plink.profile
                        	else
                                	yes 0 | head -n $lengthFam > ${fileType}.score.$plim
					for plim2 in 0.00000001 0.01 0.1; do
						numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
                                                if [ $numLine -eq 0 ]; then
	       	                                         paste ss.score.$plim2 plink.profile | awk '{print $1+$2}' > ../store${dir}/${baseName}.${score}.${chr}.${i}
							cat summStat stat.$plim2 > ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
							gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
						fi
						let i=i+1
                                        done
                        	fi
                	fi
        	done
	done


elif [ $score == "annoPred" ]; then
        cp ~/prsDatabase/makeAnnoPredSet.R .
	numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr} | wc -l`
        if [ $numLine -lt 4 ]; then
	author=`echo $baseName | cut -f1 -d'.'`
	cp ~/prsDatabase/fixSSAnnoPred.R .
	cp ~/prsDatabase/metaData .
	cp ~/prsDatabase/annopredPhenoMaker.R .
	cp /home/cbsumezey/sdk2004/ldsplits/chr${chr}.bed .

	let otherDir=dir+30
	Rscript annopredPhenoMaker.R $chr $file
	ln -s ../new.${chr}.bed our.bed
	ln -s ../new.${chr}.bim our.bim

	Rscript fixSSAnnoPred.R

	sampSize=`cat metaData | grep $author | cut -f2`
	ln -s /home/cbsumezey/sdk2004/AnnoPred/ref ref

        maxLine=`cat chr${chr}.bed | wc -l`	
	rm compBeta*
	i=1
	while [ $i -lt $maxLine ]; do
		rm test* sumstats* tier*
		cat annoSummStat | head -1 > win.ss
		lengthSS=`cat win.ss | wc -l`
		start=`cat chr${chr}.bed | head -$i | tail -1 | cut -f2`
		while [ $lengthSS -le 2500 ] && [ $i -lt $maxLine ]; do
	                stop=`cat chr${chr}.bed | head -$i | tail -1 | cut -f3`
			cat annoSummStat | head -1 > win.ss
	                cat annoSummStat | awk -v var="$start" '$5 > var {print $0}' | awk -v var="$stop" '$5 < var {print $0}' >> win.ss
			lengthSS=`cat win.ss | wc -l`
			let i=i+1
		done
                cat win.ss | cut -f2 -d' ' > extraRsids
                plink --bfile our --chr $chr --extract extraRsids --make-bed --out win

		g=1
		pGuess=0.05
		if [ $g == 1 ]; then
			taskset -c $dir,$otherDir python /home/cbsumezey/sdk2004/AnnoPred/AnnoPred.py --sumstats=win.ss --ref_gt=win --val_gt=win \
				 --coord_out=testCOORD --N_sample=${sampSize} --annotation_flag="tier3" \
				 --P=$pGuess --local_ld_prefix=testLD --out=test --temp_dir=.
		else		
			taskset -c $dir,$otherDir python /home/cbsumezey/sdk2004/AnnoPred/AnnoPred.py --sumstats=win.ss --ref_gt=win --val_gt=win \
                	        --coord_out=testCOORD --N_sample=${sampSize} --annotation_flag="tier3" \
                                --P=$pGuess --user_h2=funcEnrich --local_ld_prefix=testLD --out=test --temp_dir=.
		fi
		cat test_h2_non_inf_betas*${pGuess}* | tail -n +2 >> compBeta1.$pGuess.$g
		cat test_h2_inf_betas*${pGuess}* | tail -n +2 >> compBeta2.$pGuess.$g
		cat test_pT_non_inf_betas*${pGuess}* | tail -n +2 >> compBeta3.$pGuess.$g
		cat test_pT_inf_betas*${pGuess}* | tail -n +2 >> compBeta4.$pGuess.$g
	done

        i=1
        ls compBeta* | while read line; do
		numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${i}.gz | wc -l`
                if [ $numLine -eq 0 ]; then
                        Rscript makeAnnoPredSet.R $line
	                plink --threads 1 --bfile ../new.${chr}  --score outSummStat 2 4 8 no-sum
	                if [ -f plink.profile ]; then
	                        sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
	                        gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
	                fi
	                mv outSummStat ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
		fi
		let i=i+1
        done
	fi

elif [ $score == "LDPredFunct" ]; then
        author=`echo $baseName | cut -f1 -d'.'`
        cp ~/prsDatabase/getFunctMat.R .
        cp ~/prsDatabase/metaData .
        cp ~/prsDatabase/pToStat .

        sampSize=`cat metaData | grep $author | cut -f2`
        cp ~/prsDatabase/ldpredfunctPhenoMaker.R .
        cp ~/prsDatabase/makeLDPredFunctBeta.R .
        Rscript ldpredfunctPhenoMaker.R $chr $file

        tc=1
        for ci in 1 2; do
		numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${tc} | wc -l`
        	if [ $numLine -eq 0 ]; then

                if [ $ci == 1 ] || [ $ci == 3 ]; then ta=1; else ta=2; fi
                cp ~/prsDatabase/totalAnnot${ta}.RDS totalAnnot.RDS
                cp ~/prsDatabase/allSplitHerit$ci allH2
                herit=`cat allH2 | grep $author | cut -f2 -d' '`
                cp ~/prsDatabase/coefRes2/${author}_coef${ci}.results coefs
                Rscript getFunctMat.R $herit $author
		let otherDir=dir+30

                for ldrad in 1000; do
                        rm outCoord
                        taskset -c $dir,$otherDir python ../LDpred-funct/ldpredfunct.py --gf=../new.${chr} --pf train.phen --FUNCT_FILE=funcEnrich --coord=outCoord \
                                --ssf=functSummStat --N=${sampSize} --posterior_means=outPost --H2=${herit} --out=outVal --ld_radius=${ldrad}

                        sed 's/  */\t/g' outPost_LDpred-inf-ldscore.txt | cut -f3 | tail -n +2 > rsids
                        fgrep -w -f rsids preSummStat > shortSummStat
                        sed 's/  */\t/g' outPost_LDpred-inf-ldscore.txt | tail -n +2 | sort -k3 > betterPost
                        cat shortSummStat | cut -f1-7 > keep1 ; cat shortSummStat | cut -f9-10 > keep2 ; cat betterPost | cut -f7 > betas
                        paste keep1 betas keep2 > summStat1

                        Rscript makeLDPredFunctBeta.R
                        cat finalBeta | cut -f1 > rsids
                        fgrep -w -f rsids preSummStat > shortSummStat
                        cat finalBeta | sort -k1 > betterPost
                        cat shortSummStat | cut -f1-7 > keep1 ; cat shortSummStat | cut -f9-10 > keep2 ; cat betterPost | cut -f2 > betas
                        paste keep1 betas keep2 > summStat2

                        for i in 1 2; do
				numLine=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.${tc} | wc -l`
                		if [ $numLine -eq 0 ]; then
                                plink --threads 1 --bfile ../new.${chr}  --score summStat$i 2 4 8 no-sum
                                if [ -f plink.profile ]; then
                                        sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${tc}
                                        gzip ../store${dir}/${baseName}.${score}.${chr}.${tc}
                                fi
                                mv summStat$i ../sets/chr${chr}/${baseName}.${score}.${chr}.${tc}
                                let tc=tc+1
				fi
                        done
                done
		fi
        done


elif [ $score == "sblup" ]; then  #WARNING - HIGH MEMORY USAGE
	checkFile=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.1 | wc -l`
        if [ $checkFile -eq 0 ]; then
        cp ~/prsDatabase/makeMA.R .
        cp ~/prsDatabase/metaData .
        cp ~/prsDatabase/allSplitHerit1 allH2

        author=`echo $file | cut -f1 -d'.'`
        sampSize=`cat metaData | fgrep $author | cut -f2`
        numSnps=`cat metaData | fgrep $author | cut -f3`
        herit=`cat allH2 | fgrep $author | cut -f2 -d' '`
	if (( $(echo "$herit < 0" |bc -l) )); then herit=0.01; fi
        param=`echo "scale=4; $numSnps * (1/$herit)" | bc`

        cat preSummStat > win.ss
        cat win.ss | cut -f2 | fgrep -w -f ~/prsDatabase/hapmapSnps > extraRsids
        plink --memory 4000 --bfile ../new.$chr --chr $chr --extract extraRsids --make-bed --out win
        Rscript makeMA.R $sampSize
        gcta64 --bfile win --cojo-file ss.ma --cojo-sblup $param --cojo-wind 100 --thread-num 1 --out part
        cat part.sblup.cojo >> sblup.sblup.cojo

        cat sblup.sblup.cojo | cut -f1 | sort | uniq > rsids
        cat preSummStat | fgrep -w -f rsids | sort -k2 > summStat
        cat sblup.sblup.cojo | sort -k1 > sblup
        cat sblup | cut -f4 > newBeta
        cat summStat | cut -f1-7 > part1; cat summStat | cut -f9-10 > part2
        paste part1 newBeta part2 > finalSummStat

        plink --memory 4000 --threads 1 --bfile ../new.${chr}  --score finalSummStat 2 4 8 no-sum
        if [ -f plink.profile ]; then
                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.1
                gzip ../store${dir}/${baseName}.${score}.${chr}.1
        fi
        mv finalSummStat ../sets/chr${chr}/${baseName}.${score}.${chr}.1
        fi

elif [ $score == "sbayesr" ]; then
	checkFile=`cat storeRecord | fgrep ${baseName}.${score}.${chr}.1 | wc -l`
        if [ $checkFile -eq 0 ]; then
        cat ~/prsDatabase/maHeader > inputMA
        cat preSummStat | cut -f2,4,5,6,8,9,10,11 >> inputMA

        gctb --sbayes R --ldm /home/cbsumezey/sdk2004/ldms/ukbb.${chr}.ldm.sparse \
         --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary inputMA \
         --chain-length 10000 --burn-in 2000 --out-freq 100 --out sim_1

        cat sim_1.snpRes | sed 's/ \+/\t/g' | cut -f2-12 | tail -n +2 | sort -k2 > sbayesResult
        cat sbayesResult | cut -f2 > rsids
        cat preSummStat | fgrep -w -f rsids | sort -k2 > readySummStat
        cat sbayesResult | cut -f8 > adjBeta
        cat readySummStat | cut -f1-7 > part1; cat readySummStat | cut -f9-12 > part2
        paste part1 adjBeta part2 > summStat

        plink --memory 4000 --threads 1 --bfile ../new.${chr}  --score summStat 2 4 8 no-sum
        if [ -f plink.profile ]; then
                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.1
                gzip ../store${dir}/${baseName}.${score}.${chr}.1
        fi
        mv summStat ../sets/chr${chr}/${baseName}.${score}.${chr}.1
        fi

elif [ $score == "stackCT" ]; then
	author=`echo $file | cut -f1 -d'.'`
	checkFile=`cat storeRecord | fgrep ${baseName}.${score}.${chr} | wc -l`
        if [ $checkFile -lt 3 ]; then
	cp ~/prsDatabase/stackCT.R .
	rm subset*

	Rscript stackCT.R $author $chr
	for i in {1..3};do
		plink --memory 4000 --threads 1 --bfile ../new.${chr}  --score summStat$i 2 4 8 no-sum
	        if [ -f plink.profile ]; then
	                sed 's/ \+/\t/g' plink.profile | cut -f7 > ../store${dir}/${baseName}.${score}.${chr}.${i}
	                gzip ../store${dir}/${baseName}.${score}.${chr}.${i}
	        fi
	        mv summStat$i ../sets/chr${chr}/${baseName}.${score}.${chr}.${i}
	done

	fi
fi

rm *

cd ..
###################################################

echo done
echo $dir >> possDirs
echo $author $dir $score `date` >> dones.$chr
