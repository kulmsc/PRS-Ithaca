score=$3
chr=$2
dir=$1

infoLim=0.3
fixAmbig=true
splitLimit=1000
deciLength=500

cd dir$dir

if [ -f allRsid ]; then
	rm allRsid 
	rm allPos
fi

cp ~/prsDatabase/subsetCLUMP .
cp ~/prsDatabase/keepEID .
plim=`cat subsetCLUMP | cut -f1 -d'-' | sort -n | tail -1`

#Change here to impact INFO>0.3 and non ambigous strands
ls *gz | while read file; do
	baseName="${file::-3}"
	zcat $file | awk -v var="$infoLim" '$7 > var {print $0}' > temp ; mv temp $baseName
	if [ "$fixAmbig" = "true" ] ; then
		cat $baseName | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' > temp ; mv temp $baseName
	fi
	if [ "$score" = "noCorrect" ] || [ "$score" = "thin" ] || [ "$score" = "prune" ]  || [ "$score" = "clump" ] ; then
		cat $baseName | awk -v var="$plim" '$10 < var {print $0}' > temp ; mv temp $baseName
	fi
	cat $baseName | cut -f3 >> allPos
	cat $baseName | cut -f2 >> allRsid
	gzip -f $baseName
done


cat allRsid | sort | uniq > temp; mv temp allRsid #could here split up allRsid and run in loops
cat allPos | sort -n | uniq > temp; mv temp allPos

bgenix -g ../chr${chr}.bgen -incl-rsids allRsid > new.bgen 2> bgenLog
plink2 --bgen new.bgen --sample ../chr${chr}.sample --make-bed --out new &> plink2Log
rm new.bgen
cat new.bim | cut -f2 | uniq -d > dupIDs



if [ "$score" = "noCorrect" ] || [ "$score" = "thin" ] || [ "$score" = "prune" ] ; then
	if [ "$score" = "noCorrect" -o "$score" = "clump" ]; then
		#No Correction
		plink --bfile new --exclude dupIDs --keep-fam keepEID --make-bed --out new2 &> plink1Log
	elif [ "$score" = "thin" ]; then
		#Arbitrary thinning
		plink --bfile new --exclude dupIDs --thin 0.5 --keep-fam keepEID --make-bed --out new2 &> plink1Log 
	elif [ "$score" = "prune" ]; then
		#LD based pruning
		plink --bfile new --exclude dupIDs --indep 50 5 1.5 --keep-fam keepEID --make-bed --out new2 &> plink1Log
	fi

	ls *split.*.gz | while read summStat; do
		baseName="${summStat::-3}"
		zcat $summStat > summStat
		plink --silent --memory 5000 --bfile new2 --score summStat 2 4 8 no-sum --q-score-range subsetPVAL summStat 2 10 header --out new &> plink2Log
	
		for i in {1..3}; do
			if [ -f new.P${i}.profile ]; then
				sed 's/ \+/\t/g' new.P${i}.profile | cut -f7 > ../store${dir}/${baseName}.${i}.P
			        gzip ../store${dir}/${baseName}.${i}.P
			fi
		done
	done


elif [ "$score" = "clump" ] ; then
	cp ~/prsDatabase/fileHeader .
	cp ~/prsDatabase/decimate.py .
	cp ~/prsDatabase/fixSplit.R .
	ls *split.*.gz | while read summStat; do
		echo START BIG CHAIN WITH $summStat
                baseName="${summStat::-3}"
		zcat $summStat > totalSummStat
		cat fileHeader totalSummStat > headTotalSummStat
		ssLength=$(cat totalSummStat | wc -l)	
		if [ $ssLength -ge $splitLimit ]; then
			python decimate.py $summStat $deciLength shorter
			let newDeciLength=deciLength*2
			python decimate.py $summStat $newDeciLength longer
		fi
		i=1
		cat subsetCLUMP | cut -f1 -d'-' | while read sub; do
			pLimit=`echo $sub | cut -f1 -d-`
			r2Val=`echo $sub | cut -f2 -d-`
			#j=1
			#cat subsetCLUMP | cut -f2 -d'-' | while read r2Val; do
				echo DOING EACH P AND R
				rm totalRsids
				changeLim=`echo $pLimit*$r2Val | bc -l`
				if [[ $changeLim < 0.001 ]]; then
					addOn=longer
				else
					addOn=shorter
				fi

				if [ $ssLength -ge $splitLimit ]; then
					ls ${baseName}.[0-9].${addOn}.* | while read subSummStat; do
						echo HAS FOUND THE SUBSUMM $subSummStat
		                	        zcat $subSummStat > summStat
                        			cat fileHeader summStat > temp ; mv temp summStat
						plink --memory 5000 --bfile new --exclude dupIDs --keep-fam keepEID --clump summStat --clump-p1 $pLimit --clump-r2 $r2Val &> plink2Log
						if [ -f plink.clumped ]; then
							echo WE HAVE CLUMPED
							sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 >> totalRsids
							rm plink.clumped
						fi
					done
					if [ -f totalRsids ]; then
						echo $summStat GOING TO FIXSPLIT.R
						plink --bfile new --exclude dupIDs --r2 --ld-snp-list totalRsids --ld-window-kb 250 --ld-window 1000 --ld-window-r2 $r2Val &> plink2Log
						Rscript fixSplit.R
						ls done*
						echo DONE FIXING
					fi
				else
					plink --memory 5000 --bfile new --exclude dupIDs --keep-fam keepEID --clump headTotalSummStat --clump-p1 $pLimit --clump-r2 $r2Val &> plink2Log
					if [ -f plink.clumped ]; then
						sed -e 's/ [ ]*/\t/g' plink.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > doneRsids
                                                rm plink.clumped
					fi				
				fi			
				if [ -f doneRsids ]; then
					grep -f doneRsids totalSummStat > doneSummStat
					rm doneRsids
					plink --silent --memory 5000 --bfile new --exclude dupIDs --keep-fam keepEID --score doneSummStat 2 4 8 no-sum --out new &> plink2Log
					if [ -f new.profile ]; then
						sed 's/ \+/\t/g' new.profile | cut -f7 > ../store${dir}/${baseName}.${i}.P
						gzip ../store${dir}/${baseName}.${i}.P
						rm new.profile
					fi
				fi
				#let j++
			#done
			let i++
		done
	done



elif [ "$score" = "ldpred" ] ; then
	cp /home/cbsumezey/sdk2004/1000genomes/eur.${chr}.*  .
	cp /home/sdk2004/prsDatabase/makeLDPredStandard.py .
	cp /home/sdk2004/prsDatabase/metaData .
	cp /home/sdk2004/prsDatabase/subsetRho .
	ls *.gz | cut -f1 -d'.' | while read authorName; do
		wideName=${authorName}.ss.${chr}.wide.${dir}.gz
                splitName=${authorName}.ss.${chr}.split.${dir}.gz
                baseName="${stdName::-3}"
                python makeLDPredStandard.py $wideName
		snps=`cat summStat | wc -l`
		let ldr=snps/3000
		coord --gf=eur.$chr --ssf=summStat --N=`cat sampleSize` --out=madeCoord
		
		#NOTES!!!! That the output of ldpred will not simply be summStat.adj and should figure out how to use --local_ld_file_prefix
		cat subsetRho | cut -f2 | while read rho; do
			ldpred --cord=madeCoord --ld_radius=$ldr --PS=$rho --N=`cat sampleSize` --out summStat.adj
			plink --silent --memory 5000 --bfile new --exclude dupIDs --keep-fam keepEID --score summStat.adj 3 4 6 header --out new
			if [ -f new.profile ]; then
                                sed 's/ \+/\t/g' new.profile | cut -f7 > ../store${dir}/${baseName}.${i}.P
                                gzip ../store${dir}/${baseName}.${i}.P
                        	rm new.profile
			fi
		done


	done
fi

cd ..
