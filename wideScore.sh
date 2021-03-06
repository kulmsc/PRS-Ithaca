score=$3
chr=$2
dir=$1

infoLim=0.3
fixAmbig=true

cd dir$dir

if [ -f allRsid ]; then
	rm allRsid 
	rm allPos
fi

cp ~/prsDatabase/subsetCLUMP .
cp ~/prsDatabase/keepEID .
plim=`cat subsetCLUMP | tail -1 | cut -f1 -d'-'`

#Change here to impact INFO>0.3 and non ambigous strands
ls *.gz | grep wide | cut -f1 -d'.' | while read authorName; do
	wideName=${authorName}.ss.${chr}.wide.${dir}.gz
	splitName=${authorName}.ss.${chr}.split.${dir}.gz
	wideBase="${wideName::-3}"
	splitBase="${splitName::-3}"
	zcat $wideName | awk -v var="$infoLim" '$7 > var {print $0}' > temp ; mv temp $wideBase
	zcat $splitName | awk -v var="$infoLim" '$7 > var {print $0}' > temp ; mv temp $splitBase
	if [ "$fixAmbig" = "true" ] ; then
		cat $wideBase | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' > temp ; mv temp $wideBase
		cat $splitBase | awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' > temp ; mv temp $splitBase
	fi
	if [ "$score" = "noCorrect" ] || [ "$score" = "thin" ] || [ "$score" = "prune" ]  || [ "$score" = "clump" ] ; then
		cat $wideBase | awk -v var="$plim" '$10 < var {print $0}' > temp ; mv temp $wideBase
		cat $splitBase | awk -v var="$plim" '$10 < var {print $0}' > temp ; mv temp $splitBase
	fi
	cat $wideBase | cut -f2 >> allRsid
	cat $wideBase | cut -f3 >> allPos
	gzip -f $wideBase
	gzip -f $splitBase
done


cat allRsid | sort | uniq > temp; mv temp allRsid #could here split up allRsid and run in loops
cat allPos | sort -n | uniq > temp; mv temp allPos

bgenix -g ../chr${chr}.bgen -incl-rsids allRsid > new.bgen 2> bgenLog
plink2 --bgen new.bgen --sample ../chr${chr}.sample --make-bed --out new &> plink2Log
rm new.bgen
cat new.bim | cut -f2 | uniq -d > dupIDs

start=$(head -1 allPos)
stop=$(tail -1 allPos)
bgenix -g ../chr${chr}.bgen -incl-range ${chr}:${start}-${stop} > full.bgen 2> bgenLog
plink2 --bgen full.bgen --sample ../chr${chr}.sample --make-bed --out full
rm full.bgen
cat full.bim | cut -f2 | uniq -d > fullDups



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
	ls *.gz | grep wide | cut -f1 -d'.' | while read authorName; do
		wideName=${authorName}.ss.${chr}.wide.${dir}.gz
		splitName=${authorName}.ss.${chr}.split.${dir}.gz
		baseName="${wideName::-3}"
		zcat $wideName > summStat
		cat fileHeader summStat > temp ; mv temp summStat
		zcat $splitName | cut -f2 > origRsids
		i=1
		cat subsetCLUMP | cut -f1 -d'-' | while read pLimit; do
			j=1
			cat subsetCLUMP | cut -f2 -d'-' | while read r2Val; do
				plink --memory 5000 --bfile full --exclude fullDups --keep-fam keepEID --clump summStat --clump-p1 $pLimit --clump-r2 $r2Val
				if [ -f plink.clumped ]; then
					sed "s/^[ \t]*//" plink.clumped | sed 's/ \+ /\t/g' | sed '/^\s*$/d' | cut -f3 > rsidsEx
					grep -f rsidsEx summStat | grep -f origRsids > clump.summStat
					plink --silent --memory 5000 --bfile new --exclude dupIDs --keep-fam keepEID --score clump.summStat 2 4 8 no-sum --out new
					if [ -f new.profile ]; then
						sed 's/ \+/\t/g' new.profile | cut -f7 > ../store${dir}/${baseName}.${i}.${j}.P
						gzip ../store${dir}/${baseName}.${i}.${j}.P
						rm new.profile
					fi
					rm plink.clumped
				fi
				let j++
			done
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
