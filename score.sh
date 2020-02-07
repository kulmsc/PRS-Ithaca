score=$3
chr=$2
dir=$1

cd dir$dir

if [ -f allRsid ]; then
	rm allRsid 
fi

cp ~/prsDatabase/subsetPVAL .
cp ~/prsDatabase/keepEID .
plim=`cat subsetPVAL | tail -1 | cut -f3`

#Change here to impact INFO>0.3 and non ambigous strands
ls *split.*.gz | while read file; do
	baseName="${file::-3}"
	zcat $file | awk '$7 > 0.3 {print $0}' | awk -v var="$plim" '$10 < var {print $0}' | \
	awk '{ if (!($4=="A" && $5=="T" || $4 =="T" && $5 =="A" || $4 =="C" && $5 =="G" || $4 =="G" && $5 =="C") ) print $0; }' > temp ; mv temp $baseName ; gzip -f $baseName
	zcat $file | cut -f2 >> allRsid
done
cat allRsid | sort | uniq > temp; mv temp allRsid #could here split up allRsid and run in loops


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
else
	cp ~/prsDatabase/fileHeader .
	ls *split.*.gz | while read summStat; do
		baseName="${summStat::-3}"
		zcat $summStat > summStat
		cat fileHeader summStat > temp ; mv temp summStat
		i=1
		cat subsetPVAL | cut -f3 | while read pLimit; do	
			plink --memory 5000 --bfile new --exclude dupIDs --keep-fam keepEID --clump summStat --clump-p1 $pLimit --clump-r2 0.5
			if [ -f plink.clumped ]; then
				sed "s/^[ \t]*//" plink.clumped | sed 's/ \+ /\t/g' | sed '/^\s*$/d' | cut -f3 > rsidsEx
				grep -f rsidsEx summStat > clump.summStat
				plink --silent --memory 5000 --bfile new --exclude dupIDs --keep-fam keepEID --score clump.summStat 2 4 8 header no-sum --out new
				if [ -f new.profile ]; then
					sed 's/ \+/\t/g' new.profile | cut -f7 > ../store${dir}/${baseName}.${i}.P
					gzip ../store${dir}/${baseName}.${i}.P
					rm new.profile
				fi
				rm plink.clumped
			fi
			let i++
		done
	done
fi

cd ..
