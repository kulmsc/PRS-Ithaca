
for chr in {1..22}; do
	cd chr$chr
		for i in {4..24}; do
			ls | fgrep lassosum | fgrep -v -f ../../ldpredfunctNames | grep $i$ | while read line; do
				rm $line
			done
		done
	cd ..
done

#for i in {1..48}; do
#        if [ `ls store${i} | fgrep "$1" | wc -l` -gt 0 ]; then
#                echo store$i
#       fi
#done


#cat dupFiles | while read line; do
#        echo $line
#        for i in {1..48}; do
#                if [ -e store${i}/"$line" ]; then
#                        echo store$i
#                        rm store${i}/$line
#                        break
#                fi
#        done
#done

#cat moveLines | while read line; do
#        echo $line
#
#	rm theStores
#	for i in {1..48}; do
#        	if [ `ls store${i} | fgrep "$line" | wc -l` -gt 0 ]; then
#                	echo store$i >> theStores
#        	fi
#	done
#
#	fromStore=`head -1 theStores`
#	toStore=`tail -1 theStores`
#	#cat theStores
#	echo "mv ${fromStore}/${line}* $toStore"
#	mv ${fromStore}/${line}* $toStore
#done
