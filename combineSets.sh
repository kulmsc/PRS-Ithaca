for i in {1..8}; do
	cd chr$i
	ls | while read line; do
		author=`echo $line | cut -f1 -d'.'`
		sets=`echo $line | cut -f4 -d'.'`
		cat $line >> ../${author}.${sets}.1.set
		echo $author
		echo $sets

	done
	cd ..
done
