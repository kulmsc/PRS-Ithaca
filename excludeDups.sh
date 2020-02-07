for i in {1..28};do
	echo $i
	ls store${i}/* | cut -f2 -d'/' > temp
	for j in {1..28};do
		if [ $i != $j ]; then
			ls store${j}/* | cut -f2 -d'/' | fgrep -w -f temp | while read line; do
				rm store${j}/${line}
			done
		fi
	done
done
