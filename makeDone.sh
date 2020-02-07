
rm temp
for i in {1..48};do
        ls store${i} | while read line; do
                echo $line >> temp
        done
done

#for chr in 1 2 3 4 5 6; do
#	x=`cat temp | cut -f1,2,3,4 -d'.' | fgrep ss.$chr.annoPred | sort | uniq | wc -l`
#	yes dones | head -$x > dones.$chr
#done
