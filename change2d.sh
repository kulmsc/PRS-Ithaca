
for i in {1..48};do
  cd store$i
  ls | fgrep winnersCurse-2d | while read line;do
    first=`echo $line | cut -f1 -d'.'`
    second=`echo $line | cut -f2 -d'.' | cut -f1 -d'-'`
    third=`echo $line | cut -f3,4,5,6,7 -d'.'`
    mv $line ${first}.${second}-winnersCurse2d.${third}
  done
  cd ..
done
