addOn="3"
ls *gz | while read line; do
	IFS='.'
	read -r -a array <<< "$line"
	IFS='-'
	newLine="${array[0]}.${array[1]}.${addOn}.${array[2]}"
	mv $line $newLine
done
