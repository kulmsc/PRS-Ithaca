cd results
ls | while read line; do
	baseName="${line::-3}"
	mv $line /home/cbsumezey/sdk2004/results/${baseName}.1.gz
done
cd ..
