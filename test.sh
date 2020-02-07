rm possDirs
for i in {1..18}; do
        echo $i >> possDirs
done
echo True > possAccess



maxJobs=3
num=1
ls /home/cbsumezey/sdk2004/chr${num} | grep gz$ | while read line; do
	echo before
	./work.sh $line	&
	echo after
	sleep 1

	echo $! >> currentJobs
	ps S > currentPs
	numJobs=`grep -f currentJobs currentPs | wc -l`
	goOn=False
	while [ $goOn == "False" ]; do
		if [ $numJobs -lt $maxJobs ]; then
			echo NOW WE CAN GO
			goOn=True
		else
			echo MUST WAIT FOR ROOM TO GO
			sleep 5
			ps S > currentPs
			numJobs=`grep -f currentJobs currentPs | wc -l`
		fi
	done		
done

