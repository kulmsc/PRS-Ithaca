phaseNum=$1

cp /home/sdk2004/prsDatabase/phase${phaseNum}.eid phase.eid

maxJobs=`cat maxJobs`
penalty=0
chrStart=15
chrStop=22

if [ $# -eq 0 ]; then
        echo "no argument"
else
	echo "have argument" $1
	rm possDirs
	for (( i=1; i<=$maxJobs; i++ )); do
		echo $i >> possDirs
	done
	echo True > possAccess

	#iterating through all of the chromosomes
	for (( num=$chrStart; num<=$chrStop; num++ )); do
		echo chromosome $num
		cat allScores | while read score; do		
			echo score $score
			
			#iterating through all of the summary statistic files
			ls /home/cbsumezey/sdk2004/chr${num} | while read line; do
				echo the file is $line

				author=`echo $line | cut -f1 -d'.'`

				echo onto scoring
				./totalScore.sh $num $line $score &> logs/${line}.${num}.${score}.log &
				sleep $(( ( RANDOM % 5 )  + 1 ))

				goOn=False
				while [ $goOn == "False" ]; do
					openSlots=`cat possDirs | wc -l`
					newMaxJobs=`cat maxJobs`

					#checking to see if the maxJobs file has changed
					if [ $newMaxJobs -gt $maxJobs ]; then
						let diff=$newMaxJobs-$maxJobs
						for (( i=1; i<=$diff; i++ )); do
							let newDir=$i+$maxJobs
                					echo $newDir >> possDirs
        					done
						maxJobs=$newMaxJobs
					elif [ $newMaxJobs -lt $maxJobs ]; then
						let penalty=$maxJobs-$newMaxJobs
						maxJobs=$newMaxJobs
					fi

	
					if [ $openSlots -gt 0 ]; then
						if [ $penalty -gt 0 ]; then
							let penalty=$penalty-1
							keepGoing=True
							toOpen=`cat possAccess`
							while [ $keepGoing == "True" ]; do
								if [ $toOpen == "True" ]; then
          								echo False > possAccess
                							dir=`head -1 possDirs`
                							grep -v -w $dir possDirs > tempPoss; mv tempPoss  possDirs
                							echo True > possAccess
                							keepGoing=False
        							else
                							toOpen=`cat possAccess`
                							sleep 1
        							fi
							done
						else
							echo NOW WE CAN GO
							goOn=True
						fi
					else
						echo MUST WAIT FOR ROOM TO GO
						sleep $(( ( RANDOM % 8 )  + 1 ))
						openSlots=`cat possDirs | wc -l`
					fi
				done
			done		
		done
	done
fi
