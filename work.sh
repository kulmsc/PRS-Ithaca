echo $1

echo before
cat possDirs
#DETERMINE THE DIRECTORY #######################
toOpen=`cat possAccess`
keepGoing=True

while [ $keepGoing == "True" ]; do
        if [ $toOpen == "True" ]; then
		echo going for DIR
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
echo we are in directory $dir
##################################################

cat possDirs
echo after
sleep 8
echo $dir >> possDirs
