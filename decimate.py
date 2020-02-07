import numpy as np
import subprocess as sb
import gzip
import os
import sys

fileName=sys.argv[1]
splitName=fileName.split('.')

splitBy=sys.argv[2]
splitBy=int(splitBy)

addOn=sys.argv[3]

def gzRead(fileName,delimiter):
        totalData=[]
        with gzip.open(fileName,"r") as f:
                for line in f.read().splitlines():
                        totalData.append(line.split(delimiter))
        totalData=np.array(totalData)
        return(totalData)

def gzWrite(theStuff,fileName,delimiter):
        with gzip.open(fileName, 'wb') as f:
                for line in theStuff:
                        f.write(delimiter.join(line)+'\n')

ss=gzRead(fileName,'\t')
fileBits=ss.shape[0]/splitBy
if fileBits > 0:
	start=0
	stop=splitBy
	for i in range(fileBits):
		newName=list(splitName)
		newName.insert(5,str(i+1)+'.'+addOn)
		newName='.'.join(newName)
		subSS=ss[start:stop,:]
		gzWrite(subSS,newName,'\t')
		start+=splitBy
		stop+=splitBy
		if i+2==fileBits:
			stop=ss.shape[0]
else:
	newName=list(splitName)
	newName.insert(5,'1'+addOn)
	newName='.'.join(newName)
	gzWrite(ss,newName,'\t')
