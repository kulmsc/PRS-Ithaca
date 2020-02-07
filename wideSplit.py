import numpy as np
import subprocess as sb
import gzip
import os

np.set_printoptions(suppress=True,formatter={'float_kind':'{:0.2f}'.format})


def gzRead(fileName,delimiter):
        totalData=[]
        with gzip.open(fileName,"r") as f:
                for line in f.read().splitlines():
                        totalData.append(line.split())
        totalData=np.array(totalData)
        return(totalData)

def listRead(fileName):
        with open(fileName,"r") as f:
                rList=f.read().splitlines()
        return(rList)

def gzWrite(theStuff,fileName,delimiter):
	with gzip.open(fileName, 'wb') as f:
		for line in theStuff:
			f.write(delimiter.join(line)+'\n')


sb.call("rm *split*",shell=True)
sb.call("rm allPos",shell=True)
sb.call("ls *gz | while read line ; do zcat $line | cut -f3 >> allPos ; done",shell=True)
allPos=list(set(listRead("allPos")))
allPos=list(map(int,allPos))
allPos.sort()
sb.call("rm allPos",shell=True)
interval=len(allPos)/3
breakPos=[allPos[interval*i] for i in range(3)[1:]]
breakPos.append(999999999)

avaFiles=[f for f in os.listdir('.') if f.endswith('.gz')]

for ava in avaFiles:
	print('ava',ava)
	beginName=ava.split("gz")[0]
	readSS=gzRead(ava,"\t")
	readComp=readSS[:,2].astype('int')

	i=1
	oldPosDelim=0
	for posDelim in breakPos:
		newName=beginName+'split.'+str(i)+'.gz'
		newNameWide=beginName+'wide.'+str(i)+'.gz'
		subSS=readSS[np.logical_and(readComp<=posDelim,readComp>oldPosDelim),:]
		subSSWide=readSS[np.logical_and(readComp<=(posDelim+250000),readComp>(oldPosDelim-250000)),:]
		oldPosDelim=posDelim
		gzWrite(subSS.tolist(),newName,'\t')
		gzWrite(subSSWide.tolist(),newNameWide,'\t')
		i=i+1

