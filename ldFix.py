import numpy as np
import subprocess as sb
import gzip
import os
import sys

fName=sys.argv[1]

def gzRead(fileName,delimiter):
        totalData=[]
        with gzip.open(fileName,"r") as f:
                for line in f.read().splitlines():
                        totalData.append(line.split(delimiter))
        totalData=np.array(totalData)
        return(totalData)

def normRead(fileName,delimiter):
        totalData=[]
        with open(fileName,"r") as f:
                for line in f.read().splitlines():
        		splitLine=line.split(delimiter)
			if len(splitLine)==7: 
			       	totalData.append(splitLine)
	del totalData[0]
        totalData=np.array(totalData)
        return(totalData)

def listRead(fileName):
        with open(fileName,"r") as f:
                rList=f.read().splitlines()
        return(rList)

#read in the list of all poss rsids

possRsid=listRead("totalRsids")
ld=normRead("res.ld","\t")
ss=gzRead(fName,"\t")
toRemove=[]
for i in range(ld.shape[0]):
	snp1=ld[i,2]
	snp2=ld[i,5]
	if snp1 != snp2 and snp1 in possRsid and snp2 in possRsid:
		pos1=np.where(ss[:,1]==snp1)[0][0]
		pos2=np.where(ss[:,1]==snp2)[0][0]
		if ss[pos1,9] < ss[pos2,9]:
			index=possRsid.index(snp1)
			del possRsid[index]
		else:
			index=possRsid.index(snp2)
			del possRsid[index]

with open("doneRsids","w") as f:
	for line in possRsid:
		f.write(line+'\n')	
#now extract the good, or pruned summStat
#write pruned summStat

