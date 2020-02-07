import numpy as np
import pdb
import subprocess as sb
import gzip
import os
import sys
import argparse

ssName="preSummStat"

parser = argparse.ArgumentParser(description='A tutorial of argparse!')
parser.add_argument("--N", default=1, help="N")

args = parser.parse_args()
N = args.N

def normRead(fileName,delimiter):
        totalData=[]
        with open(fileName,"r") as f:
                for line in f.read().splitlines():
                        totalData.append(line.split())
        totalData=np.array(totalData)
        return(totalData)

def justRead(fileName,delimiter):
	rList=[]
        with open(fileName,"r") as f:
                for line in f.read().splitlines():
			rList.append(line.split(delimiter))
        return(np.array(rList))

def justWrite(theStuff,fileName,delimiter):
        with open(fileName, 'wb') as f:
                for line in theStuff:
                        f.write(delimiter.join(line)+'\n')

ss=normRead(ssName,'\t')
for i in range(ss.shape[0]):
	ss[i,0]='chr'+ss[i,0]

#STANDARD
std=ss[:,(0,2,4,3,5,8,1,9,7)]
std=np.insert(std,0,['chr','pos','ref','alt','reffrq','info','rs','pval','effalt'],axis=0)

#CUSTOM
#std=ss[:,(0,2,4,3,1,7,9)]
#std[:,5] = np.exp(std[:,5].astype("float64")).astype("str")
#std=np.insert(std, 0, ['hg19chrc', 'pos', 'a1', 'a2', 'snpid', 'or', 'p'], axis=0)

#LDPRED
#std=ss[:,(0,2,1,4,3,5,9,7,8)]
#std=np.hstack((std, np.tile(N, (std.shape[0], 1))))
#std=np.insert(std,0,['CHR','POS','SNP_ID','REF','ALT','REF_FRQ','PVAL','BETA','SE','N'],axis=0)

justWrite(std,'summStat','\t')

#meta=justRead('metaData','\t')
#author=ssName.split('.')[0]
#n=meta[meta[:,0]==author,1]
#with open('sampleSize','w') as f:
#	f.write(str(n))
