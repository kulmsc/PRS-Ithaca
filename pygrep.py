import numpy as np

def normRead(fileName,delimiter):
        totalData=[]
        with open(fileName,"r") as f:
                for line in f.read().splitlines():
                        totalData.append(line.split())
        totalData=np.array(totalData)
        return(totalData)

def justWrite(theStuff,fileName,delimiter):
        with open(fileName, 'wb') as f:
                for line in theStuff:
                        f.write(delimiter.join(line)+'\n')

ss=normRead('summStat','\t')
ssTop=ss[0,:]
rsid=normRead('doneRsids','\t')
ss2=ss[np.isin(ss[:,1],rsid[:,0]),:]
ss2=np.vstack((ssTop,ss2))

justWrite(ss2,"doneSummStat",'\t')
