'''
subroutines for read/write
based on driver.py
'''

import os
import sys

def linestolist2(linesraw,linelist):
    delimiter = ''
    for k in range(len(linelist)):
        temp = linelist[k]
        temp4 = list(temp)
        for j in range(len(temp4)):
            if temp4[j]=='\n':
                del temp4[j]
        temp5 = delimiter.join(temp4)
        temp6 = temp5.split(' ')
        linelist[k] = temp6
        
def readlog(logfile):
    f=open(logfile,'r')
    loglines=f.readlines()
    loglist=list(loglines)
    linestolist2(loglines,loglist)
    f.close()
    loglist2=[]
    for k in range(len(loglist)):
        temp=[]
        for j in range(len(loglist[k])):
                if loglist[k][j]!='':
                    temp.append(loglist[k][j])
        loglist2.append(temp)
    return loglist2

def writelog(logfile,listtowrite):
    f=open(logfile,'w')
    for item in listtowrite:
        for k in range(len(item)):
            f.write("%s " % item[k])
        f.write('\n')
    f.close()



