import array
import numpy as np
import matplotlib.pyplot as plt
import rootpy.ROOT as ROOT
from peaks import peaks

p=peaks("total.root")
yDict=p.yDict
snDict=p.snDict

normal=[]
edge=[]
invert=[]
for sn in snDict["all"]:
    y=yDict["array1000",sn]
    rEdge=len(y)-1
    aMax=y.argmax()
    aMin=y.argmin()
    if aMax>0 and rEdge>aMin:
	if aMax<aMin:
	    normal.append(sn)
	else:
	    y2=y[1:rEdge]
	    aMax=y2.argmax()
	    aMin=y2.argmin()
	    if aMax < aMin:
		normal.append(sn)
	    else:
		invert.append(sn)
    else:
	normal.append(sn)
peak={}
pict={}
for sn in normal:
    pList=[]
    y=yDict["array1000",sn]
    start=y.argmax();
    right=y[y.argmax():y.argmin()]
    left=y[0:y.argmax()+1]
    end=len(right)+1

    y2=left
    while  len(y2)!=0:
	aMax=y2.argmax()
	r=aMax
	while r==y2.argmax():
	    r-=1
	    y2=left[0:r+1]
	    if len(y2)==0:
		break
	pict["left",sn,aMax]= r
	pList.append(aMax)

    y2=right
    r=0
    thresh=y.max()*.1
    while len(y2)!=0:
	if y2.max()<.1*thresh:
	    break
	y2=y2[y2.argmax():]
	aMax=start+y2.argmax()
	if aMax not in pList:
	    pList.append(aMax)
	    pict["right",sn,aMax]= r+start
	start=aMax
	r=0
	while 0==y2.argmax():
	    r+=1
	    y2=y2[1:]
	    if len(y2)==0:
		break
    if len(pList)>1:
	print "sn :"+str(sn)
	print pList
    peak[sn]=pList


#############
