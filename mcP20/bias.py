import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import gROOT
from math import *
from array import array
from scipy import stats
import random

def findEdge(n):
    # print "finding Edge"
    N=len(n)
    sum=0
    count=0
    for i in xrange(0,N):
	if n[i]!=0:
	    sum+=float(n[i])
	    count+=1.0
    truAve= sum/count
    start=0
    end=0
    print "number of events: "+str(sum)
    print "number of events per bin "+str(truAve)
    for i in xrange(0,N):
	t= n[i]+3*sqrt(n[i])
	start=i
	if t>truAve:
	    break
    
    for i in xrange(N-1,start,-1):
	t= n[i]+3*sqrt(n[i])
	end=i
	if t>truAve:
	    break

    return [start,end]

def randList(n,size=None):
    if size is None:
	size=int(.1*n)
    ret=[]
    if size>n:
	print "size to big empyt list"
	return ret
    while len(ret)<size:
	a=random.randint(0,n-1)
	if a not in ret:
	    ret.append(a)
    return ret

# inFile=ROOT.TFile("proc_cry_0.root")
# inFile=ROOT.TFile("proc_cry_total.root")
inFile=ROOT.TFile("proc_cry_total9.root")
# inFile=ROOT.TFile("proc_cry_total3.root")
tree=ROOT.TTree()
tree=inFile.Get("photon_Data")
# y=np.zeros(1,dtype=float)
x=array("f",[0])
y=array("f",[0])
z=array("f",[0])
t=array("f",[0])
pid=array("i",[0])

tree.SetBranchAddress("x",x)
tree.SetBranchAddress("y",y)
tree.SetBranchAddress("z",z)
tree.SetBranchAddress("t",t)
tree.SetBranchAddress("pid",pid)
bias_est=[]
bias_hi=[]
bias_lo=[]
slopeList=[]
estList=[]
entries=int(tree.GetEntries())
px=[]
py=[]
pz=[]
tim=[]
cellLen=950
for i in xrange(0,entries):
    tree.GetEntry(i)
    # if abs(t[0])!=0:
    if abs(t[0])!=0 and abs(x[0])<(60.0) and abs(z[0])<(5.0) and abs(y[0])<(500.0)and abs(pid[0])==13:
	px.append(x[0])
	py.append(y[0])
	pz.append(z[0])
	tim.append(t[0])
print "number of muon events "+str(len(tim))
for q in xrange(0,1000):
    print q
    rl= randList(len(tim))
    samTime=[]
    samPos=[]
    for v in rl:
	samTime.append(tim[v])
	samPos.append(py[v])

    time=np.array(samTime)
    pos=np.array(samPos)
    slope, intercept, r_value, p_value, std_err = stats.linregress(time,pos)
# print "events: "+str(len(time))
    dy=sqrt(12*pos.var())
    dt=sqrt(12*time.var())
# print "b-a "+str(dt)
# print "dy "+str(dy)
    est=-cellLen/dt
# print "intercept "+str(intercept)
# print "r "+str(r_value)
# print "p "+str(p_value)
# print "std "+str(std_err)
############scatter###########
# plt.plot(time,pos,'.')
# plt.show()
#############time hist***********
    n, bins, patches=plt.hist(time,bins=61)
    # plt.clf()
    # plt.plot(time,pos,'.')
# print np.array(n).mean()
# print n
    ival= findEdge(n)
# print "ival "+str(ival)
# print n[ival[0]]
# print n[ival[1]]
    dt_dumb_lo= bins[ival[1]]-bins[ival[0]]
    if dt_dumb_lo!=0:
	slope_dumb_lo=-cellLen/dt_dumb_lo
    else:
	continue 
# print "dt dumb lo "+str(dt_dumb_lo)
# print "slope dume lo "+str(slope_dumb_lo)

    dt_dumb_hi= bins[ival[1]-1]-bins[ival[0]+1]
    if dt_dumb_hi!=0:
	slope_dumb_hi=-cellLen/dt_dumb_hi
    else:
	continue 
    # print "slope "+str(slope)
    # print "slope est "+str(est)
# print "dt dumb hi "+str(dt_dumb_hi)
# print "slope dume hi "+str(slope_dumb_hi)

    best=100*(slope-est)/slope
    blo=100*(slope-slope_dumb_lo)/slope
    bhi=100*(slope-slope_dumb_hi)/slope
    if (abs(blo)>200 or abs(bhi) >200 or abs(best) >200) and False :
	print "slope "+ str(slope)
	print "est "+ str(est)
	print "lo slope "+ str(slope_dumb_lo)
	print "hi slope "+ str(slope_dumb_hi)
	print "lo bias "+ str(blo)
	print "hi bias "+ str(bhi)
	print "est bias"+ str(best)
	print n
	print ival
	print "first bin "+str(n[ival[0]])
	print "second bin "+str(n[ival[1]])
	print "time0 lo "+str(bins[ival[0]])
	print "time1 lo "+str(bins[ival[1]])
	print "time0 hi  "+str(bins[ival[0]+1])
	print "time1 hi "+str(bins[ival[1]-1])
	# plt.show()
	a=1
	while(a>0 and False):
	    lo=input("enter first bin")
	    hi=input("enter second bin")
	    dt_dumb_lo= bins[hi]-bins[lo]
	    dt_dumb_hi= bins[hi-1]-bins[lo+1]
	    slope_dumb_lo=-cellLen/dt_dumb_lo
	    slope_dumb_hi=-cellLen/dt_dumb_hi
	    blo=100*(slope-slope_dumb_lo)/slope
	    bhi=100*(slope-slope_dumb_hi)/slope
	    print "lo "+ str(blo)
	    print "hi "+ str(bhi)
	    print "est "+ str(best)
	    a=input("enter neg value to cont")
    bias_est.append(best)
    bias_lo.append(blo)
    bias_hi.append(bhi)
    slopeList.append(slope)
    estList.append(est)
    plt.clf()
print "est mean " + str(np.array(bias_est).mean())
print "est std " +str(np.array(bias_est).std())
print "lo mean " +str(np.array(bias_lo).mean())
print "lo std " +str(np.array(bias_lo).std())
print "hi mean " +str(np.array(bias_hi).mean())
print "hi  std" +str(np.array(bias_hi).std())
print "slope mean " +str(np.array(slopeList).mean())
print "slope  std" +str(np.array(slopeList).std())
print "est mean " +str(np.array(estList).mean())
print "est  std" +str(np.array(estList).std())
print "sample size "+str(.1*len(time))
plt.clf()
plt.hist(bias_lo)
plt.title("long dt est")
plt.show()
plt.clf()
plt.hist(bias_hi)
plt.title("short dt est")
plt.show()
plt.clf()
plt.hist(bias_est)
plt.title("Variance est dt")
plt.show()
plt.hist(slopeList)
plt.title("Slope")
plt.show()
plt.clf()
plt.hist(estList)
plt.title("Estimated Slope")
plt.show()
