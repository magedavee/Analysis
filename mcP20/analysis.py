import numpy as np
import math
from math import *
import matplotlib.pyplot as plt
import ROOT
from ROOT import gROOT
from math import *
from array import array
from scipy import stats
from scipy.stats import norm
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import Axes3D as ax3
inFile=ROOT.TFile("proc_cry_total12.root")
# inFile=ROOT.TFile("proc_cry_total9.root")
# inFile=ROOT.TFile("cry_0.root.pro")
tree=ROOT.TTree()
tree=inFile.Get("photon_Data")
# y=np.zeros(1,dtype=float)
x=array("f",[0])
y=array("f",[0])
z=array("f",[0])
x0=array("f",[0])
y0=array("f",[0])
z0=array("f",[0])
x1=array("f",[0])
y1=array("f",[0])
z1=array("f",[0])
left=array("i",[0])
right=array("i",[0])
t=array("f",[0])
e=array("f",[0])
pid=array("i",[0])
vol=array("i",[0])

tree.SetBranchAddress("x",x)
tree.SetBranchAddress("y",y)
tree.SetBranchAddress("z",z)
tree.SetBranchAddress("x0",x0)
tree.SetBranchAddress("y0",y0)
tree.SetBranchAddress("z0",z0)
tree.SetBranchAddress("x1",x1)
tree.SetBranchAddress("y1",y1)
tree.SetBranchAddress("z1",z1)
tree.SetBranchAddress("t",t)
tree.SetBranchAddress("left",left)
tree.SetBranchAddress("right",right)
tree.SetBranchAddress("pid",pid)
tree.SetBranchAddress("vol",vol)
tree.SetBranchAddress("e",e)
entries=int(tree.GetEntries())
# entries=26255
px=[]
py=[]
pz=[]
px0=[]
py0=[]
pz0=[]
px1=[]
py1=[]
pz1=[]
dx=[]
dy=[]
dz=[]
leftList=[]
rightList=[]
tim=[]
volList=[]
cellLen=950
eList=[]
rList=[]
slopeList=[]
zCutList=[]
ratList=[]
totalList=[]
for i in xrange(0,entries):
    tree.GetEntry(i)
    # if abs(t[0])!=0 and abs(pid[0])==13:
    if abs(t[0])<30 and abs(t[0])!=0 and abs(x[0])<(80) and abs(z[0])<(80) and abs(y[0])<(500.0) and not (math.isnan(x0[0]) or math.isnan(x1[0])) and abs(pid[0])==13 :
	    if (z0[0]-z1[0])>0 and right[0]>0 and left[0]>0 :
		# if y[0]>(-57*t[0]-100) and y[0]<(-57*t[0]+100): 
		# if y[0]<(-57*t[0]-100) or y[0]>(-57*t[0]+100):
                leftList.append(left[0])
                rightList.append(right[0])
                rat=float(left[0])/float(right[0])
                total=left[0]+right[0]
                totalList.append(total)
                lr=log10(rat)
                ratList.append(lr)
                px.append(x[0])
                py.append(y[0])
                pz.append(z[0])
                px0.append(x0[0])
                py0.append(y0[0])
                pz0.append(z0[0])
                px1.append(x1[0])
                py1.append(y1[0])
                pz1.append(z1[0])
                dx.append(x0[0]-x1[0])
                dy.append(y0[0]-y1[0])
                dz.append(z0[0]-z1[0])
                tim.append(t[0])
                volList.append(vol[0])
                eList.append(e[0])
print "number of muon events "+str(len(tim))

####################Slope Est##########################
print "Slope Est"
time=np.array(tim)
posX=np.array(px)
posY=np.array(py)
posZ=np.array(pz)
energy=np.array(eList)
logRat=np.array(ratList)
slope, intercept, r_value, p_value, std_err = stats.linregress(time,posY)
dY=sqrt(12*posY.var())
dt=sqrt(12*time.var())
print "b-a "+str(dt)
print "mean "+str(time.mean())
print "std "+str(time.std())
print "dY "+str(dY)
print "mean "+str(posY.mean())
print "std "+str(posY.std())
est=-dY/dt
print "slope est "+str(est)
slopeList.append(-slope)
bias=100*(slope-est)/slope
print "bias "+str(bias)

minT=time.min()
maxT=time.max()
ran=10
minT=-ran
maxT=ran
tRange=range(int(minT),int(maxT))
posRange=[]
posRangeEst=[]
t0=time.mean()
for i in tRange:
    posRange.append(slope*float(i)+intercept)
    posRangeEst.append(est*(float(i)-t0))
####################Time v Pos##########################
print 
print "time v pos"
print "slope "+str(slope)
print "intercept "+str(intercept)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
print "dt mean " +str(time.mean())
print "dt std " +str(time.std())
print
# plt.plot(time,py,'.')
# plt.plot(tRange,posRange,'o',color='r',)
# plt.plot(tRange,posRangeEst,'*',color='c')
# plt.ylabel("postion along the cell(mm)")
# plt.xlabel("delta time(ns)")
# plt.title("Postion along the Cell vs Time")
# plt.show()

####################Pos V Rat###########################
print "pos v rat"
slope, intercept, r_value, p_value, std_err = stats.linregress(logRat,posY)
print "slope "+str(slope)
print "intercept "+str(intercept)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
# plt.plot(logRat,posY,'.')
# plt.ylabel("postion along the cell(mm)")
# plt.xlabel("log ratio of signals")
# plt.title("log10(S0/S1) vs Pos")
# plt.show()
print
####################Time v Rat##########################
print "time v rat"
slope, intercept, r_value, p_value, std_err = stats.linregress(time,logRat)
print "slope "+str(slope)
print "intercept "+str(intercept)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
n=len(ratList)
stdSum=0
for i in xrange(0,len(ratList)):
    ratEst=slope*time[i]+intercept
    stdSum+=(ratList[i]-ratEst)**2
stdErr=sqrt(stdSum/len(ratList))

ran=10
minT=-ran
maxT=ran
tRange=range(int(minT),int(maxT))
logRange=[]
for i in tRange:
    logRange.append(slope*float(i)+intercept)

# plt.figure()
# plt.plot(tRange,logRange,'*',color='c')
# plt.plot(time,logRat,'.')
# plt.errorbar(tRange, logRange, yerr=1*stdErr, fmt='o',color='r',ecolor='r',capthick=2)
# plt.ylabel("log ratio of signals")
# plt.xlabel("delta time(ns)")
# plt.title("log10(S0/S1) vs Time")
# plt.show()
print
#############time hist***********
print "time "
print "time mean: "+str(time.mean())
print "time std: "+str(time.std())
# n, bins, patches=plt.hist(time,bins=31,range=[-15,15])
# plt.title("Delta Time (ns)")
# plt.show()
print
############# posY hist***********
print "postion"
print "posY mean: "+str(posY.mean())
print "posY std: "+str(posY.std())
print
# n, bins, patches=plt.hist(posY,bins=30)
# plt.title("Muon Position (cm)")
# plt.show()
#############logRat hist***********
print "log ratio"
print "logRat mean: "+str(logRat.mean())
print "logRat std: "+str(logRat.std())
print
# n, bins, patches=plt.hist(logRat,bins=30)
# plt.title("log ratio of signals")
# plt.show()
#############energy hist***********
print "energy "
print "energy mean: "+str(energy.mean())
print "energy std: "+str(energy.std())
print
# plt.show()
mu , sigma=norm.fit(eList)
n, bins, patches=plt.hist(energy,normed=True,bins=100,range=[0,500])
gaus = mlab.normpdf( bins, mu, sigma)
plt.title("Energy (Mev)")
plt.plot(bins, gaus, 'r--', linewidth=2)
plt.show()
print "mu "+str(mu)
print "mean "+str(energy.mean())
print "sigma "+ str(sigma)
print







print "###################################"
print 
print
print "###################################"
print "Cut calculation"
#################cuts#########################

delList=[]
slope, intercept, r_value, p_value, std_err = stats.linregress(time,logRat)

print "std err"+str(stdErr)
for i in xrange(0,len(ratList)):
    rat=ratList[i]
    t=time[i]
    e=eList[i]
    ratEst=slope*t+intercept
    eCut=e>75
    ratCut=abs(rat-ratEst)>1*stdErr
    cutList=[ratCut,eCut]
    for cut in cutList:
        if cut  :
            if i not in delList:
                delList.append(i)


print "deleted events "+str(len(delList))
for i in reversed(delList):
    del leftList[i]
    del rightList[i]
    del ratList[i]
    del totalList[i]
    del px[i]
    del py[i]
    del pz[i]
    del px0[i]
    del py0[i]
    del pz0[i]
    del px1[i]
    del py1[i]
    del pz1[i]
    del dx[i]
    del dy[i]
    del dz[i]
    del tim[i]
    del volList[i]
    del eList[i]

posBin=[]
for i in xrange(0,100):
    posBin.append([])

for i in xrange(0,len(totalList)):
    y=py[i]+500
    binNum=y/10
    posBin[int(binNum)].append(totalList[i])
sTotal=[]
yIndex=[]
min=time.min()
max=time.max()
interval=max-min
for i in xrange(0,100):
   pb= posBin[i]
   sTotal.append(np.array(pb).mean())
   yIndex.append(min+i*interval/100)

print "number of muon events cut "+str(len(tim))



####################Slope Est##########################
print "Slope Est"
time=np.array(tim)
posX=np.array(px)
posY=np.array(py)
posZ=np.array(pz)
energy=np.array(eList)
logRat=np.array(ratList)
slope, intercept, r_value, p_value, std_err = stats.linregress(time,posY)
dY=sqrt(12*posY.var())
dt=sqrt(12*time.var())
print "b-a "+str(dt)
print "std "+str(time.std())
print "dY "+str(dY)
print "std "+str(posY.std())
est=-dY/dt
print "slope est "+str(est)
slopeList.append(-slope)
bias=100*(slope-est)/slope
print "bias "+str(bias)
minT=time.min()
maxT=time.max()
ran=10
minT=-ran
maxT=ran
tRange=range(int(minT),int(maxT))
posRange=[]
posRangeEst=[]
t0=time.mean()
for i in tRange:
    posRange.append(slope*float(i)+intercept)
    posRangeEst.append(est*(float(i)-t0))

####################Time v Pos##########################
print 
print "time v pos"
print "slope "+str(slope)
print "intercept "+str(intercept)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
print "dt mean " +str(time.mean())
print "dt std " +str(time.std())
print
# plt.plot(time,py,'.')
# plt.plot(tRange,posRange,'o',color='r',)
# plt.plot(tRange,posRangeEst,'*',color='c')
# plt.ylabel("postion along the cell(mm)")
# plt.xlabel("delta time(ns)")
# plt.title("Postion along the Cell vs Time")
# plt.show()

####################Pos V Rat###########################
print "pos v rat"
slope, intercept, r_value, p_value, std_err = stats.linregress(logRat,posY)
print "slope "+str(slope)
print "intercept "+str(intercept)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
# plt.plot(logRat,posY,'.')
# plt.ylabel("postion along the cell(mm)")
# plt.xlabel("log ratio of signals")
# plt.title("log10(S0/S1) vs Pos")
# plt.show()
print
####################Time v Rat##########################
print "time v rat"
slope, intercept, r_value, p_value, std_err = stats.linregress(time,logRat)
print "slope "+str(slope)
print "intercept "+str(intercept)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
# plt.plot(time,logRat,'.')
# plt.ylabel("log ratio of signals")
# plt.xlabel("delta time(ns)")
# plt.title("log10(S0/S1) vs Time")
# plt.show()
print
#############time hist***********
print "time "
print "time mean: "+str(time.mean())
print "time std: "+str(time.std())
# n, bins, patches=plt.hist(time,bins=31,range=[-15,15])
# plt.title("Delta Time (ns)")
# plt.show()
print
############# posY hist***********
print "postion"
print "posY mean: "+str(posY.mean())
print "posY std: "+str(posY.std())
print
# n, bins, patches=plt.hist(posY,bins=30)
# plt.title("Muon Position (cm)")
# plt.show()
#############logRat hist***********
print "log ratio"
print "logRat mean: "+str(logRat.mean())
print "logRat std: "+str(logRat.std())
print
# n, bins, patches=plt.hist(logRat,bins=30)
# plt.title("log ratio of signals")
# plt.show()
#############energy hist***********
print "energy "
print "energy mean: "+str(energy.mean())
print "energy std: "+str(energy.std())
print
n, bins, patches=plt.hist(energy,normed=True,bins=100,range=[0,500])
plt.title("Energy (Mev)")
mu , sigma=norm.fit(eList)
gaus = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, gaus, 'r--', linewidth=2)
print "mu "+str(mu)
print "mean "+str(energy.mean())
print "sigma "+ str(sigma)
plt.show()
print
