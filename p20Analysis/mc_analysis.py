import numpy as np
import math
import matplotlib.pyplot as plt
import ROOT
from ROOT import gROOT
from math import *
from array import array
from scipy import stats
from scipy.stats import norm
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import Axes3D as ax3

def sigma(ratList,time,slope):
    s=0
    for i in xrange(0,len(ratList)):
	rat=ratList[i]
	t=time[i]
	s+=(rat-slope*t)**2
    sig=sqrt(s/float(len(ratList)))
    return sig




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
bList=[]
mList=[]
rList=[]
slopeList=[]
estList=[]
zCutList=[]
# for zCut in xrange(1,1000,10):
# print zCut
zCut=8000
zCutList.append(float(zCut)/100.0)
ratList=[]
for i in xrange(0,entries):
    tree.GetEntry(i)
    # if abs(t[0])!=0 and abs(pid[0])==13:
    if abs(t[0])<30  and abs(t[0])!=0 and abs(x[0])<(80) and abs(z[0])<(80) and abs(y[0])<(500.0) and  not (math.isnan(x0[0]) or math.isnan(x1[0])) and abs(pid[0])==13 :
	    if (z0[0]-z1[0])>0 and right[0]>0 and left[0]>0 :
		# if y[0]>(-57*t[0]-100) and y[0]<(-57*t[0]+100): 
		# if y[0]<(-57*t[0]-100) or y[0]>(-57*t[0]+100):
		# if e[0]<50:
		    leftList.append(left[0])
		    rightList.append(right[0])
		    rat=float(left[0])/float(right[0])
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
print "t v p"
mList.append(len(tim))
time=np.array(tim)
posX=np.array(px)
posY=np.array(py)
posZ=np.array(pz)
energy=np.array(eList)
dY=sqrt(12*posY.var())
dt=sqrt(12*time.var())
print "b-a "+str(dt)
print "std "+str(time.std())
print "dY "+str(dY)
print "std "+str(posY.std())
est=-dY/dt
estList.append(-est)
print "slope est "+str(est)
slope, intercept, r_value, p_value, std_err = stats.linregress(time,posY)
slopeList.append(-slope)
print "slope "+str(slope)
bias=100*(slope-est)/slope
bList.append(bias)
print "bias "+str(bias)
print 
print "intercept "+str(intercept)
rList.append(r_value)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
print "dt mean " +str(time.mean())
print sigma(time,posY,slope)
print
print "p v r"
slope, intercept, r_value, p_value, std_err = stats.linregress(ratList,posY)
print "slope "+str(slope)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
print sigma(ratList,posY,slope)
print
print "t v r"
slope, intercept, r_value, p_value, std_err = stats.linregress(time,ratList)
print "slope "+str(slope)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)

sig= sigma(ratList,time,slope)
print sig
plt.plot(time,py,'.')
plt.ylabel("postion along the cell(mm)")
plt.xlabel("delta time(ns)")
plt.title("Postion along the Cell vs Time")
plt.show()
plt.plot(time,ratList,'.')
plt.ylabel("log ratio of signals")
plt.xlabel("delta time(ns)")
plt.title("log10(S0/S1) vs Time")
plt.show()
plt.plot(ratList,py,'.')
plt.ylabel("postion along the cell(mm)")
plt.xlabel("log ratio of signals")
plt.title("log10(S0/S1) vs Pos")
plt.show()
px1=[]
py1=[]
pz1=[]
px01=[]
py01=[]
pz01=[]
px11=[]
py11=[]
pz11=[]
dx1=[]
dy1=[]
dz1=[]
leftList1=[]
rightList1=[]
tim1=[]
volList1=[]
cellLen=950
eList1=[]
bList1=[]
mList1=[]
rList1=[]
slopeList1=[]
estList1=[]
zCutList1=[]
ratList1=[]
for i in xrange(0,len(ratList)):
    if i >= len(ratList):
	break
    rat=ratList[i]
    t=time[i]
    low=slope*t-sig
    hi=slope*t+sig
    if abs(rat)<1 and rat>low and  rat<hi:
	leftList1.append(leftList[i])
	rightList1.append(rightList[i])
	ratList1.append(ratList[i])
	px1.append(px[i])
	py1.append(py[i])
	pz1.append(pz[i])
	tim1.append(tim[i])
	eList1.append(eList[i])

print "number of muon events cut2 "+str(len(tim1))
print "t v p"
mList.append(len(tim))
time1=np.array(tim1)
posX1=np.array(px1)
posY1=np.array(py1)
posZ1=np.array(pz1)
energy1=np.array(eList1)
dY=sqrt(12*posY1.var())
dt=sqrt(12*time1.var())
print "b-a "+str(dt)
print "std "+str(time1.std())
print "dY "+str(dY)
print "std "+str(posY1.std())
est=-dY/dt
estList.append(-est)
print "slope est "+str(est)
slope, intercept, r_value, p_value, std_err = stats.linregress(time1,posY1)
slopeList.append(-slope)
print "slope "+str(slope)
bias=100*(slope-est)/slope
bList.append(bias)
print "bias "+str(bias)
print 
print "intercept "+str(intercept)
rList.append(r_value)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
print "dt mean " +str(time.mean())
print
print "p v r"
slope, intercept, r_value, p_value, std_err = stats.linregress(ratList1,posY1)
print "slope "+str(slope)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
print
print "t v r"
slope, intercept, r_value, p_value, std_err = stats.linregress(time1,ratList1)
print "slope "+str(slope)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)

############scatter###########
# plt.plot(time,posZ,'.')
# plt.show()
# plt.plot(time,posX,'.')
# plt.show()
plt.plot(time1,py1,'.')
plt.ylabel("postion along the cell(mm)")
plt.xlabel("delta time(ns)")
plt.title("Postion along the Cell vs Time")
plt.show()
plt.plot(time1,ratList1,'.')
plt.ylabel("log ratio of signals")
plt.xlabel("delta time(ns)")
plt.title("log10(S0/S1) vs Time")
plt.show()
plt.plot(ratList1,py1,'.')
plt.ylabel("postion along the cell(mm)")
plt.xlabel("log ratio of signals")
plt.title("log10(S0/S1) vs Pos")
plt.show()
# plt.plot(pz,eList,'.')
# plt.show()
# plt.plot(dz,time,'.')
# plt.show()
# plt.plot(eList,time,'.')
# plt.ylabel("postion along the cell(mm)")
# plt.xlabel("delta time(ns)")
# plt.title("Postion along the Cell vs Time")
# plt.show()
#############time hist***********
n, bins, patches=plt.hist(time1,bins=31,range=[-15,15])
# n, bins, patches=plt.hist(dx,bins=31,range=[-80,1200])
# print n
# print bins
plt.title("Delta Time (ns)")
plt.show()
n, bins, patches=plt.hist(np.array(py1),bins=30)
plt.title("Muon Position (cm)")
# plt.title("Energy (Mev)")
# mu , sigma=norm.fit(eList)
# gaus = mlab.normpdf( bins, mu, sigma)
# plt.plot(bins, gaus, 'r--', linewidth=2)
# print "mu "+str(mu)
# print "mean "+str(energy.mean())
# print "sigma "+ str(sigma)
plt.show()
n, bins, patches=plt.hist(np.array(ratList1),bins=30)
plt.title("log ratio of signals")
plt.show()
#############zcut***********
# plt.plot(zCutList,slopeList,'.')
# plt.ylabel("slope")
# plt.xlabel("zCut")
# # plt.title("Postion along the Cell vs Time")
# plt.show()
# plt.plot(zCutList,estList,'.')
# plt.ylabel("est")
# plt.xlabel("zCut")
# # plt.title("Postion along the Cell vs Time")
# plt.show()
# plt.plot(zCutList,bList,'.')
# plt.ylabel("bias")
# plt.xlabel("zCut")
# # plt.title("Postion along the Cell vs Time")
# plt.show()
# plt.plot(zCutList,mList,'.')
# plt.ylabel("count")
# plt.xlabel("zCut")
# # plt.title("Postion along the Cell vs Time")
# plt.show()
# plt.plot(zCutList,rList,'.')
# plt.ylabel("r:")
# plt.xlabel("zCut")
# # plt.title("Postion along the Cell vs Time")
# plt.show()
