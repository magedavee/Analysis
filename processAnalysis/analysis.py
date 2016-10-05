import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import gROOT
from math import *
from array import array
from scipy import stats
def cut():
    print "cuts"


inFile=ROOT.TFile("proc_cry_total.root")
# inFile=ROOT.TFile("proc_cry_total9.root")
# inFile=ROOT.TFile("proc_cry_0.root.pro")
# inFile=ROOT.TFile("cry_0.root.pro")
tree=ROOT.TTree()
tree=inFile.Get("photon_Data")
# y=np.zeros(1,dtype=float)
x=array("f",[0])
y=array("f",[0])
z=array("f",[0])
t=array("f",[0])
tl=array("f",[0])
tr=array("f",[0])
e=array("f",[0])
mX=array("f",[0])
mY=array("f",[0])
mZ=array("f",[0])
pid=array("i",[0])
vol=array("i",[0])
l=array("i",[0])
r=array("i",[0])

tree.SetBranchAddress("x",x)
tree.SetBranchAddress("y",y)
tree.SetBranchAddress("z",z)
tree.SetBranchAddress("t",t)
tree.SetBranchAddress("tl",tl)
tree.SetBranchAddress("tr",tr)
tree.SetBranchAddress("e",e)
tree.SetBranchAddress("momX",mX)
tree.SetBranchAddress("momY",mY)
tree.SetBranchAddress("momZ",mZ)
tree.SetBranchAddress("pid",pid)
tree.SetBranchAddress("vol",vol)
tree.SetBranchAddress("left",l)
tree.SetBranchAddress("rightt",r)
entries=int(tree.GetEntries())
# entries=6847
px=[]
py=[]
pz=[]
eList=[]
tim=[]
timL=[]
timR=[]
momX=[]
momY=[]
momZ=[]
volList=[]
lef=[]
righ=[]
tot=[]
rat=[]
cellLen=950
logRat=[]
for i in xrange(0,entries):
    tree.GetEntry(i)
    print r[0]
    # if r[0]>0:
    # if abs(t[0])!=0 and abs(pid[0])==13 and tl[0]<20 and tr[0] < 20:
    if abs(t[0])!=0 and abs(x[0])<(80.0) and abs(z[0])<(80.0) and abs(y[0])<(500.0)and abs(pid[0])==13:
        print r[0]
        px.append(x[0])
        py.append(y[0])
        pz.append(z[0])
        tim.append(t[0])
        timL.append(tl[0])
        timR.append(tr[0])
        volList.append(vol[0])
        momX.append(mX[0])
        momY.append(mY[0])
        momZ.append(mZ[0])
        eList.append(e[0])
        lef.append(l[0])
        righ.append(r[0])
        lr=l[0]+r[0]
        # rat=float(l[0])/float(r[0])
        # lograt=log10(rat)
        # logRat.append(lograt)
        # if lr==0:
            # print "what"
        # ra=float(l[0])/float(lr)
        # rat.append(ra)
        tot.append(lr)
print "number of muon events "+str(len(tim))
time=np.array(tim)
energy=np.array(eList)
timeL=np.array(timL)
timeR=np.array(timR)
posX=np.array(px)
posY=np.array(py)
posZ=np.array(pz)
left=np.array(lef)
right=np.array(righ)
total=np.array(tot)
ratio=np.array(rat)
dy=sqrt(12*posY.var())
dt=sqrt(12*time.var())
print "b-a "+str(dt)
print "std "+str(time.std())
print "dy "+str(dy)
print "std "+str(posY.std())
est=-dy/dt
print "slope est "+str(est)
slope, intercept, r_value, p_value, std_err = stats.linregress(time,posY)
print "slope "+str(slope)
print "bias "+str(100*(slope-est)/slope)
print 
print "intercept "+str(intercept)
print "r "+str(r_value)
print "p "+str(p_value)
print "std "+str(std_err)
############scatter###########
# plt.plot(time,posZ,'.')

# plt.show()
# plt.plot(time,posX,'.')
# plt.show()
plt.plot(posY,time,'.')
plt.ylabel("postion along the cell(mm)")
plt.xlabel("delta time(ns)")
plt.show()
# plt.plot(momZ,posZ,'.')
# plt.ylabel("postion along the cell(mm)")
# plt.xlabel("delta time(ns)")
# plt.xlabel("mom")
# plt.ylabel("postion")
# plt.title("Postion along the Cell vs Time")
# plt.show()
#############time hist***********
n, bins, patches=plt.hist(energy,bins=100,range=[0,200])
# print np.array(n).mean()
# print n
# print bins
plt.title("Delta Time (ns)")
plt.show()
# n, bins, patches=plt.hist(posY,bins=31,range=[-600,600])
# plt.title("Muon Position (cm)")
# plt.show()
# n, bins, patches=plt.hist(total,bins=31,range=[0,70000])
# plt.show()
