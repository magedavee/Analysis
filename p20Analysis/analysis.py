from initROOT import initROOT
import ROOT
from ROOT import * 
import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import stats
from scipy.stats import norm 
from array import array

def sigma(ratList,time,slope):
    s=0
    for i in xrange(0,len(ratList)):
	rat=ratList[i]
	t=time[i]
	s+=(rat-slope*t)**2
    sig=sqrt(s/float(len(ratList)))
    return sig

initROOT()

fil=TFile("/home/mage/Data/p20Data/proData.root")
tree=ROOT.TTree()
tree=fil.Get("proData")
dt=array("f",[0])
tot=array("f",[0])
left=array("i",[0])
right=array("i",[0])
tree.SetBranchAddress("dt",dt)
tree.SetBranchAddress("tot",tot)
tree.SetBranchAddress("right",right)
tree.SetBranchAddress("left",left)
dtList=[]
totList=[]
rightList=[]
leftList=[]
ratList=[]
ent= tree.GetEntries()
for i in xrange(1,ent):
    tree.GetEntry(i)
    total=(left[0]+right[0])/20
    if left[0]!=0 and right[0]!=0 and total>10000 and total<17000:
	dtList.append(dt[0])
	rightList.append(right[0])
	leftList.append(left[0])
	totList.append(total)
	rat=float(left[0])/float(right[0])
	lr=log10(rat)
	ratList.append(lr)
    else :
	print total
slope, intercept, r_value, p_value, std_err = stats.linregress(dtList,ratList)
print "slope "+str(slope)
print "r "+str(r_value)
sig=sigma(ratList,dtList,slope)
print "sigma" +str(sig)
# plt.plot(dtList,ratList,'.')
# plt.ylabel("log ratio of signals")
# plt.xlabel("delta time(ns)")
# plt.title("log10(S0/S1) vs Time")
# plt.show()
print len(dtList)
print "mean "+str(np.array(dtList).mean())
print "std "+ str( np.array(dtList).std())
mu , sigma=norm.fit(totList)
print "mu "+str(mu)
print "sigma "+ str(sigma)

n, bins, patches =plt.hist(np.array(totList),bins=50,normed=True)
# n, bins, patches =plt.hist(np.array(ratList),bins=50,normed=True)
# n, bins, patches =plt.hist(np.array(dtList),bins=50,normed=True)
# gaus = mlab.normpdf( bins, mu, sigma)
# plt.plot(bins, gaus, 'r--', linewidth=2)
# plt.xlabel("PE")
plt.show()
