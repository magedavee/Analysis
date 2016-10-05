from initROOT import initROOT
import ROOT
from ROOT import gROOT, TCanvas, TF1,TFile,TTree
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os.path
from scipy.stats import * 
from array import array
import sys
from math import *
import cPickle as pickle
from funcAna import *

##################Init###########
initROOT()
dir="../Data/p50/"
name=dir+"s004.root"
# name=dir+"p50.root"
if os.path.isfile(name):
    print "looking at "+name
else:
    print "error "+name+" does not exist"
    quit()
setup="p50_s4_setup.root"
if os.path.isfile("./"+setup):
    pmt=ROOT.PmtData(name,setup)
else:
    print "no setup "+setup
    pmt=ROOT.PmtData(name)
###############Setup varibles#########
ent=pmt.GetEntries()
dt=[]
ddt=[]
pulse=[]
c=29.9792
lat=14.4198
M=4.094877445468857

for i in xrange (0,6):
    dt.append([])
    if i<4:
        pulse.append([])

corr=readPickle('GainMatch')
print "corr"
print corr

###############Load Data #########

for i in xrange(0,ent):
    pmt.SetEntry(i)
    printPercent(i,ent)
    q,t,lr=getPara(pmt,i)
    for j in xrange(0,len(q)):
        q[j]=pmt.GetPulseIntegral(j,i)
    fixGain(q,corr)
    tot,totT,totB=getTots(q)

    if cut(t,q) :
        tt=t[0]-t[1]
        # if abs(int(tt))< numOfBin :
            # cellDtList[int(abs(tt))].append(tot)
        dt[0].append(t[0])###top cell 
        dt[1].append(t[1])###bottom cell
        dt[2].append(t[2])###left pmt
        dt[3].append(t[3])###top left bot right
        dt[4].append(t[4])###right 
        dt[5].append(t[5])###top right bot left
        ddt.append(t[0]-t[1])
ddt=np.array(ddt)

for i in xrange (0,len(dt)):
    dt[i]=np.array(dt[i])

###############Analysis #########
totMeanList=[]
# for i in cellDtList :
    # totMeanList.append(np.array(i).mean())


name=['cell_1','cell_2','left','0_3cross','right ','1_2cross']
printHist(dt)
plotHist(ddt,"ddt")

print len(ddt)
out=[]
for i in xrange(0,len(dt)):
    out.append(dt[i].mean())
    out.append(dt[i].std())
writePickle('TimeOff',out)
