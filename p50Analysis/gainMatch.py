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
c=29.9792
lat=14.4198
M=4.094877445468857
dt=[]
logRat=[]
for i in xrange(6):
    dt.append([])
    logRat.append([])
pulse=[]
for i in xrange(4):
    pulse.append([])
corr=[]
corr.append(1)
totList=[]
qTopList=[]
qBotList=[]

###############Load Data #########

###################%%%%%%%%%%%%%%#########################
for i in xrange(0,ent):
    pmt.SetEntry(i)
    printPercent(i,ent)

    q,t,lr=getPara(pmt,i)
    tot,totT,totB=getTots(q)


    if cut(t,q) :
        totList.append(tot)
        qTopList.append(totT)
        qBotList.append(totB)
        for i,d in enumerate(dt):
            d.append(t[i])
            logRat[i].append(lr[i])
        for j in xrange(0,len(pulse)):
            pulse[j].append(q[j])

###############Analysis #########
name=['cell_1','cell_2','left','0_3cross','right ','1_2cross']
comp=[0,3,4]
for i in comp :
    print '###########'
    mean,std=plotHist(logRat[i],"logRat_"+name[i])
    co=10**mean
    print "The correction factor "+str(co)
    corr.append(co)

writePickle('GainMatch',corr)
print len(pulse[0])
print '###########'

