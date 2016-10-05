from initROOT import initROOT
import ROOT
from ROOT import gROOT, TCanvas, TF1,TFile,TTree
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm 
from scipy import stats
from array import array
from math import *
import cPickle as pickle
from funcAna import *


###############load in dataCall #########
mean = []
std =[]
# infile = open('dataCal', 'r')
# while 1:
    # try:
        # mean.append(pickle.load(infile))
        # std.append(pickle.load(infile))
    # except EOFError :
        # break
# infile.close()
corr=readPickle('GainMatch')
data=readPickle('TimeOff')
for i,d in enumerate(data):
    if i % 2 == 0 :
        mean.append(d)
    else :
        std.append(d)
###############Const##################
print "mean"
print mean
c=29.9792
lat=14.4198
ptime=lat/c
M=4.094877445468857
mList=[]
layer=0.0#mean[2]-ptime
initROOT()
dir="../Data/p50/"
name=dir+"s004.root"
print "processing "+str(name)

###############Varibles##################
dt=[]
logRat=[]
for i in xrange(6):
    dt.append([])
    logRat.append([])
ddt=[]
lenCal=[]
pulse=[]
pmt=ROOT.PmtData(name,"p50_s4_setup.root")
ent=pmt.GetEntries()
tanl=[]
cosl=[]
sinl=[]
# tant=[]
theta=[[],[],[]]
tCheck=[]
dAnti=[]
print "what it should be "+str(ptime)
for i in xrange(0,ent):
    pmt.SetEntry(i)
    printPercent(i,ent)
    q,t,lr=getPara(pmt,i)
    fixGain(q,corr)
    tot,totT,totB=getTots(q)
    if cut(t,q) :
        fixTime(t,mean,ptime=ptime)
        # fixTime(t,mean)
        tt=t[0]-t[1]
        ddt.append(tt)
        l=sqrt((M*tt)**2+(lat)**2)
        t0=50.0/M-t[0]
        t2=50.0/M-t[1]
        tcheck=t[2]-l/c
        tCheck.append(tcheck)
        lenCal.append(l)
        ta=M*tt/lat
        s=M*tt/l
        c=lat/l
        tanl.append(ta)
        sinl.append(s)
        cosl.append(c)
        theta[0].append((180.0/pi)*np.arctan(ta))
        theta[1].append((180.0/pi)*np.arcsin(s))
        theta[2].append((180.0/pi)*np.arccos(c))
        for j,d in enumerate(dt):
            d.append(t[j])
            logRat[j].append(lr[j])
        for j in xrange(0,len(pulse)):
            pulse[j].append(q[j])
###############Analysis Setup #########
name=['cell_1','cell_2','left','0_3cross','right ','1_2cross']
printHist(dt)
tanl=np.array(tanl)
cosl=np.array(cosl)
sinl=np.array(sinl)
for T in theta :
    T=np.array(T)
tCheck=np.array(tCheck)
# tant=np.array(tant)
tCheck=np.array(tCheck)
###############Analysis #########
print '##############Hist#############'
plotHist(tanl,name="Hist tanl")
plotHist(cosl,name="Hist cosl")
plotHist(sinl,name="Hist sinl")
for i,T in enumerate(theta) :
    plotHist(T,name="Hist theta"+str(i))
plotHist(tCheck,name="Hist tCheck")
plotHist(lenCal,name="Hist cal lenght")
plotHist(ddt,name='Hist ddt')
