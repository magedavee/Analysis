from initROOT import initROOT
import ROOT
from ROOT import gROOT, TCanvas, TF1,TFile,TTree
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os.path
from scipy.stats import norm 
from array import array
import sys
import math
from timeit import default_timer
start= default_timer()
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
# pmt.CalIntegral(0)
# pmt.CalIntegral(1)
ent=pmt.GetEntries()
# ent =10000
nChList=[]
dtList=[]
integral=[]
totalCharge=[]
dt=[]
logSig=[]
dtBin=[]
logSigBin=[]
totAve=[]
logAve=[]
timeBin=[]
for i in xrange (0,2):
    dt.append([])
    totalCharge.append([])
    integral.append([])
    logSig.append([])
    dtList.append([])
    dtBin.append([])
    logSigBin.append([])
    totAve.append([])
    logAve.append([])
    timeBin.append([])

for i in xrange(0,ent):
    if i%10**4==0:
        print str(float(int(i*1000/ent))/10.0)+" %"
    pmt.SetEntry(i)
    
    charge=[0,0,0,0]
    for j in xrange(0,len(charge)):
        charge[j]=pmt.GetPulseIntegral(j,i)
    
    dt[0].append(pmt.DeltaT(0,1))
    dt[1].append(pmt.DeltaT(2,3))
    totalCharge[0].append(charge[0]+charge[1])
    totalCharge[1].append(charge[2]+charge[3])
    if totalCharge[0][i]>=400000 and abs(dt[0][i])<=20:
        if totalCharge[1][i]>=400000 and abs(dt[1][i])<=20:
            for cell in xrange(0,2):
                integral[cell].append(totalCharge[cell][i])
                dtList[cell].append(dt[cell][i])
                logSig[cell].append(math.log10(totalCharge[cell][i]))
    # print str(i)+" event has charge " +str(charge)
    nch=pmt.GetNCha()
    nChList.append(nch)
    # for i in xrange(0,2):
            # integral[1][1].append(charge[1])
            # tot.append(totalCharge)
            # dtList.append(dt[i])
            # tot.append

        # if totalCharge>=400000 or abs(dt[0])<=20:
            # dtList[len(dtList)].append(dt[0])


pmt.Write("p50_s4_setup.root")
for cell in xrange(0,2):
    print "events : "+str(len(dtList[cell]))
    max=np.array(dtList[cell]).max()
    min=np.array(dtList[cell]).min()
    spread=max-min
    for i in xrange(0, spread):
        dtBin[cell].append([])
        logSigBin[cell].append([])
    for i in xrange(0, len(integral[cell])):
        index=int(dtList[cell][i])
        dtBin[cell][index].append(integral[cell][i])
    for i in xrange(0, spread):
        timeBin[cell].append(i+min)
        if len(dtBin[cell][i]) !=0 :
            totAve[cell].append(np.array(dtBin[cell][i]).mean())
        else:
            totAve[cell].append(0)
        if len(logSigBin[cell][i]) !=0 :
            logAve[cell].append(np.array(logSigBin[cell][i]).mean())
        else:
            logAve[cell].append(0)
    sd=np.array(dtList[cell]).std()
    mean=np.array(dtList[cell]).mean()
    b_a=math.sqrt(12)*sd
    plt.hist(dtList[cell],bins=20)
    plt.title("Delta Time for cell "+str(cell+1))
    plt.xlabel("delta time(ns)")
    #plt.show()
    print "cell "+str(cell)
    print "b-a time "+str(b_a)
    print "mean time: "+str(mean)
    print "std time: "+str(mean)
    sd=np.array(logSig[cell]).std()
    mean=np.array(logSig[cell]).mean()
    plt.hist(logSig[cell],bins=30)
    plt.title("log signal for cell "+str(cell))
    plt.xlabel("log signal")
    #plt.show()
    print "mean log sig: "+str(mean)
    print "std log sig: "+str(mean)
    plt.hist(integral[cell],bins=100)
    plt.title("Pulse Integral of cell "+str(cell))
    plt.xlabel("pulse integral")
    #plt.show()
    mean=np.array(integral[cell]).mean()
    sd=np.array(integral[cell]).std()
    print "pulse mean: "+str(mean)
    print "pulse sd: "+str(sd)
    plt.plot(timeBin[cell],totAve[cell],'.')
    plt.xlabel("ave e vs dt")
    #plt.show()
    plt.plot(dtList[cell],logSig[cell],'.')
    plt.xlabel("dt vs log sig")
    #plt.show()
    plt.plot(timeBin[cell],logAve[cell],'.')
    plt.xlabel("dt vs log sig")
    #plt.show()

timeOfLoop=default_timer() - start
print "the program ran in "+str(timeOfLoop)+ " sec"
# plt.hist(integral[0],bins=100)
# plt.title("Pulse Integral of pmt 1")
# plt.xlabel("pulse integral")
# #plt.show()
# print "input"
# raw_input()
# plt.hist(integral[1],bins=100)
# plt.title("Pulse Integral of pmt 2")
# plt.xlabel("pulse integral")
# #plt.show()

# plt.hist(tot,bins=100)
# plt.title("Pulse Integral pmt 1+ pmt2")
# plt.xlabel("combined pulse integral")
# #plt.show()
