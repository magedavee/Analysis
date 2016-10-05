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
initROOT()
name=[]
name.append("/home/mage/Data/p20Data/root/P20_2015-04-16-10-44-25.root")
# name.append("/home/mage/Data/p20Data/root/P20_2015-04-16-09-44-16.root")
# name.append("/home/mage/Data/p20Data/root/P20_2015-04-16-08-44-07.root")
dtList=[]
integral0=[]
integral1=[]
total=[]
reduc=[]
pulseListL=[]
pulseListR=[]
logRat=[]
for j in xrange(0,1):
    print "processing "+str(name[j])
    pmt=ROOT.PmtData(name[j],"p20Setup.root")
    # event=0
    ent= pmt.GetEntries()
    # # ent=1000
    for i in xrange(0,ent):
        # print i
        pmt.SetEntry(i)
        i0= pmt.GetPulseIntegral(0,i)
        i1= pmt.GetPulseIntegral(1,i)
        dt=pmt.DeltaT(0,1)
        tot=i0+i1
        upB=160000
        lowB=100000
        iUp=i0<upB and i1 <upB
        iDn=i0>lowB  and i1 > lowB
        iCheck= iUp and iDn and tot<310000
        if pmt.GetNCha()==4 and tot >14000 and iCheck :
            total.append(tot)
            pulseListL.append(i0)
            pulseListR.append(i1)
            dtList.append(dt)
            rat=float(i0)/float(i1)
            lr=log10(rat)
            logRat.append(lr)
            if tot>240000 and tot < 300000 :
                reduc.append(tot)

            # pmt.GetTrace(0).Draw()
            # pmt.GetTrace(1).Draw("SAME")
            # print tot
            # print i0
            # print i1
            # raw_input("pause ")
############dT###################
print "number "+str(len(dtList))
mean=np.array(dtList).mean()
std=np.array(dtList).std()
print "dt mean "+ str(mean)
print "dt spread "+ str(sqrt(12)*std)
print "dt std "+ str(std)
print

# n, bins, patches =plt.hist(np.array(dtList),bins=20)
# plt.xlabel("delta time (ns)")
# plt.show()
############logRat###################
mean=np.array(logRat).mean()
std=np.array(logRat).std()
print "log mean "+ str(mean)
print "log std "+ str(std)
print

# n, bins, patches =plt.hist(np.array(logRat),bins=20)
#plt.xlabel("Log Ratio (S0/S1)")
#plt.show()
############logRat vs Dt###################
n=len(logRat)
stdSum=0
slope, intercept, r_value, p_value, std_err = stats.linregress(dtList,logRat)
for i in xrange(0,len(logRat)):
    ratEst=slope*dtList[i]+intercept
    stdSum+=(logRat[i]-ratEst)**2
stdErr=sqrt(stdSum/len(logRat))

ran=10
minT=-ran
maxT=ran
tRange=range(int(minT),int(maxT))
logRange=[]
for i in tRange:
    logRange.append(slope*float(i)+intercept)

print "slope "+str(slope)
print "intercept "+str(intercept)
print "r_value "+str(r_value)
print "stdErr "+str(stdErr)
#plt.plot(dtList,logRat,'.')
#plt.errorbar(tRange, logRange, yerr=1*stdErr, fmt='o',color='g',ecolor='g',capthick=2)
#plt.ylabel("log ratio of signals")
#plt.xlabel("delta time(ns)")
#plt.title("log10(S0/S1) vs Time")
#plt.show()
print 
############Pulse###################
mu , sigma=norm.fit(reduc)
n, bins, patches =plt.hist(np.array(total),bins=100,normed=True)
gaus = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, gaus, 'r--', linewidth=2)
plt.xlabel("Pulse Integral")
plt.show()

print "pulse mean "+ str(np.array(total).mean())
print "pulse std "+ str(np.array(total).std())
print "mu "+str(mu)
print "sigma "+str(sigma)
print 

############cut###################
print "cutting"
delList=[]
for i in xrange(0,len(logRat)):
    rat=logRat[i]
    t=dtList[i]
    ratEst=slope*t+intercept
    ratCut=abs(rat-ratEst)>1*stdErr
    cutList=[ratCut]
    for cut in cutList:
        if cut  :
            if i not in delList:
                delList.append(i)

print "deleted events "+str(len(delList))
for i in reversed(delList):
    del logRat[i]
    del dtList[i]
    del total[i]
    del pulseListL[i]
    del pulseListR[i]
reduc=[]
for tot in total:
    if tot>240000 and tot < 310000 :
        reduc.append(tot)
print ##############################
print ##############################
############dT###################
print "number "+str(len(dtList))
mean=np.array(dtList).mean()
std=np.array(dtList).std()
print "dt mean "+ str(mean)
print "dt spread "+ str(sqrt(12)*std)
print "dt std "+ str(std)
print

# n, bins, patches =plt.hist(np.array(dtList),bins=20)
#plt.xlabel("delta time (ns)")
#plt.show()
############logRat###################
mean=np.array(logRat).mean()
std=np.array(logRat).std()
print "log mean "+ str(mean)
print "log std "+ str(std)
print

# n, bins, patches =plt.hist(np.array(logRat),bins=20)
#plt.xlabel("Log Ratio (S0/S1)")
#plt.show()
############logRat vs Dt###################
n=len(logRat)
stdSum=0
slope, intercept, r_value, p_value, std_err = stats.linregress(dtList,logRat)

ran=10
minT=-ran
maxT=ran
tRange=range(int(minT),int(maxT))
logRange=[]
for i in tRange:
    logRange.append(slope*float(i)+intercept)

print "slope "+str(slope)
print "intercept "+str(intercept)
print "r_value "+str(r_value)
#plt.plot(dtList,logRat,'.')
#plt.errorbar(tRange, logRange, yerr=1*stdErr, fmt='--',color='r',ecolor='g',capthick=2)
#plt.ylabel("log ratio of signals")
#plt.xlabel("delta time(ns)")
#plt.title("log10(S0/S1) vs Time")
#plt.show()
print 
############Pulse###################
mu , sigma=norm.fit(reduc)
gaus = mlab.normpdf( bins, mu, sigma)
n, bins, patches =plt.hist(np.array(total),bins=100,normed=True)
gaus = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, gaus, 'r--', linewidth=2)
plt.xlabel("Pulse Integral")
plt.show()

print "pulse mean "+ str(np.array(total).mean())
print "pulse std "+ str(np.array(total).std())
print "mu "+str(mu)
print "sigma "+str(sigma)
