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

c=29.9792458
dtbin=41
dthalfbin=(dtbin+1)/2

global pmtLab
pmtLab=[]
pmtLab.append("pmt 1 - pmt 0 t0+d1-t0=d1 ")
pmtLab.append("pmt 3 - pmt 2 to+d1+tp+d2+d3-to-d1-tp-d2=d3 ")
pmtLab.append("pmt 2 - pmt 0 tp+d1+d2 ")
pmtLab.append("pmt 2 - pmt 1 tp+d2")
pmtLab.append("pmt 3 - pmt 0 tp+d1+d2+d3")
pmtLab.append("pmt 3 - pmt 1 tp+d2+d3")

global y
global spreadY
y=[]
spreadY=[]

def dd(y1,y2):
    result=[]

    for i in xrange(0,len(y1)) :
            result.append(y1[i]-y2[i])
    return np.array(result)
def printArray (col, show=False):
    print pmtLab[col]
    print y[:,col]
    print y[dthalfbin,col]
    print spreadY[dthalfbin,col]
    if show:
        plt.title(pmtLab[col])
        plt.plot(dtAx,c*y[:,col])
        plt.show()
def printDic(key,show=False):
    print "key is "+key
    print res[key]
    print "crossing dt "+str(res[key][20])
    if show:
        plt.plot(dtAx,c*res[key])
        plt.show()
# dd([[0]],0,2)
initROOT()
dir="../Data/p50/"
name=dir+"s004.root"
dtList=[[],[],[],[],[],[]]
integral0=[]
integral1=[]
total=[]
reduc=[]
pulseListL=[]
pulseListR=[]
logRat=[]
print "processing "+str(name)
pmt=ROOT.PmtData(name,"p50_s4_setup.root")
# event=0
ent= pmt.GetEntries()
# # ent=1000
dtHist=[]
dtAx=[]
dt_pos=[]
for i in xrange(0,dtbin):
    dtHist.append([[],[],[],[],[],[]])
    dtAx.append(i-20)

for i in xrange(0,ent):
    # print i
    # if i%10**4==0:
        # print str(float(int(i*1000/ent))/10.0)+" %"
    pmt.SetEntry(i)
    pulse=[]
    for i in range(4):
        # print pmt.GetPulseIntegral(0,i)
        pulse.append(pmt.GetPulseIntegral(0,i))
    dt=[]
    dt.append(pmt.DeltaT(1,0))
    dt.append(pmt.DeltaT(3,2))
    dt.append(pmt.DeltaT(2,0))
    dt.append(pmt.DeltaT(2,1))
    dt.append(pmt.DeltaT(3,0))
    dt.append(pmt.DeltaT(3,1))
    tot=pulse[0]+pulse[1]
    upB=180000
    lowB=90000
    pCheck=[]
    for i in range(4):
        pCheck.append(pulse[i]<upB and pulse[i]> lowB)
    dtch1=abs(dt[0])<20 and abs(dt[1])<20
    dtch2=abs(dt[2])<40 and abs(dt[3])<40 and abs(dt[4])<40 and abs(dt[5])<40  
    ddt=(dt[0] - dt[1])+dthalfbin-1
    ddtch=abs(ddt)<dtbin
    strch=abs(ddt)<4
    chach = pmt.GetNCha()==4
    # if  ddtch and chach and tot >14000 :# and dtch1 and dtch1 and dtch2 :
    dt12=[dt[0],dt[1]]
    dt_pos.append(dt12)
    # else :
        # if tot>310000 :
            # pmt.GetTrace(0).Draw()
            # pmt.GetTrace(1).Draw("SAME")
            # print tot
            # raw_input()
    if chach and tot >14000 :
        # print dtch1
        # print dtch2
        # print ddtch
        if dtch1 and ddtch and dtch2 :
            total.append(tot)
            pulseListL.append(pulse[0])
            pulseListR.append(pulse[1])
            for i in xrange(0,6):
                    dtHist[ddt][i].append(dt[i])
                    dtList[i].append(dt[i])
    rat=float(pulse[0])/float(pulse[1])
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
dt_pos=np.array(dt_pos)
print
print dt_pos
mean=dt_pos[:,0].mean()
std=dt_pos[:,0].std()
print "number of events : "+str(len(dt_pos[:,0]))
print "mean 1: "+str(mean)
print 
print "std 1: "+str(std)
print "spread 1: "+str(sqrt(12)*std)
plt.hist(dt_pos[:,0])
plt.show()
print "#################################"
mean=dt_pos[:,1].mean()
std=dt_pos[:,1].std()
print "mean 2: "+str(mean)
print "std 2: "+str(std)
print "spread 2: "+str(sqrt(12)*std)
plt.hist(dt_pos[:,1])
plt.show()
print
for i in xrange(0,6):
    print pmtLab[i]
    print "number "+str(len(dtList[i]))
    mean=np.array(dtList[i]).mean()
    std=np.array(dtList[i]).std()
    print "dt mean "+ str(mean)
    print "dt spread "+ str(sqrt(12)*std)
    print "dt std "+ str(std)
    print

    # n, bins, patches =plt.hist(np.array(dtList[i]),bins=20)
    # plt.title("delta time (ns) "+str(i))
    # plt.xlabel("Mean: "+str(mean)+" std: "+str(std))
    # plt.show()
n, bins, patches =plt.hist(np.array(total),bins=20)
plt.show()
for i in xrange(0,len(dtHist)):
    y.append([])
    spreadY.append([])
    for j in xrange(0,len(dtHist[i])):
        y[i].append(0)
        spreadY[i].append(0)
        if len(dtHist[i][j])!=0:
            # print "dtHist number for hist :"+str(j)+" and ddt = "+str(i-20)
            # print len(dtHist[i][j])
            y[i][j]=np.array(dtHist[i][j]).mean()
            spreadY[i][j]=sqrt(12)*np.array(dtHist[i][j]).std()
for i in xrange(0,6):
    # print pmtLab[i]
    ct=[]
    for j in xrange(0,dtbin):
        ct.append(c*y[j][i])
    # plt.plot(dtAx,ct)
    # plt.title(pmtLab[i])
    # plt.show()
y=np.array(y)
spreadY=np.array(spreadY)
# resdthalfbin0=dd(y[:,2],y[:,5]) 
# res310=dd(y[:,4],y[:,3]) 
# res3120=dd(y[:,2],y[:,3]) 
# res3210=dd(y[:,2],y[:,1]) 
# res320=dd(y[:,2],y[:,4]) 
global res
res={}
res['2010']=dd(y[:,2],y[:,0]) 
res['3010']=dd(y[:,4],y[:,0]) 
# plt.plot(dtAx,res210)
# print res210
# plt.show()
# plt.plot(dtAx,res310)
# print res310
# plt.show()
# plt.plot(dtAx,res3120)
# print res3120
# plt.show()
d1=y[dthalfbin,0]
d3=y[dthalfbin,1]
printArray(0,True)
printArray(1,True)
printArray(2,True)
print "2_0 - 1_0 to+tp+d1+d2 -to-d1=tp+d2"
printDic('2010')
print pmtLab[3]
print "3_0 - 1_0 to+tp+d1+d2 +d3 -to-d1=tp+d2+d3"
printDic('3010')
lat=15.0
tp=lat/c
print "tp= "+str(tp)
d2=res['2010'][dthalfbin]-tp
print "d1 = "+str(d1)
print "d2 = "+str(d2)
print "d3 = "+str(d3)
# for i in xrange(0,6):
    # print pmtLab[i]
    # y=[]
    # dtL=dtHist[i]
    # print len(dtL)
    # for j in xrange(0, len(dtL)):
        # if len(dtL[j]) != 0:
            # y.append(np.array(dtL[i]).mean())
        # else :
            # y.append(0)
        # print "len y  "+str(len(y))
        # plt.plot(y,dtHist,'o')
        # plt.show()

# ############logRat###################
# mean=np.array(logRat).mean()
# std=np.array(logRat).std()
# print "log mean "+ str(mean)
# print "log std "+ str(std)
# print

# # n, bins, patches =plt.hist(np.array(logRat),bins=20)
# #plt.xlabel("Log Ratio (S0/S1)")
# plt.show()
# ############logRat vs Dt###################
# n=len(logRat)
# stdSum=0
# slope, intercept, r_value, p_value, std_err = stats.linregress(dtList,logRat)
# for i in xrange(0,len(logRat)):
    # ratEst=slope*dtList[i]+intercept
    # stdSum+=(logRat[i]-ratEst)**2
# stdErr=sqrt(stdSum/len(logRat))

# ran=10
# minT=-ran
# maxT=ran
# tRange=range(int(minT),int(maxT))
# logRange=[]
# for i in tRange:
    # logRange.append(slope*float(i)+intercept)

# print "slope "+str(slope)
# print "intercept "+str(intercept)
# print "r_value "+str(r_value)
# print "stdErr "+str(stdErr)
# #plt.plot(dtList,logRat,'.')
# #plt.errorbar(tRange, logRange, yerr=1*stdErr, fmt='o',color='g',ecolor='g',capthick=2)
# #plt.ylabel("log ratio of signals")
# #plt.xlabel("delta time(ns)")
# #plt.title("log10(S0/S1) vs Time")
# plt.show()
# print 
# ############Pulse###################
# mu , sigma=norm.fit(reduc)
# n, bins, patches =plt.hist(np.array(total),bins=100,normed=True)
# gaus = mlab.normpdf( bins, mu, sigma)
# plt.plot(bins, gaus, 'r--', linewidth=2)
# plt.xlabel("Pulse Integral")
# plt.show()

# print "pulse mean "+ str(np.array(total).mean())
# print "pulse std "+ str(np.array(total).std())
# print "mu "+str(mu)
# print "sigma "+str(sigma)
# print 

# ############cut###################
# print "cutting"
# delList=[]
# for i in xrange(0,len(logRat)):
    # rat=logRat[i]
    # t=dtList[i]
    # ratEst=slope*t+intercept
    # ratCut=abs(rat-ratEst)>1*stdErr
    # cutList=[ratCut]
    # for cut in cutList:
        # if cut  :
            # if i not in delList:
                # delList.append(i)

# print "deleted events "+str(len(delList))
# for i in reversed(delList):
    # del logRat[i]
    # del dtList[i]
    # del total[i]
    # del pulseListL[i]
    # del pulseListR[i]
# reduc=[]
# for tot in total:
    # if tot>240000 and tot < 310000 :
        # reduc.append(tot)
# print ##############################
# print ##############################
# ############dT###################
# print "number "+str(len(dtList))
# mean=np.array(dtList).mean()
# std=np.array(dtList).std()
# print "dt mean "+ str(mean)
# print "dt spread "+ str(sqrt(12)*std)
# print "dt std "+ str(std)
# print

# # n, bins, patches =plt.hist(np.array(dtList),bins=20)
# #plt.xlabel("delta time (ns)")
# plt.show()
# ############logRat###################
# mean=np.array(logRat).mean()
# std=np.array(logRat).std()
# print "log mean "+ str(mean)
# print "log std "+ str(std)
# print

# # n, bins, patches =plt.hist(np.array(logRat),bins=20)
# #plt.xlabel("Log Ratio (S0/S1)")
# plt.show()
# ############logRat vs Dt###################
# n=len(logRat)
# stdSum=0
# slope, intercept, r_value, p_value, std_err = stats.linregress(dtList,logRat)

# ran=10
# minT=-ran
# maxT=ran
# tRange=range(int(minT),int(maxT))
# logRange=[]
# for i in tRange:
    # logRange.append(slope*float(i)+intercept)

# print "slope "+str(slope)
# print "intercept "+str(intercept)
# print "r_value "+str(r_value)
# #plt.plot(dtList,logRat,'.')
# #plt.errorbar(tRange, logRange, yerr=1*stdErr, fmt='--',color='r',ecolor='g',capthick=2)
# #plt.ylabel("log ratio of signals")
# #plt.xlabel("delta time(ns)")
# #plt.title("log10(S0/S1) vs Time")
# plt.show()
# print 
# ############Pulse###################
# mu , sigma=norm.fit(reduc)
# gaus = mlab.normpdf( bins, mu, sigma)
# n, bins, patches =plt.hist(np.array(total),bins=100,normed=True)
# gaus = mlab.normpdf( bins, mu, sigma)
# plt.plot(bins, gaus, 'r--', linewidth=2)
# plt.xlabel("Pulse Integral")
# plt.show()

# print "pulse mean "+ str(np.array(total).mean())
# print "pulse std "+ str(np.array(total).std())
# print "mu "+str(mu)
# print "sigma "+str(sigma)
