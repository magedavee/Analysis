from funcAna import *
from sql_functions import *
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
totList=[]
qTopList=[]
qBotList=[]
reTotList=[]
reQTopList=[]
reQBotList=[]
pulse=[]
pmt=ROOT.PmtData(name,"p50_s4_setup.root")
ent=pmt.GetEntries()
tanl=[]
cosl=[]
sinl=[]
# tant=[]
theta=[]
tCheck=[]
dAnti=[]
for i in xrange(6):
    dAnti.append([])
    for j in xrange(6):
        dAnti[i].append([])
print "what it should be "+str(ptime)
for i in xrange(0,ent):
    pmt.SetEntry(i)
    printPercent(i,ent)
    q,t,lr=getPara(pmt,i)
    fixGain(q,corr)
    tot,totT,totB=getTots(q)
    if cut(t,q) :
        # fixTime(t,mean,ptime=ptime)
        fixTime(t,mean)
        tt=t[0]-t[1]
        ddt.append(tt)
        totList.append(tot)
        qTopList.append(totT)
        qBotList.append(totB)
        l=sqrt((M*tt)**2+(lat)**2)
        t0=50.0/M-t[0]
        t2=50.0/M-t[1]
        tcheck=t[2]-l/c
        tCheck.append(tcheck)
        lenCal.append(l)
        tanl.append(M*tt/lat)
        sinl.append(M*tt/l)
        cosl.append(lat/l)
        theta.append((180.0/pi)*np.arcsin(M*tt/l))
        # tant.append(t[2]/tt)
        for j in xrange(6):
            for k in xrange(6):
                dAnti[j][k].append(t[j]-t[k])
        for j,d in enumerate(dt):
            d.append(t[j])
            logRat[j].append(lr[j])
        for j in xrange(0,len(pulse)):
            pulse[j].append(q[j])
print "Analizing"
###############Analysis Setup #########
name=['cell_1','cell_2','left','0_3cross','right ','1_2cross']
# printHist(dt,logRat)
tanl=np.array(tanl)
cosl=np.array(cosl)
sinl=np.array(sinl)
theta=np.array(theta)
tCheck=np.array(tCheck)
# tant=np.array(tant)
totList=np.array(totList)
qTopList=np.array(qTopList)
qBotList=np.array(qBotList)
tCheck=np.array(tCheck)
for i,l in enumerate(lenCal):
    reTotList.append(totList[i]/l)
    reQTopList.append(qTopList[i]/l)
    reQBotList.append(qBotList[i]/l)
for i,l in enumerate(lenCal):
    qt=qTopList[i]
    qb=qBotList[i]
    rqt=reQTopList[i]
    rqb=reQBotList[i]
    t=ddt[i]
    # insertEvent(i,qt,qb,rqt,rqb,l,t)
reTotList=np.array(reTotList)
reQTopList=np.array(reQTopList)
reQBotList=np.array(reQBotList)
###############Analysis #########
print '##############Hist#############'
plotHist(tanl,name="Hist tan")
plotHist(cosl,name="Hist cos")
plotHist(sinl,name="Hist sin",xlabel='Sin theta')
plotHist(theta,name="Hist theta",xlabel='theta')
plotHist(tCheck,name="Hist tCheck")
# plotHist(tant,name="Hist tant")
plotHist(lenCal,name="Hist cal lenght",xlabel='Track Length(cm)')
plotHist(reTotList,name='Pulse integral over length',xlabel='dE/dx (ch/MeV)')
plotHist(totList,name='Pulse integral')
plotHist(reQTopList,name='Hist tot over length')
plotHist(reQBotList,name='Hist tot over length')
plotHist(ddt,name='Hist delta dt',xlabel='dt1 - dt2 (ns)')
for i,d in enumerate(dt):
    plotHist(d,name='dt '+name[i],xlabel='delta time (ns)')
for i,d in enumerate(logRat):
    plotHist(d,name="log ratio "+name[i],xlabel='log ratio log(S1.S0)')
print '##############Reg#############'
for i in xrange(len(dt)):
    for j in xrange(i+1,len(dt)):
        plotReg(dt[i],dt[j],'Reg '+name[i]+"_"+name[j])
plotReg(dt[3],ddt,name[i]+" vs ddt",xlabel='dt (ns)',ylabel='ddt (ns)')
plotReg(dt[3],tanl,name[i]+"vs tan",xlabel='dt (ns)',ylabel='tan')
plotReg(dt[3],cosl,name[i]+"vs cos",xlabel='dt (ns)',ylabel='tan')
plotReg(dt[3],sinl,name[i]+"vs sin",xlabel='dt (ns)',ylabel='tan')
plotReg(dt[3],theta,name[i]+"vs theta",xlabel='dt (ns)',ylabel='tan')
for i,d in enumerate(dt):
    plotReg(d,logRat[i],name[i]+"log ratio vs dt",xlabel="dt (ns)",ylabel='log ratio')
for i in xrange(len(dt)):
    plotReg(dt[i],ddt,'Reg '+name[i]+"_ddt")
    plotReg(dt[i],tanl,'Reg '+name[i]+"_tanl")
    plotReg(dt[i],cosl,'Reg '+name[i]+"_cosl")
    plotReg(dt[i],sinl,'Reg '+name[i]+"_sinl")
    plotReg(dt[i],theta,'Reg '+name[i]+"_theta")
    plotReg(dt[i],tCheck,'Reg '+name[i]+"_tCheck")
    plotReg(lenCal,dt[i],name="Reg lenCal vs "+name[i])
    plotReg(logRat[i],ddt,'Reg '+name[i]+"log_ddt")
    plotReg(logRat[i],tanl,'Reg '+name[i]+"log_tanl")
    plotReg(logRat[i],cosl,'Reg '+name[i]+"log_cosl")
    plotReg(logRat[i],sinl,'Reg '+name[i]+"log_sinl")
    plotReg(logRat[i],theta,'Reg '+name[i]+"log_theta")
    plotReg(logRat[i],tCheck,'Reg '+name[i]+"log_tCheck")
    plotReg(lenCal,logRat[i],name="Reg lenCal vs log"+name[i])
    # plotReg(logRat[i],tant,'Reg '+name[i]+"_tant")
plotReg(lenCal,qTopList,name='Reg renormQtop')
plotReg(lenCal,qBotList,name='Reg renormQbot')
plotReg(totList,ddt,name="Reg tot vs ddt")
plotReg(logRat[0],totList,name="Reg tot vs logRat ")
plotReg(tanl,ddt,name="Reg tanl vs ddt")
plotReg(cosl,ddt,name="Reg cosl vs ddt")
plotReg(sinl,ddt,name="Reg sinl vs ddt")
plotReg(tCheck,ddt,name="Reg sinl vs ddt")
plotReg(theta,ddt,name="Reg theta vs ddt")
plotReg(tCheck,theta,name='Reg tChech vs theta')
plotReg(tCheck,cosl,name='Reg tChech vs cos')
plotReg(tCheck,sinl,name='Reg tChech vs sin')
plotReg(tCheck,tanl,name='Reg tChech vs tan')
plotReg(tCheck,lenCal,name='Reg tChech vs len')
res=[]
timBin=[]
time=[]
print '##############Slice#############'
plotSlice(qTopList,dt[0],name="Slice_q_top")
plotSlice(qBotList,dt[0],name="Slice_q_bot")
plotSlice(totList,dt[0],name="Slice_q_tot")
plotSlice(reQTopList,dt[0],name="Slice_re_q_top")
plotSlice(reQBotList,dt[0],name="Slice_re_q_bot")
plotSlice(reTotList,dt[0],name="Slice_re_top")
plotSlice(qTopList,logRat[0],name="log_Slice_q_top ",bins=21,start=-1,size=.1)
plotSlice(qBotList,logRat[0],name="log_Slice_q_bot",bins=21,start=-1,size=.1)
plotSlice(totList,logRat[0],name="log_Slice_q_tot",bins=21,start=-1,size=.1)
plotSlice(reQTopList,logRat[0],name="log_Slice_re_q_top",bins=21,start=-1,size=.1)
plotSlice(reQBotList,logRat[0],name="log_Slice_re_q_bot",bins=21,start=-1,size=.1)
plotSlice(reTotList,logRat[0],name="log_Slice_re_top",bins=21,start=-1,size=.1)
plot2dHist(qTopList,dt[0], name="2D_q_top_dt")

for i in xrange(6):
    for j in xrange(i+1,6):
        print'##############'
        n=name[i] +" vs "+name[j]
        res.append(plotReg(dt[i],dt[j],n))
        n=name[j] +" vs "+name[i]
        res.append(plotReg(dt[j],dt[i],n))
for i in xrange(len(dt)) :
    print'##############'
    res.append(plotReg(ddt,dt[i],"ddt vs "+name[i]))
    res.append(plotReg(dt[i],ddt,name[i]+" vs ddt"))
print "number of events "+str(len(ddt))
