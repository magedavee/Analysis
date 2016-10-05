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
##################plot regession of two var###########
def plotReg(lt1,lt2,name='plot',xlabel=None,ylabel=None):
    print name
    if type(lt1) is list and type(lt2) is list :
        lt1=np.array(lt1)
        lt2=np.array(lt2)
    else :
        if type(lt1) is not np.ndarray and type(lt2) is not np.ndarray:
            return 
    m, b, r_value, p_value, std_err = stats.linregress(lt1,lt2)
    print "m:  "+str(m)
    print "b:  "+str(b)
    print "r:  "+str(r_value)
    print "p:  "+str(p_value)
    print "e:  "+str(std_err)
    plt.plot(lt1,lt2,'.')
    if xlabel is not None :
        plt.xlabel(xlabel)
    if ylabel is not None :
        plt.ylabel(ylabel)
    plt.title(name)
    plt.savefig("./figures/"+name,bbox_inches='tight')
    plt.clf()
    return [m,b,r_value,p_value,std_err]
##################plot regession of hist of varible###########
def plotHist(l,name,bins=100,xlabel=None,ylabel=None):
    if type(l) is not np.ndarray:
        if type(l)==list :
            l=np.array(l)
        else :
            return
    print name+" "+str(l.mean())+" "+str(l.std())
    plt.hist(l,bins)
    if xlabel is not None :
        plt.xlabel(xlabel)
    if ylabel is not None :
        plt.ylabel(ylabel)
    plt.title(name)
    plt.savefig("./figures/"+name,bbox_inches='tight')
    plt.clf()
    return l.mean(),l.std()
##################plot two var by binning on var###########
def plotSlice(pro,sel,name="",bins=61,start=-30,size=1,xlabel=None,ylabel=None):
    print name
    x=[]
    xBin=[]
    for i in xrange(bins):
        xBin.append([])
    for i,p in enumerate(pro):
        b=int((sel[i]-start)/size)
        if b <bins and b > 0 :
            xBin[b].append(pro[i])
    proAve=[]
    yerr=[]
    for i in xrange(bins):
        num=len(xBin[i])
        npx=np.array(xBin[i])
        proAve.append(npx.mean())
        x.append(i*size+start)
        if num !=0 :
            error=3.0*npx.std()/sqrt(num)
            yerr.append(error)
        else :
            yerr.append(0)
    # plt.figure()
    plt.errorbar(x,proAve,yerr=yerr)
    if xlabel is not None :
        plt.xlabel(xlabel)
    if ylabel is not None :
        plt.ylabel(ylabel)
    plt.title(name)
    plt.savefig("./figures/"+name,bbox_inches='tight')
    plt.clf()
    # plt.show()
    return proAve

def plot2dHist(x,y,name='',norm=None,bins=40,xlabel=None,ylabel=None):
    plt.hist2d(y, x, bins=bins)
    if xlabel is not None :
        plt.xlabel(xlabel)
    if xlabel is not None :
        plt.ylabel(ylabel)
    plt.colorbar()
    plt.title(name)
    plt.savefig("./figures/"+name,bbox_inches='tight')
    plt.clf()


##############check if list###############
def checkList(l):
    if type(l) is not np.ndarray:
        if type(l)==list :
            return np.array(l)
        else :
            print "Opps not a list or np.ndarray"
            return
    else :
        return l
    
##############percent print###############
def printPercent(i,ent) :
    inval=int(.1*ent)
    if i%inval==0:
        print str(float(int(i*1000/ent))/10.0)+" %"
############## pickle file io###############

def readPickle(name):
    if type(name) is not str :
        print 'whats name?'
        return
    print "processing "+name
    data=[]
    infile = open(name, 'r')
    while True :
        try:
            p=pickle.load(infile)
            data.append(p)
        except EOFError :
            break
    infile.close()
    return data

def writePickle(name,arr):
    print "Now creating "+name
    f=open(name,"wb")
    for c in arr :
        print str(c)+" is being written"
        pickle.dump(c,f)
    f.close()
##############cut###############
def cut (t,q):
    tot=0
    totT=0
    totB=0
    for i,w in enumerate(q):
        tot+=w
        if i <2 :
            totT+=w
        else:
            totB+=w
    lb=40000
    if abs(t[2])> 30 or abs(t[3])>30 or abs(t[4])> 30 or abs(t[5])>30 :
        return False

    if abs(t[0])> 30 or abs(t[1])>30 :
        return False
    if q[0]==0 or q[2]==0:
        cond = False
    if not ((q[0]> lb and q[2]> lb) and  (q[3]> lb and  q[1]>lb)) :
        return False
    
    if  tot <1000000:
        return False
    return True
##############Sets varibles###############
def getTots(q) :

    tot=0
    totT=0
    totB=0
    for i,w in enumerate(q):
        tot+=w
        if i <2 :
            totT+=w
        else:
            totB+=w
    return [tot,totT,totB]
def getPara(pmt,i):
    
    t=[0,0,0,0,0,0]
    q=[0,0,0,0]
    lr=[0,0,0,0,0,0]
    t[0]=pmt.DeltaT(0,1)
    t[1]=pmt.DeltaT(2,3)
    t[2]=pmt.DeltaT(0,2)
    t[3]=pmt.DeltaT(0,3)
    t[4]=pmt.DeltaT(1,3)
    t[5]=pmt.DeltaT(1,2)
    for j in xrange(0,len(q)):
        q[j]=float(pmt.GetPulseIntegral(j,i))
    lr[0]=log10(q[0]/q[1])
    lr[1]=log10(q[2]/q[3])
    lr[2]=log10(q[0]/q[2])
    lr[3]=log10(q[0]/q[3])
    lr[4]=log10(q[1]/q[3])
    lr[5]=log10(q[1]/q[2])
    return q,t,lr

##############print hist info plots###############
def printHist(dt=None,logRat=None, pulse=None):
    print '#################print hist'
    name=['cell_1','cell_2','left','0_3cross','right ','1_2cross']
    if  dt is not None :
        dt=checkList(dt)
    if  logRat is not None :
        logRat=checkList(logRat)
    if  pulse is not None :
        pulse=checkList(pulse)
    for i in xrange (0,len(dt)):
        if  dt is not None:
            plotHist(dt[i],"Hist dt_"+name[i])
        if  logRat is not None :
            plotHist(logRat[i],"Hist logRat_"+name[i])
        if  logRat is not None and dt is not None :
            plotReg(logRat[i],dt[i],"Reg log_vs_dt_"+name[i])
        if  pulse is not None :
            if i<len(pulse):
                print pulse[i].mean()

##############Corrections###############
def fixTime(t,mean,ptime=0) :
    if type(t) is not list :
        print "not at list"
        return
    t[0]-=mean[0]
    t[1]-=mean[1]
    t[2]+=(ptime-mean[2])
    t[3]+=(ptime-mean[2]-mean[1])
    t[4]+=(ptime-mean[2]-mean[1]+mean[0])
    t[5]+=mean[0]+ptime-mean[2]
    # t[2]+=(ptime-mean[3])
    # t[3]+=(ptime-mean[3]-mean[1])
    # t[4]+=(ptime-mean[3]-mean[1]+mean[0])
    # t[5]+=mean[0]+ptime-mean[3]

def fixGain(q,corr) :
    for j in xrange(len(q)):
        q[j]=corr[j]*q[j]
