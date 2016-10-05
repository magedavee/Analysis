from initROOT import initROOT
import ROOT
from ROOT import gROOT, TCanvas, TF1
import numpy as np
import matplotlib.pyplot as plt
initROOT()
# pmt=ROOT.PmtData("/home/mage/Data/p20Data/root/P20_2015-04-16-10-44-25.root")
pmt=ROOT.PmtData("/home/mage/Data/proc_cry_0.root.pro")
event=0
ent= pmt.GetEntries() 
pmt.CalIntegral(0)
pmt.CalIntegral(1)
integral0=[]
integral1=[]
pulseList=[]
dtList=[]
for i in xrange(0,ent):
    pmt.SetEntry(i)
    i0= pmt.GetPulseIntegral(0,i)
    i1= pmt.GetPulseIntegral(1,i)
    if pmt.GetNCha() >=4 and i0 >10000:
	ch1=pmt.GetPulse(0)
	ch2=pmt.GetPulse(1)
	peak1=np.array(ch1).min()
	peak2=np.array(ch2).min()
	if  peak1 < - 50 or peak2<-50:
	    # pmt.GetTrace(0).Draw()
	    # pmt.GetTrace(1).Draw("SAME")
	    # input("pause ")
	    tr1=pmt.GetTrace(0)
	    tr2=pmt.GetTrace(1)
	    t1=0
	    t2=0
	    for i in xrange(0,1100):
		if tr1.Eval(i) < .5*peak1:
		    t1=i
		    break
	    for i in xrange(0,1200):
		if tr2.Eval(i) < .5*peak2:
		    t2=i
		    break
	    dt=t1-t2
	    integral0.append(i0)
	    integral1.append(i1)
	    if abs(dt)<50:
		dtList.append(dt)
	    pulse=np.array(pmt.GetPulse(0))
	    pulseList.append(pulse)
	# plt.plot(pulse)
	# plt.show()
print len(dtList)
print np.array(dtList).mean()
print np.array(dtList).var()
plt.hist(np.array(dtList),bins=50)
plt.xlabel("Delta Time (ns)")
plt.show()
# pmt.GetIntegral(0).Draw()
# input()
# pmt.GetIntegral(1).Draw()
# input()
    # pmt.SetEntry(i)
    # print i
    # print pmt.GetNCha()
    # if pmt.GetNCha() >= 4:
	# pmt.GetTrace(0).Draw()
	# pmt.GetTrace(1).Draw("SAME")
	# input("pause ")
# while event >=0:
    # pmt.GetTrace(0).Draw()
    # # input("pause ")
    # pmt.GetTrace(1).SetLineColor(4);
    # pmt.GetTrace(1).Draw("SAME")
    # event=int(raw_input("promt: "))
    # pmt.SetEntry(event)
#############
