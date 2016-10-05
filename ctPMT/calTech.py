import array
import numpy as np
import matplotlib.pyplot as plt
import rootpy.ROOT as ROOT
from peaks import peaks

p=peaks("total.root")
speMean=[]
speSigma=[]
mpeSpeMean=[]
spe_mpeRatio=[]
mpeMean=[]
mpeSigma=[]
PE=[]
mean_sig2=[]
for sn in p.snDict['spe']:
    meanS=p.fitDict["spe","mean",sn]
    sigS=p.fitDict["spe","sigma",sn]
    meanM=p.fitDict["mpe","mean",sn]
    rat=meanM/meanS
    speMean.append(meanS)
    mpeSpeMean.append(meanM)
    speSigma.append(sigS)
    spe_mpeRatio.append(rat)
sigSpeMean=np.array(speSigma).mean()
meanMpeSpeRatio= np.array(spe_mpeRatio).mean()
aveSpeMean=np.array(speMean).mean()
sig_spe=sigSpeMean/aveSpeMean
norm=sig_spe*sig_spe+1
print aveSpeMean
chek=[]
for sn in p.snDict['mpe']:
    sigMpe=p.fitDict["mpe","sigma",sn]
    meanMpe=p.fitDict["mpe","mean",sn]
    rat=sigMpe/(aveSpeMean*norm)
    PE.append(rat*rat)
    rat2=meanMpe/sigMpe
    mean_sig2.append(rat2*rat2)
    ratio=rat/rat2
    chek.append(ratio*ratio)
p.plotHist("all",320)
raw_input()
#spe=[]
#print p.snDict["spe"]
#for sn in p.snDict["spe"]:
    #print str(sn)+" "+str(p.fitDict["spe","mean",sn])
    #spe.append(p.fitDict["spe","mean",sn])
#plt.hist(spe)
#plt.show()
#print np.array(PE).mean()
#print np.array(mean_sig2).mean()
#print np.array(chek).mean()
#plt.hist(mean_sig2)
#plt.show()
#plt.hist(chek)
#plt.show()
#p.plotParameter("sig2","sig2")
#p.plotParameter("mpe/spe","mpe/spe")
#p.plotParameter("spe","mean")
#p.plotParameter("spe","sigma")
#p.plotParameter("spe","ratio")
#p.plotParameter("mpe","mean")
#p.plotParameter("mpe","sigma")
#p.plotParameter("best","pe")
#p.plotParameter("mpe","ratio")
#p.plotParameter("back","mean")
#p.plotParameter("back","sigma")
#p.plotParameter("back","ratio")
############
