import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import ROOT
from ROOT import gROOT
from math import *
#this script sums over the cells
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-build/lib/libEventLib.so")
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-Sec/include/Output/Event.hh")
histCell=ROOT.TH2D("Cell Ionization Hits","Cell Ionization Hits",14,0,14,10,0,10)
sceFile=ROOT.TFile("TotalPhoton.root")
tree=ROOT.TTree()
tree=sceFile.Get("Muon Calibration")
Nx=14.0
Ny=10.0
lat=144.198
c=3*10**2
X=np.zeros(1,dtype=float)
Y=np.zeros(1,dtype=float)
Z=np.zeros(1,dtype=float)
L=np.zeros(1,dtype=float)
E=np.zeros(1,dtype=float)
YVar=np.zeros(1,dtype=float)
LVar=np.zeros(1,dtype=float)
RVar=np.zeros(1,dtype=float)
Dt=np.zeros(1,dtype=float)
Length=np.zeros(1,dtype=float)
Right=np.zeros(1,dtype=int)
Left=np.zeros(1,dtype=int)
Vol=np.zeros(1,dtype=int)
Evt=np.zeros(1,dtype=int)
p50=[]
tree.SetBranchAddress("x",X)
tree.SetBranchAddress("y",Y)
tree.SetBranchAddress("z",Z)
tree.SetBranchAddress("yVar",YVar)
tree.SetBranchAddress("dt",Dt)
tree.SetBranchAddress("right",Right)
tree.SetBranchAddress("leftt",Left)
tree.SetBranchAddress("vol",Vol)
tree.SetBranchAddress("evt",Evt)
tree.SetBranchAddress("E",E)
#geant 
tree.SetBranchAddress("length",Length)
#py
tree.SetBranchAddress("l",L)
tree.SetBranchAddress("lVar",LVar)
tree.SetBranchAddress("rVar",RVar)
entries=int(tree.GetEntries())
dtList=[]
xList=[]
yList=[]
yyList=[]
zList=[]
ySDList=[]
ySDList=[]
eList=[]
estEList=[]
geoLen=[]
geantLen=[]
totalList=[]
lSDList=[]
rSDList=[]
lVarList=[]
rVarList=[]
y_SDList=[]
dSDList=[]
ratList=[]
nxList=[]
nyList=[]
evtList=[]
volList=[]
skip=[]
evtMem=-1
evtCount=0
throwout=0
simLen=[]
fulLen=[]
#chech to make sure events are numbered correclty
for i in xrange(0,entries):
    tree.GetEntry(i)
    evt=Evt[0]
    if evt!=evtMem :
	if evtCount<10 and evtCount>0:
	    skip.append(evtMem)
	evtMem=evt
	evtCount=1
    else:
	evtCount+=1

#copy good events to new buffers
for i in xrange(0,entries):
    tree.GetEntry(i)
    if Evt[0] not in skip:
	# if E[0]>10**-100 or True:
	lSDList.append(LVar[0])
	rSDList.append(RVar[0])
	rVarList.append(RVar[0]**2)
	lVarList.append(LVar[0]**2)
	dSDList.append(LVar[0]-RVar[0])
	vol=Vol[0]
	evtList.append(Evt[0])
	volList.append(Vol[0])
	nx=int(vol)%int(Nx)
	ny= int(vol)/int(Nx)
	histCell.Fill(nx,ny)
	# print "event: "+str(Evt[0])+" nx:"+ str(nx)+ " ny:"+ str(ny)+"rec: "+str(nx+ny*Nx)+" vol: "+str(vol)
	nxList.append(nx)
	nyList.append(ny)
	dt=Dt[0]
	right=Right[0]
	left=Left[0]
	x=X[0]
	y=Y[0]
	z=Z[0]
	yVar=YVar[0]
	ySDList.append(yVar)
	ySDList.append(yVar**2)
	y_SDList.append(1.3*(y+500)/(c*RVar[0]))
	xList.append(x)
	yList.append(y)
	#yyList.append(sqrt(y+500))
	zList.append(z)
	# print z
	dtList.append(dt)
	total=int(left)+int(right)
	totalList.append(total)
        if ny%2 == 0 :
            p50.append(0)
        else :
            p50.append(1)
	# print str(i)+" total: "+str(L[0])+" right: "+str(Length[0])
	if total!=0:
	    ratio=float(right)/float(total)
	    ratList.append(ratio)
	eList.append(E[0])
	geoLen.append(L[0])
	geantLen.append(Length[0])
    else:
	throwout+=1
# p=np.polyfit(np.array(lSDList),np.array(yList),2)
m, b, r_value, p_value, std_err = stats.linregress(np.array(dtList),np.array(yList))
# print np.array(dEdx).mean()

yphotonList=[]
for i in xrange(0,len(dtList)):
    yphotonList.append(m*dtList[i]+b)
evtMem=-1
x0=0
y0=0
z0=0
x1=0
y1=0
z1=0
estList=[]
errList=[]
lenTotList=[]
truLenList=[]
measureDeDx=[]
trueDeDx=[]
dyList=[]
thetaList=[]
thetaList2=[]
#calculate length energy ect
for i in xrange(0,len(evtList)):
    if evtMem!=evtList[i]: 
	if evtMem>-1:
	    ny1=0
	    dxRecon=(x0-x1)
	    dyRecon=(y0-y1)
	    dzRecon=(z0-z1)
	    dyList.append(dyRecon)
	    theta=atan(dyRecon/dzRecon)
	    # theta=atan(dyRecon/(Ny*lat))
	    # print "theta: "+str(theta)
	    thetaList.append(theta)
	    estLen=sqrt((dxRecon)**2+(dyRecon)**2+(dzRecon)**2)/10.0
	    # print "   dx: "+str(x0-x1)+" dy: "+str(y0-y1)+" dz: "+str(z0-z1)
	    dxTrue=(pos0[0]-pos1[0])
	    dyTrue=(pos0[1]-pos1[1])
	    dzTrue=(pos0[2]-pos1[2])
	    truLen=sqrt((dxTrue)**2+(dyTrue)**2+(dzTrue)**2)/10.0
	    # print "true: "+str(truLen)
	    # print "   dx: "+str(pos0[0]-pos1[0])+" dy: "+str(pos0[0]-pos1[0])+" dz: "+str(pos0[0]-pos1[0])
	    truLenList.append(truLen)
	    estList.append(estLen)
            fulLen.append(geantCount)
            simLen.append(geoCount)
	    # print "geantCount: "+str(geantCount)
	    # print "geoCount: "+str(geoCount)
	    # print evtList[i]
	    # print "x0: "+str(x0)+" y0: "+str(y0)+" z0: "+str(z0)+" ny: "+str(nyList[i])
	    # print "x1: "+str(x1)+" y1: "+str(y1)+" z1: "+str(z1)+" ny: "+str(ny1)
	    # print "pos0 x: "+str(pos0[0])+" pos0 y: "+str(pos0[1])+" pos0 z: "+str(pos0[2])
	    # print "pos1 x: "+str(pos1[0])+" pos1 y: "+str(pos1[1])+" pos1 z: "+str(pos1[2])
	    # print "dx recon: "+str(dxRecon)+" dy: "+str(dyRecon)+" dz: "+str(dzRecon)
	    # print "dx true: "+str(dxTrue)+" dy: "+str(dyTrue)+" dz: "+str(dzTrue)
	    # print "recon: "+str(estLen)
	    # print "true: "+str(truLen)
	    # print "lenth: "+str(geoCount)
	    # # print "x: "+str(float(nxList[i])*lat)+" y: "+str(yphotonList[i])+" z: "+str( float(nyList[i])*lat)
	    # print ""
    

	geantCount=geantLen[i]
	geoCount=geoLen[i]
	evtMem=evtList[i]
	x0=float(nxList[i]-Nx/2.0)*lat+lat/2
	z0=float(nyList[i])*lat-(lat*Ny)/2.0
	y0=yphotonList[i]
	pos0=[xList[i],yList[i],zList[i]]
    else:
	x1=float(nxList[i]-Nx/2.0)*lat+lat/2
	z1=float(nyList[i])*lat-(lat*Ny)/2.0
	ny1=nyList[i]
	y1=yphotonList[i]
	geantCount+=geantLen[i]
	geoCount+=geoLen[i]
	pos1=[xList[i],yList[i],zList[i]]
	# print evtList[i]
dxRecon=(x0-x1)
dyRecon=(y0-y1)
dzRecon=(z0-z1)
dyList.append(dyRecon)
theta=atan(dyRecon/lat)
estLen=sqrt((dxRecon)**2+(dyRecon)**2+(dzRecon)**2)/10.0
estList.append(estLen)
dxTrue=(pos0[0]-pos1[0])
dyTrue=(pos0[1]-pos1[1])
dzTrue=(pos0[2]-pos1[2])
truLen=sqrt((dxTrue)**2+(dyTrue)**2+(dzTrue)**2)/10.0
truLenList.append(truLen)

thetaList.append(theta)
for i in thetaList:
    for j in xrange(0,10):
	thetaList2.append(i)
lenSingleList=[]
measuredDeDxSingle=[]
phoPerLenListEvent=[]
phoPerLenList=[]
eTotList=[]
phoTotList=[]
# for i in xrange(0,len(thetaList)):
    # theta=thetaList[i]
    # for j in xrange(0,10):
	# print str(evtList[j+i*10])+" "+str(nyList[j+i*10])
    # print " "

for i in xrange(0,len(thetaList)):
    theta=thetaList[i]
    eTot=0
    phoTot=0
    l=lat/(10*cos(theta))
    for j in xrange(0,10):
	e=eList[j+i*10]
	eTot+=e
	lenSingleList.append(l)
	measuredDeDxSingle.append(e/l)
	photon=totalList[j+i*10]
	phoTot+=photon
	phoPerLenList.append(photon/l)
	# print "length: "+str(l)
	# print "photon: "+str(photon)
    # estList[i]+=l
    # truLenList[i]+=l
    # print estList[i]/l
    eTotList.append(eTot)
    phoTotList.append(phoTot)
    phoPerLenListEvent.append(phoTot/(estList[i]))
    measureDeDx.append(eTot/estList[i])
    trueDeDx.append(eTot/truLenList[i])
    err=(float(estList[i])-(truLenList[i]))/float(truLenList[i])
    errList.append(err)
    # phoPerLenList.append(phoTot/(l*10))
    # print eTot
    # print phoTot
    # print ""

eErr=[]
for i in estList:
    eest=1.916*i
    estEList.append(eest)
for i in xrange(len(estEList)):
        err=100.0*((estEList[i]-eTotList[i])/eTotList[i])
        eErr.append(err)


print "length error:"+str(np.array(errList).mean())
print "recon dedx event"+str(np.array(measureDeDx).mean())
print "dedx using true length"+str(np.array(trueDeDx).mean())
print "recon dedx Single:"+str(np.array(measuredDeDxSingle).mean())
print "y/tvar: "+str(np.array(y_SDList).mean())
print "slope "+str(m)
print "intercept "+str(b)
print "r value :"+str(r_value)
print "p value :"+str(p_value)
print "std error :"+str(std_err)


########plot regresion delta time y pos###############
# fig=plt.figure()
# ax = fig.add_subplot(111)
# plt.plot(np.array(dtList),np.array(yList),'.')
# equ=[m,b]
# x=np.arange(np.array(dtList).min(),np.array(dtList).max())
# eq="best fit:"+str(round(m,2))+"deltaTime+"+str(round(b,2))
# plt.plot(x,m*x+b,'-')
# plt.title("y postion VS delta time")
# plt.ylabel("y postion(mm)")
# plt.xlabel("delta time(ns)")
# ax.text(2, 12,eq, fontsize=12)

########delta time################ 
print "Delta time"
plt.hist(dtList,bins=100)
plt.title("Delta T(ns)")
plt.show()
print " mean: "+str(np.array(dtList).mean())
print " std: "+str(np.array(dtList).std())
#######ypos################ 
print "Delta time"
plt.hist(yList,bins=100)
plt.title("Postion Along Cell (cm)")
plt.show()
print " mean: "+str(np.array(yList).mean())
print " std: "+str(np.array(yList).std())
###############dt var left vs len ########
#x=np.arange(np.array(lSDList).min(),np.array(lSDList).max(),.001)
#plt.plot(np.array(lSDList),np.array(yList),'.')
#plt.plot(x,p[0]*x**2+p[1]*x+p[2],"red")
#eq=str(p[0])+"x^2+"+str(p[1])+"x"+str(p[2])
#plt.text(4,-600,eq,fontsize=12)
#plt.ylabel("y postion(mm)")
#plt.xlabel("delta time standard deviation(ns)")
#plt.title("position along the cell vs spread of photon times")

###############dt var right vs len ########
#plt.plot(np.array(rSDList),np.array(yList),'.')
#plt.ylabel("y postion(mm)")
#plt.xlabel("delta time standard deviation(ns)")
#plt.title("position along the cell vs spread of photon times")

###############dt var right vs len ########
#plt.plot(np.array(rSDList),np.array(lSDList),'.')
#plt.ylabel("y postion(mm)")
#plt.xlabel("delta time standard deviation(ns)")
#plt.title("position along the cell vs spread of photon times")

# ###########plot error in length###############
# print
# print "error lenght"
# plt.hist(errList,bins=100)
# plt.title("% difference of true length and reconstructed length")
# plt.show()
# print "errList mean: "+str(np.array(errList).mean())
# print "errList std: "+str(np.array(errList).std())
# # ##############plot recon DeDx#####################
# print
# print "Measured DeDx"
# plt.hist(measureDeDx,bins=100,range=[0,4])
# plt.title("measured dE/dx Event")
# print "dE/dx mean: "+str(np.array(measureDeDx).mean())
# print "dE/dx std: "+str(np.array(measureDeDx).std())
# plt.show()
# # # ##################True length DeDX################
# print
# print "true DeDc"
# plt.hist(trueDeDx,bins=100,range=[0,4])
# plt.title("true dE/dx")
# print "true dE/dx mean: "+str(np.array(trueDeDx).mean())
# print "true de/dx std: "+str(np.array(trueDeDx).std())
# perDiff=(np.array(measureDeDx).mean()-np.array(trueDeDx).mean())/np.array(trueDeDx).mean()
# print "% difference "+str(perDiff) 
# plt.show()

# # ##################CalEnergy################
# print
# print "CalEnergy"
# plt.hist(estEList,bins=60,range=[240,300])
# plt.title("Measured Energy")
# print "est Energy: "+str(np.array(estEList).mean())
# print "est std: "+str(np.array(estEList).std())
# plt.show()
# # ##################TrueEnergy################
# print
# print "True energy"
# plt.hist(eTotList,bins=60,range=[240,300])
# plt.title("True Energy")
# print "true Energy: "+str(np.array(eTotList).mean())
# print "true std: "+str(np.array(eTotList).std())
# # plt.show()
# # ##################Energy Error################
# print
# print "Energy Error"
# plt.hist(eErr,bins=100)
# plt.title("% Error")
# eEr=100*(np.array(estEList).mean()-np.array(eTotList).mean())/np.array(eTotList).mean()
# print "error Energy: "+str(np.array(eErr).mean())
# print "std: "+str(np.array(eErr).std())
# print "cali error "+str(eEr)
# # plt.show()
# ##################ful Length ################
# print
# print "ful Length     "
# plt.hist(fulLen,bins=100)
# plt.title("Lenght according to Geant4")
# print "full len: "+str(np.array(fulLen).mean())
# print "std: "+str(np.array(fulLen).std())
# plt.show()
# ##################sim Length ################
# print
# print "sim Length     "
# plt.hist(simLen,bins=100)
# plt.title("Lenght as piecswise")
# print "full len: "+str(np.array(simLen).mean())
# print "std: "+str(np.array(simLen).std())
# plt.show()
# ##################true Length ################
# print
# print "true Length     "
# plt.hist(truLenList,bins=100)
# plt.title("True Lenght")
# print "True len: "+str(np.array(truLenList).mean())
# print "std: "+str(np.array(truLenList).std())
# plt.show()
# ##################est Length ################
# print
# print "est  Length     "
# plt.hist(estList,bins=100)
# plt.title("Calculated Length")
# print "est len: "+str(np.array(estList).mean())
# print "std: "+str(np.array(estList).std())
# plt.show()
##############root plot of cells hit ###########################
# histCell.Draw("colz")

#######################angle distro#############################
# plt.hist(thetaList,bins=50,range=[0,pi/2])
# plt.title("angle of decent")

################Scatter length vs energy Single##############
# plt.plot(np.array(lenSingleList),np.array(eList),'.')
# plt.title("Photon count vs cell path length")
# plt.xlabel("Photon count")
# plt.ylabel("Length(cm)")
################Scatter length vs Energy Event##############
# plt.plot(np.array(estList),np.array(eTotList),'.')
# plt.title("Photon count vs cell path length")
# plt.xlabel("Photon count")
# plt.ylabel("Length(cm)")
################Scatter length vs photon Single##############
# plt.plot(np.array(lenSingleList),np.array(totalList),'.')
# plt.title("Photon count vs cell path length")
# plt.xlabel("Photon count")
# plt.ylabel("Length(cm)")

# ###############Scatter length vs photon in event##############
# plt.plot(np.array(estList),np.array(phoTotList),'.')
# plt.title("Photon count vs event path length")
# plt.xlabel("Photon count")
# plt.ylabel("Length(cm)")

################Scatter energy vs photon in cell################
#plt.plot(np.array(totalList),np.array(eList),'.')
#plt.title("Photon count vs Deposited Energy")
#plt.xlabel("Photon count")
#plt.ylabel("Energy deposted into a single cell(MeV)")

################scatter energy vs photon in event##############
#plt.plot(np.array(phoTotList),np.array(ETotList),'.')
#plt.title("photon count vs deposited energy")
#plt.xlabel("photon count")
#plt.ylabel("Energy deposted entire event(MeV)")

################photo in pmt#########
# plt.hist(totalList,bins=100)
# plt.xlabel("Photon count")
# plt.title("photon detected by pmt")


##############plot recon DeDx single#####################
# plt.hist(measuredDeDxSingle,bins=100)
# plt.xlabel("dE/dx single cell (MeV/cm)")

###############photon per cm single ##############
#plt.hist(phoPerLenList,bins=50)
#plt.title("Photon count per cm single cell")
#plt.xlabel("Photon count")

###############photon per cm Event ##############
#plt.hist(phoPerLenListEvent,bins=50)
#plt.title("Photon count per cm Event")
#plt.xlabel("Photon count")

# plt.show()

