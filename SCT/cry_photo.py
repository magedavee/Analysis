from initROOT import initROOT
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import gROOT
from ROOT import * 
from math import *
#this script sums each vetex to the cell.
# initROOT()
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-build/lib/libEventLib.so")
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-Sec/include/Output/Event.hh")
histCell=ROOT.TH2D("Cell Ionization Hits","Cell Ionization Hits",14,0,14,10,0,10)
histCell0=ROOT.TH2D("Cell Ionization Hits","Cell Ionization Hits",14,0,14,10,0,10)
pmtCell=ROOT.TH2D("Cell Ionization Hits","Cell Ionization Hits",14,0,14,10,0,10)
pmt1Cell=ROOT.TH2D("pmt Hits","Cell Ionization Hits",14,0,14,10,0,10)
pmt1Hist=ROOT.TH2D("pmtpos","pmt hits",100,-1200,1200,100,-1200,1200)
Nx=14
Ny=10
lat=144.198
goodEvent=[]
count=0
tot=0
print "start photo"
dtList=[]
lvarList=[]
rvarList=[]
vartList=[]
ratList=[]
countList=[]
rightList=[]
leftList=[]
lenList=[]
lengList=[]
posxList=[]
posyList=[]
poszList=[]
energyList=[]
spreadList=[]
volList=[]
nxList=[]
nyList=[]
evtList=[]
num=0
k=22
p50=[[],[]]
for K in xrange(1,2):
    # if k!=1 and k!=9:
    cryFile=ROOT.TFile("../cry_SCE"+str(k)+".root")
    tree=ROOT.TTree()
    print type(tree)
    tree=cryFile.Get("PG4")
    print type(tree)
    ion=ROOT.IoniClusterEvent()
    sp=ROOT.SecondaryParticleEvent()
    print type(tree)
    tree.SetBranchAddress("ScIoni",ion)
    tree.SetBranchAddress("SecParticle",sp)
    entries=int(tree.GetEntries())
    print "entries: "+str(entries)
    # entries=10
    percent=-1
    for i in xrange(1,entries):
	pmtL=[]	
	pmtR=[]	
	cells=[]	
	for j in xrange(0,Nx*Ny):
	    pmtL.append([])
	    pmtR.append([])
	    cells.append([])
	percentUp=int(float(i)/float(entries)*100)
	if percentUp != percent:
	    percent=percentUp
	    # if percent%10==0:
	    print "\n file "+str(k)+" "+str(percentUp)+" %"
	# print "event "+str(i)+" out of  "+str(entries)+" for file "+str(k)
	tree.GetEntry(i)
	num=i
	clus=ion.nIoniClusts
	det=sp.nParticlesDet
	t=0;
##start time
	# if i ==0:
	t=ion.clusts.At(0).t
##particle level
	totlen=0
	xCell=[]
	yCell=[]
	for j in xrange(0,clus):
	    vert=ion.clusts.At(j)
	    x=vert.x[0]
	    y=vert.x[2]
	    z=vert.x[1]
	    time=vert.t-t
	    e=vert.E
	    pid=vert.PID
	    vol=vert.vol
	    length=vert.l
	    # print clus
	    # print pid
	    # print vol
	    if vol>=0 and (pid==13 or pid==11):
		nx=vol%Nx
		ny=vol/(Nx)
		tup=[x,y,z,t,e,length]

		# print "vol: "+str(vol)
		cells[vol].append(tup)
                if nx%2==0 :
                    p50[0].append(tup)    
                else :
                    p50[1].append(tup)    

		# print "ion x: "+str(x)+" nx: "+str(nx)
		# print "ion z: "+str(z)+" ny: "+str(ny)
		# print "vol: "+str(vol)
		if (ny not in yCell):
		    yCell.append(ny)
		if (nx not in xCell):
		    xCell.append(nx)
		histCell0.Fill(nx,ny)
	
	if len(yCell) ==10 and len(xCell)==1:
	    # print "event: "+str(i)
	    for j in xrange(0,det):
		photo=sp.particles.At(j)
		x=photo.x[0]
		y=photo.x[1]
		z=photo.x[2]+450
		pid=photo.PID
		if pid == 0:
		    if y<0:
			pmt1Hist.Fill(x,z)
			nx=floor(x/lat)+Nx/2.0
			ny=floor(z/lat)+Ny/2.0
			pmtNum=int(nx+Nx*ny)
			pmt1Cell.Fill(nx,ny)
			pmtL[pmtNum].append(photo.t-t)
			# print "pmt x: "+str(x)+" nx: "+str(nx)
			# print "pmt z: "+str(z)+" ny: "+str(ny)
			# print "pmt num: "+str(pmtNum)
		    else:
			# pmt1Hist.Fill(x,z)
			pmt1Cell.Fill(nx,ny)
			nx=floor(x/lat)+Nx/2.0
			ny=floor(z/lat)+Ny/2.0
			pmtNum=int(nx+Nx*ny)
			pmtR[pmtNum].append(photo.t-t)

	    for j in xrange(0,Nx*Ny):
		if len(cells[j])>0 and (len(pmtR[j])>500 and len(pmtL[j])>500):
		    left=len(pmtL[j])
		    right=len(pmtR[j])
		    total=left+right
		    rat=float(left)/float(total)
		    xList=[]
		    yList=[]
		    zList=[]
		    z=0
		    tlList=[]
		    trList=[]
		    pos0=0
		    pos1=0
		    totlen=0
		    E=0
		    for h in xrange(0,len(cells[j])):
			if h==0:
			    pos0=[cells[j][h][0],cells[j][h][1],cells[j][h][2]]
			if h==(len(cells[j])-1):
			    pos1=[cells[j][h][0],cells[j][h][1],cells[j][h][2]]
			xList.append(cells[j][h][0])
			yList.append(cells[j][h][1])
			zList.append(cells[j][h][2])
			E+=cells[j][h][4]
			totlen+=cells[j][h][5]
		    l=sqrt((pos0[0]-pos1[0])**2+(pos0[1]-pos1[1])**2+(pos0[2]-pos1[2])**2)
		    # print "l: "+str(l)
		    # print "length: "+str(totlen)
		    for h in xrange(0,len(pmtR[j])):
			trList.append(pmtR[j][h])
		    for h in xrange(0,len(pmtL[j])):
			tlList.append(pmtL[j][h])
		    xPos=np.array(xList,'f')
		    yPos=np.array(yList,'f')
		    zPos=np.array(zList,'f')
		    tl=np.array(tlList,'f')
		    tr=np.array(trList,'f')
		    dt=tl.mean()-tr.mean()
		    lvarList.append(tl.std())
		    rvarList.append(tr.std())
		    dtList.append(dt)
		    posxList.append(xPos.mean())
		    posyList.append(yPos.mean())
		    poszList.append(zPos.mean())
		    energyList.append(E)
		    spreadList.append(yPos.std())
		    
		    nx=j%Nx
		    ny=j/(Nx)
		    # print "nx: "+str(nx)
		    # print "ny: "+str(ny)
		    pmtCell.Fill(nx,ny)
		    nxList.append(nx)
		    nyList.append(ny)
		    volList.append(j)
		    rightList.append(right)
		    leftList.append(left)
# pythgroan
		    lenList.append(l)
#geant
		    lengList.append(totlen)
		    evtList.append(num+k*400)



#Data output
outTree= ROOT.TTree("Muon Calibration", "Muon Calibration")
x=np.zeros(1,dtype=float)
y=np.zeros(1,dtype=float)
z=np.zeros(1,dtype=float)
#py
l=np.zeros(1,dtype=float)
#geant
leng=np.zeros(1,dtype=float)
lVar=np.zeros(1,dtype=float)
rVar=np.zeros(1,dtype=float)
yVar=np.zeros(1,dtype=float)
ratio=np.zeros(1,dtype=float)
dt=np.zeros(1,dtype=float)
E=np.zeros(1,dtype=float)
rightCount=np.zeros(1,dtype=int)
leftCount=np.zeros(1,dtype=int)
vol=np.zeros(1,dtype=int)
NX=np.zeros(1,dtype=int)
NY=np.zeros(1,dtype=int)
evt=np.zeros(1,dtype=int)
outTree.Branch("E",E,"E/D")
#py
outTree.Branch("l",l,"l/D")
#geant
outTree.Branch("length",leng,"leng/D")
outTree.Branch("lVar",lVar,"lVar/D")
outTree.Branch("rVar",rVar,"rVar/D")
outTree.Branch("x",x,"x/D")
outTree.Branch("y",y,"y/D")
outTree.Branch("z",z,"z/D")
outTree.Branch("yVar",yVar,"yVar/D")
outTree.Branch("dt",dt,"dt/D")
outTree.Branch("right",rightCount,"rightCount/I")
outTree.Branch("leftt",leftCount,"leftCount/I")
outTree.Branch("vol",vol,"vol/I")
outTree.Branch("nx",NX,"NX/I")
outTree.Branch("ny",NY,"NY/I")
outTree.Branch("evt",evt,"evt/I")

for i in xrange(0,len(posxList)):
    # print "pos: "+str(posxList[i])
    x[0]=posxList[i]
    y[0]=posyList[i]
    z[0]=poszList[i]
    # print "pos var: "+str(spreadList[i])
    yVar[0]=spreadList[i]
    # print "dt: "+str(dtList[i])
    dt[0]=dtList[i]
    vol[0]=volList[i]
    NX[0]=int(nxList[i])
    NY[0]=int(nyList[i])
    histCell.Fill(NX,NY)
    rightCount[0]=rightList[i]
    leftCount[0]=leftList[i]
    E[0]=energyList[i]
#py
    l[0]=lenList[i]
#geant
    leng[0]=lengList[i]
    lVar[0]=lvarList[i]
    rVar[0]=rvarList[i]
    evt[0]=evtList[i]
    outTree.Fill()
out = ROOT.TFile("Photon"+str(k)+".root", "recreate")
outTree.Write()
out.Close()
plt.hist(np.array(lenList), bins=400)
# plt.show()
print "done"
# cell hit
# histCell.Draw("colz")
# histCell0.Draw("colz")
#photon hit
pmt1Cell.Draw("colz")
# pmtCell.Draw("colz")
# pmt1Hist.Draw("colz")
raw_input()
