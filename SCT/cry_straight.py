import numpy as np
import pylab as P
import ROOT
from ROOT import gROOT
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-build/lib/libEventLib.so")
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-Sec/include/Output/Event.hh")
histCell=ROOT.TH2D("Cell Ionization Hits","Cell Ionization Hits",14,0,14,10,0,10)
pmt1Hist=ROOT.TH2D("d","d",100,-1200,1200,100,-1200,1200)
pmt2Hist=ROOT.TH2D("d","d",100,-1200,1200,100,-1200,1200)
cellHitX=ROOT.TH1D("cellH","Cell's Hit",15,0,15)
cellHitY=ROOT.TH1D("cellH","Cell's Hit",11,0,11)
Nx=14
Ny=10
goodEvent=[]
count=0
tot=0
for k in xrange(0,1):
    # if k!=1 and k!=9:
    cryFile=ROOT.TFile("../cry_SCE"+str(k)+".root")
    tree=ROOT.TTree()
    tree=cryFile.Get("PG4")
    ion=ROOT.IoniClusterEvent()
    sp=ROOT.SecondaryParticleEvent()
    tree.SetBranchAddress("ScIoni",ion)
    tree.SetBranchAddress("SecParticle",sp)
    entries=int(tree.GetEntries())
    entries=10
    for i in xrange(0,entries):
	tree.GetEntry(i)
	clust=ion.nIoniClusts
	clus=ion.nIoniClusts
	E=ion.EIoni
	e=0
	ylist=[]
	xlist=[]
	t=0;
	tot+=1
	if i ==0:
	    t=ion.clusts.At(0).t
	for j in xrange(0,clus):
	    vert=ion.clusts.At(j)
	    time=vert.t-t
	    vol=vert.vol
	    pid=vert.PID
	    if vol>=0 and pid==13:
		coord=vert.x
		x=vol%Nx
		y=vol/(Nx)
		vol2=x+(Nx)*y
		# print "x: "+str(x)+" x coord " +str(coord[0])+" y: "+str(y)+"cm  z coord "+str(coord[2])+"cm"+" vol :"+str(vol)
		# print "x: "+str(x)+ " y: "+str(y)+" vol :"+str(vol)+" vol2: "+str(vol2)
		# print "x coord "+str(coord[0])+"cm  y coord "+str(coord[1])+"cm  z coord "+str(coord[2])+"cm"+" time "+str(time)+"ns"
		# print ""
		if (y not in ylist):
		    ylist.append(y)
		if (x not in xlist):
		    xlist.append(x)
	if len(ylist) ==10 and len(xlist)==1:
	    count+=1
	    for j in xrange(0,10):
		histCell.Fill(xlist[0],ylist[j])
	    # histCell.Fill(x,y)
	    # for j in xrange(0,clus):
		# coord=vert.x
		# x=vol%Nx
		# y=vol/(Nx)
		# vol=vert.vol
		# pid=vert.PID
		# # if vol>-1 and pid ==13:
		# print "x :"+str(x)
		# print "y :"+str(y)
		# histCell.Fill(x,y)
	    # # print "event "+str(i)+" for file "+str(k)
	    # print "xlist: "+str(xlist)
	    # print "ylist: "+str(ylist)
	cellHitY.Fill(len(ylist))
	cellHitX.Fill(len(xlist))
print count
print tot
histCell.Draw("colz")
raw_input()
