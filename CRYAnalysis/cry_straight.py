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
for k in xrange(0,4):
    # if k!=1 and k!=9:
    cryFile=ROOT.TFile("../cry_"+str(k)+".root")
    tree=ROOT.TTree()
    tree=cryFile.Get("PG4")
    ion=ROOT.IoniClusterEvent()
    tree.SetBranchAddress("ScIoni",ion)
    entries=int(tree.GetEntries())
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
		histCell.Fill(x,y)
		if (y not in ylist):
		    ylist.append(y)
		if (x not in xlist):
		    xlist.append(x)
	if len(ylist) ==10 and len(xlist)==1:
	    count+=1
print count
print tot
print float(count)/float(tot)
histCell.Draw("colz")
raw_input()
