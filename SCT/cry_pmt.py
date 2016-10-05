import numpy as np
import pylab as P
import ROOT
from ROOT import gROOT
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-build/lib/libEventLib.so")
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-Sec/include/Output/Event.hh")
# histCell=ROOT.TH2D("Cell Ionization Hits","Cell Ionization Hits",10,0,10,14,0,14)
# pmt1Hist=ROOT.TH2D("d","d",100,-1200,1200,100,-1200,1200)
Nx=14
Ny=10
print "start"
for k in xrange(0,32):
    # if k!=1 and k!=9:
    cryFile=ROOT.TFile("../cry_"+str(k)+".root")
    tree=ROOT.TTree()
    tree=cryFile.Get("PG4")
    sp=ROOT.SecondaryParticleEvent()
    ion=ROOT.IoniClusterEvent()
    tree.SetBranchAddress("SecParticle",sp)
    tree.SetBranchAddress("ScIoni",ion)
    entries=int(tree.GetEntries())
    for i in xrange(0,entries):
	print "event "+str(i)+" for file "+str(k)
	tree.GetEntry(i)
	clust=ion.nIoniClusts
	det=sp.nParticlesDet
	clus=ion.nIoniClusts
	E=ion.EIoni
	e=0
	for j in xrange(0,clus):
	    vert=ion.clusts.At(j)
	    x=vert.x[0]
	    vol=vert.vol
	    x=vol/Nx
	    y=vol%(Ny)
	    # histCell.Fill(x,y)
	    e+=vert.E

	percent=-1
	for j in xrange(0,det):
	    percentUp=int(float(j)/float(det)*100)
	    if percentUp != percent:
		percent=percentUp
		if percent%10==0:
		    print str(percentUp)+" %"
	    photo=sp.particles.At(j)
	    x=photo.x[0]
	    y=photo.x[1]
	    z=photo.x[2]
	    #print "x: "+str(x)+" y: "+str(y)+" z: "+str(z)
	    # if y<0:
		# pmt1Hist.Fill(x,z)
pmt1Hist.Draw("colz")
