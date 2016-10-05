import numpy as np
import pylab as P
import ROOT
from ROOT import gROOT
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-build/lib/libEventLib.so")
gROOT.ProcessLine(".L /home/mage/PROSPECT/PROSPECT-G4-Sec/include/Output/Event.hh")
Nx=14
Ny=10
EList=[]
XList=[]
YList=[]
ZList=[]
PXList=[]
PYList=[]
PZList=[]
FileEventDict={}
SingleColumnEventOut=open("SingleColumnEvent.out","w+")
SingleColumnEventOut.write("file, event\n")
outputTemp="%(file)i, %(event)i\n"
for k in xrange(0,16):
    print "processing file "+str(k)
    cryFile=ROOT.TFile("../cry_"+str(k)+".root")
    tree=ROOT.TTree()
    tree=cryFile.Get("PG4")
    sp=ROOT.SecondaryParticleEvent()
    ion=ROOT.IoniClusterEvent()
    prim=ROOT.ParticleEvent()
    tree.SetBranchAddress("SecParticle",sp)
    tree.SetBranchAddress("ScIoni",ion)
    tree.SetBranchAddress("Prim",prim)
    entries=int(tree.GetEntries())
    for i in xrange(0,entries):
	# print "event "+str(i)+" for file "+str(k)
	tree.GetEntry(i)
	clust=ion.nIoniClusts
	det=sp.nParticlesDet
	part=prim.nParticles
	clus=ion.nIoniClusts
	ylist=[]
	xlist=[]
	for j in xrange(0,clus):
	    vert=ion.clusts.At(j)
	    x=vert.x[0]
	    vol=vert.vol
	    pid=vert.PID
	    if vol>=0 and pid==13:
		coord=vert.x
		x=vol%Nx
		y=vol/(Nx)
		vol2=x+(Nx)*y
		if (y not in ylist):
		    ylist.append(y)
		if (x not in xlist):
		    xlist.append(x)
	if len(ylist) ==10 and len(xlist)==1:
	    FileEventDict["file"]=k
	    FileEventDict["event"]=i
	    output= outputTemp %FileEventDict
	    SingleColumnEventOut.write(output)
	    for j in xrange(0,part):
		vert=prim.particles.At(j)
		pid=vert.PID
		if pid==13:
		    XList.append(vert.x[0])
		    YList.append(vert.x[1])
		    ZList.append(vert.x[2])
		    PXList.append(vert.p[0])
		    PYList.append(vert.p[1])
		    PZList.append(vert.p[2])
		    EList.append(vert.E)
    cryFile.Close()
outTree= ROOT.TTree("Single_Column_Events", "single column events")
x=np.zeros(1,dtype=float)
y=np.zeros(1,dtype=float)
z=np.zeros(1,dtype=float)
px=np.zeros(1,dtype=float)
py=np.zeros(1,dtype=float)
pz=np.zeros(1,dtype=float)
E=np.zeros(1,dtype=float)
outTree.Branch("x",x,"x/d")
outTree.Branch("y",y,"y/d")
outTree.Branch("z",z,"z/d")
outTree.Branch("px",px,"px/d")
outTree.Branch("py",py,"py/d")
outTree.Branch("pz",pz,"pz/d")
outTree.Branch("E",E,"x/d")
SCE=len(XList)
preFile=open("premac.mac")
pre=preFile.read()
tempFile=open("template.mac")
temp=tempFile.read()
macDict={}
total=[]
for i in xrange(0,24):
    total.append(pre %(i,i))
for i in xrange(0,SCE):
# for i in xrange(0,24):
    x[0]=XList[i]
    y[0]=YList[i]
    z[0]=ZList[i]
    px[0]=PXList[i]
    py[0]=PYList[i]
    pz[0]=PZList[i]
    E[0]=EList[i]
    outTree.Fill()
    macDict["x"]=XList[i]
    macDict["y"]=YList[i]
    macDict["z"]=ZList[i]
    macDict["px"]=PXList[i]
    macDict["py"]=PYList[i]
    macDict["pz"]=PZList[i]
    macDict["E"]=EList[i]
    total[i%24]+=temp %macDict
for i in xrange(0,24):
    macFile=open("./cry_SCE"+str(i)+".mac",'w+')
    macFile.write(total[i])
    macFile.close()
out = ROOT.TFile("SingleColEvtPrim.root", "recreate")
outTree.Write()
out.Close()
SingleColumnEventOut.close()
