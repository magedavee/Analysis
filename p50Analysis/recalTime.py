from initROOT import initROOT
import ROOT
from ROOT import gROOT, TCanvas, TF1,TFile,TTree
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os.path
from scipy.stats import norm 
from array import array
import sys
import math
import time

##################Functions###########

##################Init###########
start = time.time()
initROOT()
dir="../Data/p50/"
name=dir+"s004.root"
# name=dir+"p50.root"
if os.path.isfile(name):
    print "looking at "+name
else:
    print "error "+name+" does not exist"
    quit()
setup="p50_s4_setup.root"
if os.path.isfile("./"+setup):
    pmt=ROOT.PmtData(name,setup)
else:
    print "no setup "+setup
    pmt=ROOT.PmtData(name)
pmt.CalTime()
pmt.Write("p50_s4_setup.root")
end = time.time()
print(end - start)
