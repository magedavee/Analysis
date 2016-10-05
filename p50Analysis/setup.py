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
import time
start = time.time()
initROOT()
dir="../Data/p50/"
name=dir+"s004.root"
pmt=ROOT.PmtData(name)
pmt.Write("p50_s4_setup.root")
end = time.time()
print(end - start)
