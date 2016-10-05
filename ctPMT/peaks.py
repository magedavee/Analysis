import numpy as np
import matplotlib.pyplot as plt
import rootpy.ROOT as ROOT

class peaks:
    def __init__(self,filename):
	f=ROOT.TFile(filename)
	self.keys=f.GetListOfKeys()
	self.histDict={}
	self.yDict={}
	self.nameDict={}
	self.keyDict={}
	self.voltDict={}
	self.snDict={}
	self.fitDict={}
	self.fitListDict={}
	self.fitListDict["spe","mean"]=[]
	self.fitListDict["spe","sigma"]=[]
	self.fitListDict["mpe","mean"]=[]
	self.fitListDict["mpe","sigma"]=[]
	self.peakDict={}
	self.pict={}
	allName={}
	allHist={}
	allKeys={}
	revName={}
	snList=[]
	for key in self.keys:
	    name=key.GetName()
	    b=name[-41:]
	    volt=int(b[-10:-6])
	    sn=int(b[-14:-11])
	    hists=key.ReadObj()
	    allHist[sn,volt]=hists
	    allName[name]=[sn,volt]
	    allName[sn,volt]=name
	    allKeys[sn,volt]=key
	    allKeys[key]=[sn,volt]
	    if sn not in snList:
		snList.append(sn)
		voltList=np.array(volt,int)
		self.voltDict[sn]=voltList
	    else:
		a=np.append(self.voltDict[sn],volt)
		self.voltDict[sn]= a
	maxHist={}
	maxName={}
	for sn in snList:
	   voltList=self.voltDict[sn]
	   if sn ==134:
	       maxHist[sn]=allHist[sn,1600]
	   else:
	       if sn == 5:
		   maxHist[sn]=allHist[sn,1500]
	       else:
		   if sn==403:
		       maxHist[sn]=allHist[sn,1600]
		   else:
		       volt=voltList.max()
		       maxHist[sn]=allHist[sn,volt]
	   name=allName[sn,volt]
	   maxName[name]=[sn]
	   maxName[sn]=name
	    
	self.histDict['all']=allHist
	self.histDict['max']=maxHist
	self.nameDict['all']=allName
	self.nameDict['max']=maxName
	self.keyDict['all']=allKeys
	self.snDict['all']=snList
	speHist=[]
	speSN=[]
	mpeHist=[]
	mpeSN=[]
	backHist=[]
	backSN=[]
	self.histDict["spe"]=speHist
	self.snDict['spe']=speSN
	self.histDict["mpe"]=mpeHist
	self.snDict['mpe']=mpeSN
	self.histDict["back"]=speHist
	self.snDict["back"]=backSN
	noSpeHist=[]
	noSpeSN=[]
	badSN=[]
	badHist=[]
	self.snDict["bad mpe"]=badSN
	self.histDict["bad"]=badHist
	self.histDict["no spe"]=noSpeHist
	self.snDict["no spe"]=noSpeSN
	for sn in snList:
	    self.__reduceRes__(sn,1000)
	    self.__reduceRes__(sn,10000)
	for sn in snList:
	    self.__findPeaks(sn)
	    self.__fit__(sn)
	#for sn in self.snDict["spe"]:
	 #   print "sn: "+str(sn)
	 #   print "spe "+str(self.fitDict["spe","mean",sn])
	 #   print "mpe "+str(self.fitDict["mpe","mean",sn])

#########Private increase bin size by 20###########
    def __reduceRes__(self,sn,size=1000):
	hist=self.histDict['max'][sn]
        y=[]
        x=[]
        for i in xrange(0,size,20):
            tot=0
            for j in xrange(i,i+20):
                tot=tot+hist.GetBinContent(j)
            y.append(tot/20)
            x.append(i)
	    tag="array"+str(size)
	yRay=np.array(y)
        self.yDict[tag,sn]=yRay
##########good for making plots##########################	
#	plt.plot(x,y)
#	max=self.voltDict[sn].max()
#	name=str(sn)+"_@_"+str(max)+"_"+str(size)
#	plt.title(name)
#	plt.xlabel("charge recoreded")
#	plt.savefig("./figs/"+str(size)+"/"+name, bbox_inches='tight')
#	plt.clf()
#	plt.show()

#########Private function to find Peaks###########
    def __findPeaks(self,sn,):
	y=self.yDict["array1000",sn]
	aMax=y.argmax()
	aMin=y.argmin()
	pList=[]
	pic=[]
	if aMax<aMin:
	    start=y.argmax()
	    right=y[y.argmax():y.argmin()]
	    left=y[0:y.argmax()+1]

	    y2=left
	    while  len(y2)!=0:
		aMax=y2.argmax()
		r=aMax
		while r==y2.argmax():
		    r-=1
		    y2=left[0:r+1]
		    if len(y2)==0:
			break
		pic.append(r+1)
		pList.append(aMax)
	    self.pict[sn]=pic

	    y2=right
	    r=0
	    thresh=y.max()*.1
	else:
	    if sn == 277:
		pList.append(4)
		pList.append(13)
		pic.append(5)
		pic.append(0)
		self.pict[sn]=pic

	self.peakDict[sn]=pList

############fit gausian#############################################
    def __fit__(self,sn):
	pList=self.peakDict[sn]
	hist=self.histDict['max'][sn]
	if sn==186 or sn==430: 
	    print str(sn)+" skip"
	else:
	    if len(pList)>1:
		pic= self.pict[sn]
		start=0;
		mid=(pic[0]-1)*20
		if sn==277:
		    mid=100
		    end=400
		else:
		    end=self.yDict["array1000",sn].argmin()*20
		try:
		    res=hist.Fit("gaus","S","",mid,end)
		    meanSpe=hist.GetFunction("gaus").GetParameter(1)
		    sigmaSpe=hist.GetFunction("gaus").GetParameter(2)
		    self.fitDict["spe","mean",sn]=meanSpe
		    self.fitDict["spe","sigma",sn]=sigmaSpe
		    self.fitListDict["spe","mean"].append(meanSpe)
		    self.fitListDict["spe","sigma"].append(sigmaSpe)
		    self.histDict["spe"].append(hist)
		    self.snDict['spe'].append(sn)
		    print sn
		    if int(res)!=0:
			print str(sn)+" was not fitted"
		except TypeError:
			print str(sn)+" was not fitted"
		try:
		    res=hist.Fit("gaus","S","",start,mid)
		    if int(res)==0:
			meanBack=hist.GetFunction("gaus").GetParameter(1)
			sigmaBack=hist.GetFunction("gaus").GetParameter(2)
			self.fitDict["back","mean",sn]=meanBack
			self.fitDict["back","sigma",sn]=sigmaBack
			self.histDict["back"].append(hist)
			self.snDict['back'].append(sn)
		    else:
			print str(sn)+" background was not fitted"
		except TypeError:
			print str(sn)+" was not fitted"
	    else:
		if sn ==211:
		    start=200
		    end=400
		else:
		    if sn ==236:
			start=100
			end=300
		    else:
			start=0
			end=self.yDict["array1000",sn].argmin()*20
	    try:
		res=hist.Fit("gaus","S","",start,end)
		meanBack=hist.GetFunction("gaus").GetParameter(1)
		sigmaBack=hist.GetFunction("gaus").GetParameter(2)
		self.fitDict["back","mean",sn]=meanBack
		self.fitDict["back","sigma",sn]=sigmaBack
		self.histDict["no spe"].append(hist)
		self.snDict["no spe"].append(sn)
		self.histDict["back"].append(hist)
		self.snDict['back'].append(sn)
		if not(sn in self.snDict['spe']):
		    if meanBack>130 and sn != 157 and sn != 42:
			self.fitDict["spe","mean",sn]=meanBack
			self.fitDict["spe","sigma",sn]=sigmaBack
			self.fitListDict["spe","sigma"].append(sigmaBack)
			self.fitListDict["spe","mean"].append(meanBack)
			self.histDict["spe"].append(hist)
			self.snDict["spe"].append(sn)
		if int(res)!=0:
		    print str(sn)+" background was not fitted"
	    except TypeError:
		    print str(sn)+" background  was not fitted"

	 
	    start=self.yDict["array1000",sn].argmin()*20
	    end=9999
	    try:
		res=hist.Fit("gaus","S","",start,end)
		meanMpe=hist.GetFunction("gaus").GetParameter(1)
		sigmaMpe=hist.GetFunction("gaus").GetParameter(2)
		self.fitDict["mpe","mean",sn]=meanMpe
		self.fitDict["mpe","sigma",sn]=sigmaMpe
		self.fitListDict["mpe","sigma"].append(sigmaMpe)
		self.fitListDict["mpe","sigma"].append(meanMpe)
		if sn == 277:
		    res=hist.Fit("gaus","S","",2000,8000)
		else:
		    while int(res)!=0 and start< end:
			start+=10
			res=hist.Fit("gaus","S","",start,end)
		print sn
		meanMpe=hist.GetFunction("gaus").GetParameter(1)
		sigmaMpe=hist.GetFunction("gaus").GetParameter(2)
		self.fitDict["mpe","mean",sn]=meanMpe
		self.fitDict["mpe","sigma",sn]=sigmaMpe
		self.histDict["mpe"].append(hist)
		self.snDict['mpe'].append(sn)
		if start >= end:
		    print str(sn)+" mpe was not fitted"
		    self.snDict["bad mpe"].append(sn)
		    self.histDict["bad"].append(hist)
	    except TypeError:
		    print str(sn)+" mpe  was not fitted"


	
###############Plot fit Parameter##################
    def plotParameter(self,peak,parameter,toFile=False):
	paraList=[]
	num=0
	if parameter == "ratio":
	    for sn in self.snDict[peak]:
		mean=self.fitDict[peak,"mean",sn]
		sigma=self.fitDict[peak,"sigma",sn]
		rat=mean/sigma
		paraList.append(rat)
		plt.title("ratio mean/sigma "+str(peak))
		num+=1
		name="ratio"+str(peak)
	else:
	    if parameter == "mpe/spe":
		print "mpe/spe"
		for sn in self.snDict["spe"]:
		    if sn != 186: 
			meanSpe=self.fitDict["spe","mean",sn]
			meanMpe=self.fitDict["mpe","mean",sn]
			rat=meanMpe/meanSpe
			#print str(sn)+": "+str(rat)
			#print meanMpe
			#print meanSpe
			paraList.append(rat)
			plt.title("ratio MPE/SPE ")
			num+=1
			name="ratio_"+str("MPE_SPE")

	    else:
		if parameter == "sig2":
		    temp=[]
		    for sn in self.snDict["spe"]:
			if sn != 186: 
			    mean=self.fitDict["spe","mean",sn]
			    temp.append(mean)
		    meanSpe=np.array(temp).mean()
		    for sn in self.snDict["mpe"]:
			if sn != 186: 
			    sigMpe=self.fitDict["mpe","sigma",sn]
			    pE=(sigMpe*sigMpe)/(meanSpe*meanSpe)
			    if pE>=15:
				print str(sn)+": "+str(pE)
				paraList.append(pE)
				plt.title("ratio MPE Sig^2/Spe Sig^(2) ")
				num+=1
			    name="ratio_"+str("MPESig_SPESig")
		else:
		    if parameter =='pe':
			if peak =='best':
			    plt.title("Photoelectrons")
			    name="PhotoElect"
			    for sn in self.snDict["all"]:
				if sn==430 or sn==186:
				    continue
				if sn in self.snDict["spe"]:
				    PE=self.fitDict["mpe","mean",sn]/self.fitDict["spe","mean",sn]
				else:
				    aveSpe=np.array(self.fitListDict["spe","mean"]).mean()
				    aveSigma=np.array(self.fitListDict["spe","sigma"]).mean()
				    mpeSigma=self.fitDict["mpe","sigma",sn]
				    PE=mpeSigma**2/(aveSpe**2+aveSigma**2)
				paraList.append(PE)
				num+=1
			
		    else:
			for sn in self.snDict[peak]:
			    paraList.append(self.fitDict[peak,parameter,sn])
			    plt.title(str(peak)+" "+str(parameter))
			    num+=1
			    name=str(peak)+"_"+str(parameter)
	xy=(10,-15)
	entry= num
	para=np.array(paraList)
	mean= para.mean()
	std= para.std()
	lable="mean "+str(mean)
	plt.annotate(lable, xy=(0, 1),  xycoords="figure points",xytext=(0, 25), textcoords='figure points')
	lable="sigma "+str(std)
	plt.annotate(lable, xy=(0, 1),  xycoords="figure points",xytext=(200, 25), textcoords='figure points')
	lable="entry "+str(num)
	plt.annotate(lable, xy=(0, 1),  xycoords="figure points",xytext=(400, 25), textcoords='figure points')
	plt.hist(paraList)
	if toFile:
	    plt.savefig("./figs/"+name, bbox_inches="tight")
	plt.show()
	plt.clf()
	print str(peak)+" "+str(parameter)
	print "mean: "+str(mean)
	print "std: "+str(std)

	print "num :"+str(num)
	i=0
	return [mean,std,num]
	#for sn in self.snDict[peak]:
	    #print str(sn)+": " +str(paraList[i])
	    #i+=1
###############Plot Hist using Root################
    def plotHist(self,cat,sn):
	self.histDict[cat][sn].Draw()
####################End of File####################
