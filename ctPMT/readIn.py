#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import rootpy.ROOT as ROOT


class readIn:
    ############Constructor#########################
    def __init__(self, filename):
        self.pdfs = {}
        self.voltDict={}
        self.peakDict={}
        self.nameList=[]
        self.yDict={}
        self.yDictTot={}
        self.maxIvals={}
        self.minIvals={}
        f = ROOT.TFile(filename)
        self.keyDict={}
        self.rejected=[]
        self.rejected2=[]
        self.rejectedMpe=[]
        self.mean=[]
        self.sigma=[]
        self.mean2=[]
        self.sigma2=[]
        self.meanMpe=[]
        self.sigmaMpe=[]
        self.ratio=[]
        self.ratioSpe=[]
        self.ratioBack=[]
        #hist keys
        self.keys = f.GetListOfKeys()
        #list that contains hist
        self.histList=[]
        for key in self.keys:
            hists = key.ReadObj()
            self.keyDict[key.GetName()]=key
            self.histList.append(hists)
            ######list of key name
            self.nameList.append(key.GetName())
            self.pdfs[key.GetName()]=hists
        ######Same as nameList but only for hist over 1700
        self.label=[]
        ############serial number
        self.sn=[]
        self.snDict={}
        for i in range(0,len(self.nameList)):
            b=self.nameList[i][-41:]
            c= b[-10:-6]
            d=b[-14:-11]
            volt=int(c)
            self.voltDict[self.nameList[i]]=volt
            if volt >= 1700:
                self.sn.append(int(d))
                self.snDict[self.nameList[i]]=int(d)
                self.label.append(self.nameList[i])

#########Get a histogram of spe peaks#############
    def getSpes(self):
        print "Getting SPE Peaks"

        for lab in self.label:
            self.__reduceRes(lab)
        for lab in self.label:
            print "finding peak "+str(lab)
            self.__spilt(lab)
            self.__fitSpePeek(lab)
        for name in self.rejected:
            print name+" not fitted"
        for name in self.rejected2:
            print name+" not fitted again"
        for lab in self.label:
            self.__findMpe(lab)
        for name in self.rejectedMpe:
            print name+" not fitted"
        plt.hist(self.sigmaMpe,50)
        for i in xrange(0,len(self.mean2)):
            m=self.mean2[i]
            s=self.sigma2[i]
            rat=m/s
            self.ratioSpe.append(rat)
        for i in xrange(0,len(self.mean)):
            m=self.mean[i]
            s=self.sigma[i]
            rat=m/s
            self.ratioBack.append(rat)
	xy=(10,-5), 
	xycoords='axes points'
	ave=np.array(self.sigmaMpe).mean()
	print "sigma mpe" + str(ave)
        plt.title("sigma MPE")
        plt.show()
	ave=np.array(self.meanMpe).mean()
	print "mean Mpe"+str( ave)
        plt.hist(self.meanMpe,50)
        plt.title("mean MPE")
        plt.show()
	ave=np.array(self.ratio).mean()
	print "ratio2 mean sigma"+str(ave)
        plt.hist(self.ratio,50)
        plt.title("ratio mean/sigma MPE")
        plt.show()
        plt.title("mean SPE")
	ave=np.array(self.mean2).mean()
	print "mean SPE"+str( ave)
        plt.hist(self.mean2,50)
        plt.show()
        plt.title("sigma SPE")
	ave=np.array(self.sigma2).mean()
	print"sigma SPE" + str(ave)
        plt.hist(self.sigma2,50)
        plt.show()
        plt.title("mean ratio mean/sigma SPE")
	ave=np.array(self.ratioSpe).mean()
	print "ratio"+str( ave)
        plt.hist(self.ratioSpe,50)
        plt.show()
        plt.title("mean background")
	ave=np.array(self.mean).mean()
	print "mean background"+str(ave)
        plt.hist(self.mean,50)
        plt.show()
        plt.title("sigma background")
	ave=np.array(self.sigma).mean()
	print "sigma background"+str(ave)
        plt.hist(self.sigma,50)
        plt.show()
        plt.title("mean ratio mean/sigma background")
	ave=np.array(self.ratioBack).mean()
	print "ratio2"+str( ave)
        plt.hist(self.ratioBack,50)
        plt.show()

#########Find MPE peak#############################
    def __findMpe(self,name):
        hist=self.pdfs[name]
        hist.Fit("gaus","","",1000,10000)
        mean=hist.GetFunction("gaus").GetParameter(1)
        sigma=hist.GetFunction("gaus").GetParameter(2)
        if mean <0:
            hist.Fit("gaus","","",500,2000)
            mean=hist.GetFunction("gaus").GetParameter(1)
            sigma=hist.GetFunction("gaus").GetParameter(2)
        self.meanMpe.append(mean)
        self.sigmaMpe.append(sigma)
	rat=mean/sigma
        self.ratio.append(rat)
        return


#########Private increase bin size by 20###########
    def __reduceRes(self,name):
        hist=self.pdfs[name]
        y=[]
        x=[]
        for i in xrange(0,1000,20):
            tot=0
            for j in xrange(i,i+20):
                tot=tot+hist.GetBinContent(j)
            y.append(tot/20)
            x.append(i)
        self.yDict[name]=y

#########Private increase bin size by 20###########
    def __reduceResTot(self,name):
        hist=self.pdfs[name]
        y=[]
        x=[]
        for i in xrange(0,10000,20):
            tot=0
            for j in xrange(i,i+20):
                tot=tot+hist.GetBinContent(j)
            y.append(tot/20)
            x.append(i)
        self.yDictTot[name]=y

#########Private this was spilt because matplot bug###########
    def __spilt(self,name):
        y=self.yDict[name]
        yMins=[]
        yMaxs=[]
        self.__findPeaks(y,'max',1,yMins,yMaxs)
        self.maxIvals[name]=yMaxs
        self.minIvals[name]=yMins

#########Private function to find Peaks###########
    def __findPeaks(self,y,minmax,start,yMins,yMaxs):
        while minmax =='min' or minmax == 'max':
            if minmax is 'min':
                res=self.__findMin(y,start)
                if res !=0:
                    yMins.append(res)
                    start=res[-1]
                    minmax='max'
                else:
                    minmax='0'
            else:
                res=self.__findMax(y,start)
                if res !=0:
                    yMaxs.append(res)
                    start=res[-1]
                    minmax='min'
                else:
                    minmax='0'


#########Find max ##############################################
    def __findMax(self,y,start):
        if start >= len(y)-1:
            print "start hit end"
            return 0
        yStart=y[start]
        i=start+1
        end=i

        while y[i]>yStart:
            i+=1
            end=i
            if i>=len(y):
                break
        a = np.array(y[start:end], int)
        res=[a.max(),a.argmax(),start,end]
        start=end
        return res

#########Find min ##############################################
    def __findMin(self,y,start):
        if start >= len(y)-1:
            print "start hit end"
            return 0
        yStart=y[start]
        min=yStart
        index=start
        count =0
        for j in xrange(start+1,len(y)):
            diff=min-y[j]
            if diff > -1:
                min=y[j]
                index=j
                count =0
            else:
                if index !=start:
                    count+=1
            if count==5:
                break
        if index !=start:
            end=index
        else:
            end=50
        a = np.array(y[start:end], int)
        res=[a.min(),a.argmin(),start,end]
        start=end
        return res
#########fitting gaussian ######################################
    def __fitSpePeek(self,name):
        iVal=self.maxIvals[name][0]
        start=iVal[2]*20
        end=iVal[3]*20
        hist=self.pdfs[name]
        try:
            hist.Fit("gaus","","",start,end)
            self.mean.append(hist.GetFunction("gaus").GetParameter(1))
            self.sigma.append(hist.GetFunction("gaus").GetParameter(2))
        except TypeError:
                print name+" was not fitted"
                self.rejected.append(name)
        if len(self.minIvals[name]) ==2:
            iVal=self.maxIvals[name][1]
            start=iVal[2]*20
            end=iVal[3]*20
            try:
                hist.Fit("gaus","","",start,end)
                self.mean2.append(hist.GetFunction("gaus").GetParameter(1))
                self.sigma2.append(hist.GetFunction("gaus").GetParameter(2))
            except TypeError:
                    print name+" was not fitted"
                    self.rejected2.append(name)

        return

#########Juest Sanity check delete after########################
    def __showYDict(self,name):
            y=self.yDictTot[name]
            plt.title(name)
            plt.plot(y)
            plt.show()

#########Juest Sanity check delete after########################
    def __showYDictTot(self,name):
            y=self.yDictTot[name]
            plt.title(name)
            plt.plot(y)
            plt.show()

#########Get the numth hist########################
    def getHist(self,num):
        if type(num) is not int:
            print "second arguement needs to be type int"
            return
        if num>=len(self.pdfs) or num <0:
            print str(num)+" is outside range "+str(len(self.pdfs)-1)
            return
        return self.histList[num]

###########Get the number of Hists##################
    def getNumOfHist(self):
        return len(self.histList)

###########Get the number of Hists##################
    def drawHist(self,num):
        self.histList[num].Draw()
        input("")

##################################################
