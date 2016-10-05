from math import *
import numpy as np
import matplotlib.pyplot as plt
from random import *
t=[]
u=[]
for i in xrange(0,10**4):
    theta=uniform(-pi/2+.01,pi/2-.01)
    time=1/cos(theta)
    if time<100:
	u.append(theta)
	t.append(cos(time))
plt.hist(t,bins=100,range=[0,40])
print np.array(t).std()
meanT=np.array(t).mean()
Var=0
print len(t)
for i in t:
    Var+=((i-meanT)**2)/len(t)
print sqrt(Var)
plt.show()
