import sys
import math
import numpy 
import random
import scipy
from scipy import stats
from scipy import var
############################################

par1="alpharan.dat"
mu1= 0.0  
std1 = 1.0 
mrand = numpy.random.normal(mu1, std1,5000)  
outdata=open(par1,'w')
for i in xrange(5000):
    muran=mrand[i]
    print >>outdata,mrand[i]

par2="deltaran.dat"
mu2= 0.0  
std2 = 1.0 
mrand2 = numpy.random.normal(mu2, std2,5000)  
outdata2=open(par2,'w')
for j in xrange(5000):
    muran2=mrand2[j]
    print >>outdata2,mrand2[j]

par3="mdmran.dat"
mu3= 0.0  
std3 = 1.0 
mrand3 = numpy.random.normal(mu3, std3,5000)  
outdata3=open(par3,'w')
for k in xrange(5000):
    muran3=mrand3[k]
    print >>outdata3,mrand3[k]

par4="msran.dat"
mu4= 0.0  
std4 = 1.0 
mrand4 = numpy.random.normal(mu4, std4,5000)  
outdata4=open(par4,'w')
for l in xrange(5000):
    muran4=mrand4[l]
    print >>outdata4,mrand4[l]

par5="rdmran.dat"
mu5= 0.0  
std5 = 1.0 
mrand5 = numpy.random.normal(mu5, std5,5000)  
outdata5=open(par5,'w')
for m in xrange(5000):
    muran5=mrand5[m]
    print >>outdata5,mrand5[m]

par6="rstran.dat"
mu6= 0.0  
std6 = 1.0 
mrand6 = numpy.random.normal(mu6, std6,5000)  
outdata6=open(par6,'w')
for m in xrange(5000):
    muran6=mrand6[m]
    print >>outdata6,mrand6[m]
