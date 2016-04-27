# This code generates the initial random numbers from normal
# distribution as needed to simulate the random walk in the
# parameter space.
import sys
import math
import numpy 
import random
############################################

rmean=0
rstd=1
def genran(x):
    outdata=open(x,'w')
    varrand=numpy.random.normal(rmean, rstd,5000)
    for i in xrange(5000):
        print >> outdata,varrand[i]
    outdata.close()

# Writing out random number for 6 variables: alpha, delta proper motions,
# mass and radius of the dark matter halo, mass and radius of stellar component.

filenames=["alpharan.dat","deltaran.dat","mdmran.dat","msran.dat","rdmran.dat","rsran.dat"]

for new in filenames:
    genran(new)


