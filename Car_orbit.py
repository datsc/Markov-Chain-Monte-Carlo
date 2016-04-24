#######################################################################
# Calculating the final position and velocity in the Dehnen potential #
#######################################################################
# Takes the input from the cm code and outputs to a final file #
###############################################################
###############################################################
# Run with python Ugur_orbit1.py.py carinainital.dat carinafinal.dat #
###############################################################
#! /usr/bin/env python

from Orbit_Int import Deh_model,kpc_Myr2km_s,dt_frac
import sys

if len(sys.argv) < 3 :
    print >> sys.stderr, "Usage:", sys.argv[0], "Inputfile Outputfile "
    exit()

try :
    filename=sys.argv[1]
except :
    print filename, "could not be opened!"
    exit()


output=sys.argv[2]
initf=open(filename,'r')
line=initf.readline()
d=line.split()
x=float(d[0])
y=float(d[1])
z=float(d[2])
u=float(d[3])
v=float(d[4])
w=float(d[5])
dt_frac=float(d[6])#fraction of orbital time (t_orbit*dt_frac: I think this is the time step)
time=float(d[7])
initf.close()

print 'if you get', 'GalPot ERROR: Trying to construct Disks from a closed std::istream',' trying to construct it elsewhere change the pot4a address in the main python code: here!'
MM = Deh_model('/home/uural/hannilib/GalPot/pot.4a')


ans = MM.calc_orbit2(x,y,z,u,v,w,dt_frac,time)
#this is the result of calc_orbit2 in the Orbit_Int.py

vxugur=(-1)*ans[4]*kpc_Myr2km_s 
vyugur=(-1)*ans[5]*kpc_Myr2km_s 
vzugur=(-1)*ans[6]*kpc_Myr2km_s 
outorb=open(output,'w')
print >> outorb,ans[1],ans[2],ans[3],vxugur,vyugur,vzugur
outorb.close 


