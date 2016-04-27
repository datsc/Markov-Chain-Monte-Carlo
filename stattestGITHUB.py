############################################################################
# CALCULATING THE CHI-SQUARES OF THE OBSERVABLES IN THE SIMULATIONS        # 
#
# Run with python stattestGITHUB.py gensimdispnoisy1.dat gensimbnoisy1.dat simdisp.dat simb.dat simchi.dat   #
#
# Modifiables:)!! 1)filename, 2)file format for the simulations,3) bin number   #      
#################################################################################

import sys
import math
from pytipsy import ifTipsy,TipsyHeader,TipsyDarkParticle,tipsypos,TipsyStarParticle
import scipy.stats 
import numpy
from numpy import arange

############################################
#  0  #    Takes two arguments and n bins  #
############################################
if len(sys.argv) <=6:
    print 'The chi2 code did not work since it did not receive the arguments'
#elif len(sys.argv)==3:
obssigma=sys.argv[1]
obssb=sys.argv[2]
simulationsigma=sys.argv[3]
simulationsb=sys.argv[4]
simulationgrad=sys.argv[5]
output=sys.argv[6]
print simulationsigma,simulationsb,output
binsigma=12
binsb=12 #11

###################################################
# 1 #       Here we read the observational data   #
###################################################
#############################
# 1.0 # Velocity gradient #
#############################
gro=[]
gro=[5.38]
dgro=[]
dgro=[1.58]#[0.145504436275] real sim value without noise so for this gro/2 is used which is almost the same with the Munoz2008 value of sqrt(2.5^2+2.5^2)
gradobs=numpy.array(gro) #10.0
dvgradobs=numpy.array(dgro)#5.0


#############################
# 1.1 # Velocity dispersion #
#############################

infile=open(obssigma,'r')
# this file now excludes the first line which was 0.0857038 5.8263606266e+00 5.9693324744e-01

sigma_obs_list=[]
dsigma_obs_list=[]
rkpc_list=[]

def file2listobs(filename,listarkpc_obs,listasigma_obs,listadsigma_obs):
    for line in filename:
        d1=line.split()
        listarkpc_obs+=[float(d1[0])]
        listasigma_obs+=[float(d1[1])]
        listadsigma_obs+=[float(d1[2])]        
file2listobs(infile,rkpc_list,sigma_obs_list,dsigma_obs_list)

sigmaobs=numpy.array(sigma_obs_list)
dsigmaobs=numpy.array(dsigma_obs_list)
rkpc=numpy.array(rkpc_list)
infile.close
#############################
# 1.2 #  Surface brightness #
#############################
infile2=open(obssb,'r')
# this file now excludes the first line which was 16.5024 1.35193 0.0857030.085703

sb_obs_list=[]
dsb_obs_list=[]
rkpcb_obs_list=[]
lrkpcb_obs_list=[]

def file2list(filename,listasb_obs,listadsb_obs,listarkpcb,listalrkpcb):
    for line in filename:
       d2=line.split()
       listasb_obs+=[float(d2[0])]
       listadsb_obs+=[float(d2[1])]
       listarkpcb+=[float(d2[2])]
       listalrkpcb+=[float(d2[3])]     

file2list(infile2,sb_obs_list,dsb_obs_list,rkpcb_obs_list,lrkpcb_obs_list)

sbobs=numpy.array(sb_obs_list)#/1.0e5
dsbobs=numpy.array(dsb_obs_list)#/1.0e5

rkpcb=numpy.array(rkpcb_obs_list)
infile2.close

######################################
# 2 #             SIMULATION DATA    #
######################################
simgrad=open(simulationgrad,'r')

print simgrad
grad_sim_list=[]
dgrad_sim_list=[]

def file2listsim1(filename,listgrad,listdgrad):
    for line in filename:
        dg=line.split()
        listgrad+=[float(dg[0])]
        listdgrad+=[float(dg[1])]
        
file2listsim1(simgrad,grad_sim_list,dgrad_sim_list)

gradsim=numpy.array(grad_sim_list)
dgradsim=numpy.array(dgrad_sim_list)

simgrad.close

######################################

simsigma=open(simulationsigma,'r')
print simsigma
sigma_sim_list=[]
dv_sim_list=[]
rkpc_sim_list=[]

def file2listsim1(filename,listarkpc_sim,listasigma_sim,listadv_sim):
    for line in filename:
        d3=line.split()
        listarkpc_sim+=[float(d3[0])]
        listasigma_sim+=[float(d3[1])]
        listadv_sim+=[float(d3[2])]
        
file2listsim1(simsigma,rkpc_sim_list,sigma_sim_list,dv_sim_list)

sigmasim=numpy.array(sigma_sim_list)
dvsim=numpy.array(dv_sim_list)
rkpcsim=numpy.array(rkpc_sim_list)

simsigma.close

######################################

simsb=open(simulationsb,'r')
sb_sim_list=[]
dsb_sim_list=[]
rkpcsb_sim_list=[]
rlkpcsb_sim_list=[]
buflist=[]
buflist2=[]
def file2listsim2(filename,listasb_sim,listadsb_sim,listarkpcsb_sim,listarlkpcsb_sim,listbuf,listbuf2):
    for line in filename:
        d=line.split()
        listasb_sim+=[float(d[0])]
        listadsb_sim+=[float(d[1])]
        listarkpcsb_sim+=[float(d[2])]
        listarlkpcsb_sim+=[float(d[3])]
        listbuf+=[float(d[4])]
        listbuf2+=[float(d[5])]

file2listsim2(simsb,sb_sim_list,dsb_sim_list,rkpcsb_sim_list,rlkpcsb_sim_list,buflist,buflist2)

sbsim=numpy.array(sb_sim_list)
dsbsim=numpy.array(dsb_sim_list)
rkpcscsim=numpy.array(rkpcsb_sim_list)
rlkpcscsim=numpy.array(rlkpcsb_sim_list)

######################################
# EXCLUDING THE INNERMOST BIN        #
######################################

sbsimshort=numpy.empty(11)
dsbsimshort=numpy.empty(11)
sigmasimshort=numpy.empty(11)
dvsimshort=numpy.empty(11)
for t in xrange(11):
     tt=t+1
     sbsimshort[t]=sbsim[tt]
     dsbsimshort[t]=dsbsim[tt]
     sigmasimshort[t]=sigmasim[tt]
     dvsimshort[t]=dvsim[tt]

simsb.close

##########################################
#  3 #  CALCULATING THE CHI SQUARE       #
##########################################
#bin=binsb+binsigma

out=open(output,'w')
chisigma=sum((sigmasimshort-sigmaobs)**2./(dvsimshort**2.+dsigmaobs**2.))/(binsigma-1)
chisb=sum((sbsimshort-sbobs)**2./(dsbsimshort**2+dsbobs**2.))/(binsb-1)
chigrad=(gradsim-gradobs)**2/(dgradsim**2+dvgradobs**2)/1.0
chi_2=chisigma+chisb+chigrad[0]
likelihood=math.exp(-chi_2/2.0)/2.5066 #sqrt(2*pi)
if math.isnan(chi_2):
    print >>out,'3000.0 3000.0 3000.0 0.0'
else:
    print >>out,chi_2,chisb,chisigma,chigrad[0],likelihood
out.close()
