############################################################################
# Calculates the L.O.S velocity gradient
############################################################################
# Run with python grad.py pa.dat grad.dat 100000 
# Needs outputs of the PA code and the center of mass from the SB code.
############################################################################

import sys
from pytipsy import ifTipsy,TipsyHeader, TipsyDarkParticle,tipsypos, TipsyStarParticle
import math
import numpy

# Takes the arguments        
######################

if len(sys.argv) <=3:       print 'Not enough input files'
simulationsigma=sys.argv[1]  
output=sys.argv[2]   
n_d=sys.argv[3]
#####################

simsigma=open(simulationsigma,'r')

print simsigma
x_list=[]
y_list=[]
z_list=[]
m_list=[]
eps_list=[]
vx_list=[]
vy_list=[]
vz_list=[]
vlos_list=[]
rmaj_list=[]
rmajamin_list=[]
rmin_list=[]

print 'aa'
def file2listsim1(filename,listm,listx,listy,listz,listvx,listvy,listvz,listeps,listrmaj,listrmajamin,listrmin,listvlos):
    for line in filename:
        d3=line.split()
        listm+=[float(d3[0])]
        listx+=[float(d3[1])]
        listy+=[float(d3[2])]
        listz+=[float(d3[3])]
        listvx+=[float(d3[4])]
        listvy+=[float(d3[5])]
        listvz+=[float(d3[6])]
        listeps+=[float(d3[7])]
        listrmaj+=[float(d3[8])]
        listrmajamin+=[float(d3[9])]
        listrmin+=[float(d3[10])]
        listvlos+=[float(d3[11])]

       
file2listsim1(simsigma,m_list,x_list,y_list,z_list,vx_list,vy_list,vz_list,eps_list,rmaj_list,rmajamin_list,rmin_list,vlos_list)

m=numpy.array(m_list)
eps=numpy.array(eps_list)
x=numpy.array(x_list)
y=numpy.array(y_list)
z=numpy.array(z_list)
vx=numpy.array(vx_list)
vy=numpy.array(vy_list)
vz=numpy.array(vz_list)
vlos=0.655*numpy.array(vlos_list)
rmaj=numpy.array(rmaj_list)
rmin=numpy.array(rmin_list)
print 'vel1',vx[0],vy[0],vz[0],vlos

#Sorting with the position on the major axis
##############################################

simsigma.close()
rm=rmaj[rmaj.argsort()];
rminn=rmin[rmaj.argsort()];
vl=vlos[rmaj.argsort()];
x_s=x[rmaj.argsort()];
y_s=y[rmaj.argsort()];
z_s=z[rmaj.argsort()];
vx_s=vx[rmaj.argsort()];
vy_s=vy[rmaj.argsort()];
vz_s=vz[rmaj.argsort()];
eps=eps[rmaj.argsort()];
m_s=m[rmaj.argsort()];
vll=vlos[rmaj.argsort()];
print 'vel2',vx_s[0],vy_s[0],vz_s[0],vl

# Reading  the centre of mass     #
######################################

com=open("com_file.dat",'r')

line=com.readline()
d=line.split()
xc=float(d[0])
yc=float(d[1])
zc=float(d[2])
vxc=float(d[3])
vyc=float(d[4])
vzc=float(d[5])
dlosc=float(d[6])
xveclosc=float(d[7])
yveclosc=float(d[8])
zveclosc=float(d[9])
com.close()

print 'centre of mass read by the meanvel.py',xc,yc,zc,vxc,vyc,vzc,dlosc,xveclosc,yveclosc,zveclosc


#  Relative distance and velocity of each star from the center of Carina and    
############################################################################
# ic=star to center #
#####################
xic=x_s-xc
yic=y_s-yc
zic=z_s-zc
vxic=(vx_s-vxc)*0.655
vyic=(vy_s-vyc)*0.655
vzic=(vz_s-vzc)*0.655
rsic=(xic**2+yic**2+zic**2)**0.5
vsic=(vxic**2+vyic**2+vzic**2)**0.5
print 'vel3',vxic,vyic,vzic
q=0

#########m#
# Finding the components of the e2 vector on the plane
#########################################################
zvec2=0.0    
rat=-yveclosc/xveclosc 
yvec2=1/(1+rat**2)**0.5;
xvec2=(1-yvec2**2)**0.5;
xvec3=zvec2*yveclosc-zveclosc*yvec2
yvec3=xvec2*zveclosc-xveclosc*zvec2
zvec3=yvec2*xveclosc-yveclosc*xvec2

xf=xic*xveclosc+yic*yveclosc+zic*zveclosc 
yf=xic*xvec2+yic*yvec2+zic*zvec2  
zf=xic*xvec3+yic*yvec3+zic*zvec3  
print 'vxic[i],xveclosc,vyic[i],yveclosc,vzic[i],zveclosc',vxic[0],xveclosc,vyic[0],yveclosc,vzic[0],zveclosc

vlos_calc=vxic*xveclosc+vyic*yveclosc+vzic*zveclosc
avlos=(vlos_calc**2)**0.5
print 'vlos_calc',numpy.mean(vlos_calc)#/vll
r2d=(yf**2+zf**2)**0.5
r3d=(xf**2+yf**2+zf**2)**0.5

rbinmin=numpy.empty(2)
rbinmax=numpy.empty(2)
vbarm=numpy.empty(2)
dvbarm=numpy.empty(2)
vbarkm=numpy.empty(2)
dvbarkm=numpy.empty(2)
nbinm=numpy.empty(2)

for j in xrange(2):
    if j==0:
        rbinmin[j]=-1.48085
        rbinmax[j]=-0.504713
    else:
        rbinmax[j]=1.48085
        rbinmin[j]=0.504713
    vbarm[j]=0.0
    nbinm[j]=0
    for f in xrange(100000):
        if(rm[f]>rbinmin[j] and rm[f]<rbinmax[j]):
            vbarm[j]=vbarm[j]+vlos_calc[f]
            nbinm[j]=nbinm[j]+1
    vbarkm[j]=vbarm[j]/nbinm[j]
    dvbarkm[j]=vbarkm[j]/(nbinm[j])**0.5
    print 'here',nbinm[j],rbinmax[j],j
r_ave=numpy.array(bin)
r_ave=[-0.9927,0.9927]
r_plt=[-1.0,1.1]
vgrad=vbarkm[1]-vbarkm[0]
dvgrad=abs(dvbarkm[1])+abs(dvbarkm[0])

fgrad=open(output,'w')
print>> fgrad,vgrad,dvgrad
fgrad.close() 

#radial vel calc for the centre
rr=(xc**2+yc**2+zc**2)**0.5
xr_hat=(xc-8.5)/rr
yr_hat=yc/rr
zr_hat=zc/rr
vxcr=vxc*xr_hat
vycr=vyc*yr_hat
vzcr=vzc*zr_hat
vrr=vxcr+vycr+vzcr
print rr,vxc,vxcr,vrr

