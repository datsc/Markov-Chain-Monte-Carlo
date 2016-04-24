###############################################################
# Converting to PKDGRAV units and halo+bulge files     #
###################################################################
# Run with python units1public.py 200000 100000 100000 oroutcirc
###################################################################
# INPUT
# NEMO UNITS:
#   G = 1
#   unit_t = Gyrs
#   unit_v = .9777746km/s
#   unit_l = kpc
#   unit_m = 222288  Msun
#
# OUTPUT
# OUR PKDGRAV UNITS
#   G = 1.0
#   unit_t = 1.4908Gyrs
#   unit_v = 0.655km/s       
#   unit_l = kpc
#   unit_m = 10^5 Msun
#######################################################
##############################
# Importing python libraries #
##############################
import sys
from pytipsy import ifTipsy,TipsyHeader, TipsyDarkParticle,tipsypos, TipsyStarParticle
import pytipsy
import math
import numpy 
############################################

n=int(sys.argv[1]) 
n_d=int(sys.argv[2])
n_s=int(sys.argv[3])
orbp=sys.argv[4]


OBSdata=open(orbp,'r')
xco_list=[]
yco_list=[]
zco_list=[]
vxco_list=[]
vyco_list=[]
vzco_list=[]

def file2listtheo(filename,listxco,listyco,listzco,listvxco,listvyco,listvzco):
    for line in filename:
        d=line.split()
        listxco+=[float(d[0])]
        listyco+=[float(d[1])]
        listzco+=[float(d[2])]
        listvxco+=[float(d[3])]
        listvyco+=[float(d[4])]
        listvzco+=[float(d[5])]

file2listtheo(OBSdata,xco_list,yco_list,zco_list,vxco_list,vyco_list,vzco_list)

xco=numpy.array(xco_list)
yco=numpy.array(yco_list)
zco=numpy.array(zco_list)
vxco=numpy.array(vxco_list)
vyco=numpy.array(vyco_list)
vzco=numpy.array(vzco_list)
xco=xco[0]
yco=yco[0]
zco=zco[0]
vxco=vxco[0]
vyco=vyco[0]
vzco=vzco[0]
OBSdata.close()


halof=open('haloo','r')

m_d_list=[]
x_d_list=[]
y_d_list=[]
z_d_list=[]
vx_d_list=[]
vy_d_list=[]
vz_d_list=[]
def file2list(filename,listadm,listadx,listady,listadz,listadvx,listadvy,listadvz):
    for line in filename:
        d=line.split()
        listadm+=[float(d[0])]
        listadx+=[float(d[1])]
        listady+=[float(d[2])]
        listadz+=[float(d[3])]
        listadvx+=[float(d[4])]
        listadvy+=[float(d[5])]
        listadvz+=[float(d[6])]
file2list(halof,m_d_list,x_d_list,y_d_list,z_d_list,vx_d_list,vy_d_list,vz_d_list)
md=numpy.array(m_d_list)
xd=numpy.array(x_d_list)
yd=numpy.array(y_d_list)
zd=numpy.array(z_d_list)
vxd=numpy.array(vx_d_list)
vyd=numpy.array(vy_d_list)
vzd=numpy.array(vz_d_list)

bulgef=open('bulgeo','r')
 
m_b_list=[]
x_b_list=[]
y_b_list=[]
z_b_list=[]
vx_b_list=[]
vy_b_list=[]
vz_b_list=[]
def file2list(filename,listabm,listabx,listaby,listabz,listabvx,listabvy,listabvz):
    for line in filename:
        d=line.split()
        listabm+=[float(d[0])]
        listabx+=[float(d[1])]
        listaby+=[float(d[2])]
        listabz+=[float(d[3])]
        listabvx+=[float(d[4])]
        listabvy+=[float(d[5])]
        listabvz+=[float(d[6])]
file2list(bulgef,m_b_list,x_b_list,y_b_list,z_b_list,vx_b_list,vy_b_list,vz_b_list)
mb=numpy.array(m_b_list)
xb=numpy.array(x_b_list)
yb=numpy.array(y_b_list)
zb=numpy.array(z_b_list)
vxb=numpy.array(vx_b_list)
vyb=numpy.array(vy_b_list)
vzb=numpy.array(vz_b_list)

#SORTING
rdic=(xd**2+yd**2+zd**2)**0.5
ric=(xb**2+yb**2+zb**2)**0.5
r_dd=rdic[rdic.argsort()]
m_dd=md[rdic.argsort()]
x_dd=xd[rdic.argsort()]
y_dd=yd[rdic.argsort()]
z_dd=zd[rdic.argsort()]
vx_dd=vxd[rdic.argsort()]
vy_dd=vyd[rdic.argsort()]
vz_dd=vzd[rdic.argsort()]
r_sd=ric[ric.argsort()]
m_sd=mb[ric.argsort()]
x_sd=xb[ric.argsort()]
y_sd=yb[ric.argsort()]
z_sd=zb[ric.argsort()]
vx_sd=vxb[ric.argsort()]
vy_sd=vyb[ric.argsort()]
vz_sd=vzb[ric.argsort()]

#SHIFT 
xds=x_dd+xco
yds=y_dd+yco
zds=z_dd+zco
vxds=vx_dd+vxco*1.524827
vyds=vy_dd+vyco*1.524827
vzds=vz_dd+vzco*1.524827

xbs=x_sd+xco
ybs=y_sd+yco
zbs=z_sd+zco
vxbs=vx_sd+vxco*1.524827
vybs=vy_sd+vyco*1.524827
vzbs=vz_sd+vzco*1.524827

outh=open('halo','w')
outb=open('bulge','w')


print >> outh ,n_d,0.0
print 'number of particles n=',n

for i in xrange(n_d):
    print >> outh,m_dd[i],xds[i],yds[i],zds[i],vxds[i],vyds[i],vzds[i]
outh.close()

print >> outb, n_s,0.0
for j in xrange(n_s):
    print >> outb,m_sd[j],xbs[j],ybs[j],zbs[j],vxbs[j],vybs[j],vzbs[j]
outb.close()
