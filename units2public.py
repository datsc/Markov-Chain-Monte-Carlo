##########################################################
# Shifting to new orbit
##########################################################
# Run with python units2public.py 200000 100000 100000 orbpxpast0.dat comeq5.dat
# Inputs are haloo and bulgeo  # Outputs are halo and bulge
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
orbeq=sys.argv[5]

EQdata=open(orbeq,'r')
xcoeq_list=[]
ycoeq_list=[]
zcoeq_list=[]
vxcoeq_list=[]
vycoeq_list=[]
vzcoeq_list=[]
rc_list=[]

def file2listtheo(filename,listxcoeq,listycoeq,listzcoeq,listvxcoeq,listvycoeq,listvzcoeq,listrccoeq):
    for line in filename:
        d=line.split()
        listxcoeq+=[float(d[0])]
        listycoeq+=[float(d[1])]
        listzcoeq+=[float(d[2])]
        listvxcoeq+=[float(d[3])]
        listvycoeq+=[float(d[4])]
        listvzcoeq+=[float(d[5])]
	listrccoeq+=[float(d[6])]

file2listtheo(EQdata,xcoeq_list,ycoeq_list,zcoeq_list,vxcoeq_list,vycoeq_list,vzcoeq_list,rc_list)

xcoeq=numpy.array(xcoeq_list)
ycoeq=numpy.array(ycoeq_list)
zcoeq=numpy.array(zcoeq_list)
vxcoeq=numpy.array(vxcoeq_list)
vycoeq=numpy.array(vycoeq_list)
vzcoeq=numpy.array(vzcoeq_list)
rccoeq=numpy.array(rc_list)
xcoeq=xcoeq[0]
ycoeq=ycoeq[0]
zcoeq=zcoeq[0]
vxcoeq=vxcoeq[0]
vycoeq=vycoeq[0]
vzcoeq=vzcoeq[0]
rccoeq=rccoeq[0]
EQdata.close()
print 'xco of the relaxationrun',xcoeq,ycoeq,zcoeq,vxcoeq,vycoeq,vzcoeq

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
print 'xco in the code',xco,yco,zco,vxco,vyco,vzco

deltaxc=xcoeq-xco
deltayc=ycoeq-yco
deltazc=zcoeq-zco
deltavxc=vxcoeq-vxco
deltavyc=vycoeq-vyco
deltavzc=vzcoeq-vzco
print 'deltaxc',deltaxc
#################
# FIRST HALO    #
#################
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

#################
# FIRST BULGE    #
#################
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
print 'max, min',max(xd), min(xd)
#SORTING
rdic=((xd-xcoeq)**2+(yd-ycoeq)**2+(zd-zcoeq)**2)**0.5
ric=((xb-xcoeq)**2+(yb-ycoeq)**2+(zb-zcoeq)**2)**0.5
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

#####################
# CENTRE
#####################

xds=x_dd-xcoeq+xco
yds=y_dd-ycoeq+yco
zds=z_dd-zcoeq+zco
#vxds=vx_dd-vxcoeq*1.524827+vxco*1.524827
vxds=vx_dd-vxcoeq*1.5246+vxco*1.5246
vyds=vy_dd-vycoeq*1.5246+vyco*1.5246
vzds=vz_dd-vzcoeq*1.5246+vzco*1.5246
print 'x_dd before and after shift',vxds[0],vx_dd[0]

xbs=x_sd-xcoeq+xco
ybs=y_sd-ycoeq+yco
zbs=z_sd-zcoeq+zco
vxbs=vx_sd-vxcoeq*1.5246+vxco*1.5246
vybs=vy_sd-vycoeq*1.5246+vyco*1.5246
vzbs=vz_sd-vzcoeq*1.5246+vzco*1.5246

#############################
# WRITING THE NEW FILES     #
#############################
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
