###################################################################
# Reads in tipsy files and writes out ascii
###################################################################
# Run with python binaryfile.bin asciifile.dat
###################################################################
###########
# TIPSY   #
###########
#TIPSY FILES:
#   header
#   gas
#   dark
#   stars
#
#FOR EACH PARTICLE:
#   mass
#   pos[0],pos[1],pos[2] (x,y,z)
#   vel[0],vel[1],vel[2] (V_x,V_y,V_z)
#   phi (potential)
#   density
##############################
# Importing python libraries #
##############################
import sys
from pytipsy import ifTipsy,TipsyHeader, TipsyDarkParticle,tipsypos, TipsyStarParticle
import pytipsy
import math
import numpy 

filename=sys.argv[1]
output=sys.argv[2]  

h   = TipsyHeader() 
inp = ifTipsy(filename, "standard") 
inp.read(h) 
print "h=",h
print "nBodies=", h.h_nBodies

print "nDark=", h.h_nDark
print "nStars=", h.h_nStar

n_s=h.h_nStar
n_d=h.h_nDark
print 'n_d=',n_d
print 'n_s=',n_s
d=TipsyDarkParticle()
s=TipsyStarParticle()
m_d_list=[]
eps_d_list=[]
m_s_list=[]
eps_s_list=[]
x_d_list=[]
y_d_list=[]
z_d_list=[]
vx_d_list=[]
vy_d_list=[]
vz_d_list=[]
x_s_list=[]
y_s_list=[]
z_s_list=[]
vx_s_list=[]
vy_s_list=[]
vz_s_list=[]
inp.seekg(tipsypos(tipsypos.dark, 0))
for i in xrange(n_d): 
    inp.read(d)
    x_d_list+=[d.pos[0]]
    y_d_list+=[d.pos[1]]
    z_d_list+=[d.pos[2]]
    vx_d_list+=[d.vel[0]]
    vy_d_list+=[d.vel[1]]
    vz_d_list+=[d.vel[2]]
    m_d_list+=[d.mass]
    eps_d_list+=[d.eps]
inp.seekg(tipsypos(tipsypos.star, 0))
for j in xrange(n_s): 
    inp.read(s)
    x_s_list+=[s.pos[0]]
    y_s_list+=[s.pos[1]]
    z_s_list+=[s.pos[2]]
    vx_s_list+=[s.vel[0]]
    vy_s_list+=[s.vel[1]]
    vz_s_list+=[s.vel[2]]
    m_s_list+=[s.mass]
    eps_s_list+=[s.eps]

x_d=numpy.array(x_d_list)
y_d=numpy.array(y_d_list)
z_d=numpy.array(z_d_list)
m_di=numpy.array(m_d_list)
vx_d=numpy.array(vx_d_list)
vy_d=numpy.array(vy_d_list)
vz_d=numpy.array(vz_d_list)
x_s=numpy.array(x_s_list)
y_s=numpy.array(y_s_list)
z_s=numpy.array(z_s_list)
m_si=numpy.array(m_s_list)
vx_s=numpy.array(vx_s_list)
vy_s=numpy.array(vy_s_list)
vz_s=numpy.array(vz_s_list)
eps_s=numpy.array(eps_s_list)
eps_d=numpy.array(eps_d_list)

out=open(output,'w')
for i in xrange(0,n_d):
    print >>out, m_di[i],x_d[i],y_d[i],z_d[i],vx_d[i],vy_d[i],vz_d[i],eps_d[i]
for j in xrange(0,n_s):
    print >>out, m_si[j],x_s[j],y_s[j],z_s[j],vx_s[j],vy_s[j],vz_s[j],eps_s[j]
out.close()


