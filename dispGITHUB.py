############################################################################
# Calculates and plots the projected velocity dispersion profile with R and sigma
# using the tipsy input.
# Run with python dispGITHUB.py final.00060 simdisp.dat simdisp.png 
#####
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
#######################################################
#######################################################
# READ AND OUTPUTS PKDGRAV
# UNITS:
# READING
#   G = DOUBLE(1.)                      ;Gravitational constant
#   unit_t = DOUBLE(0.0098)             ;Unit time (Giga-years)
#   unit_v = DOUBLE(100)                ;Unit velocity (km/s)
#   unit_l = DOUBLE(3.086e19)           ;Unit length (metres [1kpc])
#   unit_m = DOUBLE(2.325e9)            ;Unit mass (Msun)
#
#######################################################
##############################
# Importing python libraries #
##############################
from pytipsy import ifTipsy,TipsyHeader,TipsyDarkParticle,tipsypos,TipsyStarParticle
import sys
import pytipsy
import math
import numpy 
import scipy
from scipy import special

###################################
# Takes the arguments
###################################
#if len(sys.argv) <=1:
 #   filename=raw_input("File name: ")
if len(sys.argv) <=2:
    print 'Not enough input files'
filename=sys.argv[1]  #the final output of pkdgrav
output=sys.argv[2]   #the bulge in ascii
outfig1=sys.argv[3]   #bulge figure


######################################
# Assigning sun's coordinates        #                     
######################################
xsun=8.5;ysun=0.0 ;zsun=0.0; vxsun=0.0; vysun=200.0; vzsun=0.0; sol_m=1e5*1.98892e30; G=6.67428*10**(-11)
#sol_m=2.325e9*1.989e30

###################################################
# Reading data from tipsy file and making arrays  #
###################################################
h   = TipsyHeader()                      # creates a location for the header and all the rest
inp = ifTipsy(filename, "standard")      # opens the file
inp.read(h)                              # reads the header

n_s=h.h_nStar
n_d=h.h_nDark
n_b=h.h_nBodies
print 'n_d=',n_d; print 'n_s=',n_s; print 'n_b=',n_b

d=TipsyDarkParticle()
s= TipsyStarParticle()
m_d_list=[];    m_s_list=[];
x_d_list=[];    x_s_list=[]; 
y_d_list=[];    y_s_list=[]; 
z_d_list=[];    z_s_list=[]; 
vx_d_list=[];   vx_s_list=[];
vy_d_list=[];   vy_s_list=[];
vz_d_list=[];   vz_s_list=[];

#DM
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

#STARS
inp.seekg(tipsypos(tipsypos.star, 0))
for i in xrange(n_s): 
    inp.read(s)
    x_s_list+=[s.pos[0]]
    y_s_list+=[s.pos[1]]
    z_s_list+=[s.pos[2]]
    vx_s_list+=[s.vel[0]]
    vy_s_list+=[s.vel[1]]
    vz_s_list+=[s.vel[2]]
    m_s_list+=[s.mass]

#MAKING EMPTY LISTS
x_d=numpy.array(x_d_list);                x_s=numpy.array(x_s_list);       
y_d=numpy.array(y_d_list);                y_s=numpy.array(y_s_list);       
z_d=numpy.array(z_d_list);                z_s=numpy.array(z_s_list);       
vx_d=numpy.array(vx_d_list);              vx_s=numpy.array(vx_s_list);     
vy_d=numpy.array(vy_d_list);              vy_s=numpy.array(vy_s_list);     
vz_d=numpy.array(vz_d_list);              vz_s=numpy.array(vz_s_list);     
m_di=numpy.array(m_d_list);               m_si=numpy.array(m_s_list);

######################################
# Reading  the centre of mass     #
######################################

com=open("com_file.dat",'r')

#def file2listtheo(filename,listxc,listyc,listzc,listvxc,listvyc,listvzc,listdlos,listxveclos,listyveclos,listzveclos):
line=com.readline()
d=line.split()
xc=float(d[0])
yc=float(d[1])
zc=float(d[2])
vxc=float(d[3])
vyc=float(d[4])
vzc=float(d[5])
dlos=float(d[6])
xveclosc=float(d[7])
yveclosc=float(d[8])
zveclosc=float(d[9])
com.close()

print 'centre of mass',xc,yc,zc,vxc,vyc,vzc
##################################################################################
#  Calculating the distance and the unit vectors of the line of sight (c=center) #
##################################################################################
#dlosc=((xc-xsun)**2+(yc-ysun)**2+(zc-zsun)**2)**0.5
#xveclosc=(xc-xsun)/dlosc
#yveclosc=(yc-ysun)/dlosc
#zveclosc=(zc-zsun)/dlosc


#####################################################################################
#  Relative distance and velocity of each star from the center of Carina and        #  
#####################################################################################
# ic=star to center #
#####################
xic=x_s-xc
yic=y_s-yc
zic=z_s-zc
vxic=vx_s-vxc
vyic=vy_s-vyc
vzic=vz_s-vzc
rsic=(xic**2+yic**2+zic**2)**0.5
vsic=(vxic**2+vyic**2+vzic**2)**0.5

# dic=star to center #
#####################
#print x_d,y_d,z_d
xdic=x_d-xc
ydic=y_d-yc
zdic=z_d-zc
vxdic=vx_d-vxc
vydic=vy_d-vyc
vzdic=vz_d-vzc
rdic=(xdic**2+ydic**2+zdic**2)**0.5
vdic=(vxdic**2+vydic**2+vzdic**2)**0.5

#############################
#  Projected distance(star) #
#############################
r_si=(rsic**2-((xic*xveclosc)+(yic*yveclosc)+(zic*zveclosc))**2)**0.5
vlossi=(vxic*xveclosc+vyic*yveclosc+vzic*zveclosc)
v_bars=sum(vlossi)/n_d

#############################
#  Projected distance(dm)   #
#############################
r_di=(rdic**2-((xdic*xveclosc)+(ydic*yveclosc)+(zdic*zveclosc))**2)**0.5
vlosdi=(vxdic*xveclosc+vydic*yveclosc+vzdic*zveclosc)
v_bard=sum(vlosdi)/n_d

########################
#  Sorting the arrays  #
########################
r_s=r_si[r_si.argsort()]
m_s=m_si[r_si.argsort()]
r_d=r_di[r_di.argsort()]
m_d=m_di[r_di.argsort()]
xs=x_s[r_si.argsort()]
ys=y_s[r_si.argsort()]
zs=z_s[r_si.argsort()]
vxs=vx_s[r_si.argsort()]
vys=vy_s[r_si.argsort()]
vzs=vz_s[r_si.argsort()]
ms=[r_si.argsort()]
vloss=vlossi[r_si.argsort()]

############################
#  Calculating bin size    #
############################
bin=12 

maxx=max(r_s) 
minn=min(r_s)
maxxd=max(r_d) 
minnd=min(r_d)

lrmax=math.log10(maxx)
lrmin=math.log10(minn)
dr=(lrmax-lrmin)/float(bin)
lrmaxd=math.log10(maxxd)
lrmind=math.log10(minnd)
drd=(lrmaxd-lrmind)/float(bin)

############################
# Calculating sigma_p      #
############################
o=0
p=0
n=numpy.zeros(bin)
v_bark=numpy.zeros(bin)
vbark=numpy.zeros(bin)
sigma=numpy.zeros(bin)
rn=numpy.zeros(bin)
vp=numpy.zeros(bin)
v_d=numpy.zeros(n_d)
v_s=numpy.zeros(n_s)
dsigma=numpy.zeros(bin)

arcmin=9.155397604e17/3.085e19


#defining the bin radii from Car_disp_obs.dat# by chosing the limits half-way through the bins,(eg.(104.4+56)/2=80.2 
# OLD
#r_obs=[0.0802,0.1229,0.1531,0.18325,0.2126,0.2429,0.2866,0.36125,0.59635,maxx]
#r_obs=[56.4,105.1,138.9,169.7,199.9,229.0,262.2,321.4,426.1,828.8]old!#r_n=numpy.array(r_obs)/1000.0
##############
r_l=numpy.zeros(bin)
r_obs=[0.0857038,0.114725,0.151126,0.174134,0.197232,0.220339,0.249179,0.284005,0.330016,0.391652,0.504713,1.48085] #this is the max r in the bin
r_ave=[0.059246905306,0.10235957160,0.13178707523,0.16221119215,0.18688856515,0.21132620291,0.23558852524,0.26508560119,0.30666115969,0.35824041858,0.43622672025,0.73376536189]
r_l=numpy.array(r_obs)
r_ave_obs=numpy.array(r_ave)
for k in xrange(bin):
    n[k]=0
    while r_s[o]<=r_l[k]: 
        vbark[k]+=vloss[o]
        o+=1
        n[k]+=1
        if o==n_s:
            break
    v_bark[k]=vbark[k]/(n[k])
print 'inside',v_bark,vloss
print r_l


for q in xrange(bin):
    n[q]=0
    while r_s[p]<=r_l[q]:     
        v_s[p]=(vloss[p]-v_bark[q])**2
        n[q]+=1
        vp[q]+=v_s[p]
        p+=1
        if p==n_s:
            break
        sigma[q]=(vp[q]/n[q])**0.5 
    dsigma[q]=sigma[q]/numpy.sqrt(n[q])
oran=sum(v_bark)/bin

print ' n particles',n
print 'DONT forget that at the very end you divide the sigma by 1.5 for converting back from pkdgrav units!!!!'
veldat=open(output,'w')
kpc=3.085*10**19#m
for q in xrange(bin):
    #print >>veldat,r_l[q],r_l[q],sigma[q],sigma[q]*0.655,dsigma[q]*0.655
    print >>veldat,r_l[q],sigma[q]*0.655875368,dsigma[q]*0.655875368
print 'data in veldisp.dat is  r[kpc],r[m],v[0.655km/s] and v[km/s]'


######################################
# Including the theoretical data     #
######################################

Theodata=open('obsdisp146noisy180.dat','r')#ONLY for MODEL

rth1_list=[]
rth2_list=[]
vth1_list=[]
dvth2_list=[]

def file2listtheo(filename,listrth1,listvth1,listdvth2):
    for line in filename:
        d=line.split()
        listrth1+=[float(d[0])]
        listvth1+=[float(d[1])]
        listdvth2+=[float(d[2])]

file2listtheo(Theodata,rth1_list,vth1_list,dvth2_list)  

rth1=numpy.array(rth1_list)  #This is the limit, in the figure I use average
vth1=numpy.array(vth1_list)
dvth2=numpy.array(dvth2_list)
Theodata.close()


################
# Plotting     #
################
from pylab import figure,plot,show,title,legend,xlabel,ylabel,savefig,semilogx, loglog,errorbar
outf1=open(outfig1,'w')
figure()
title('Projected velocity dispersion profile -General!')
xlabel('r (kpc)')
ylabel('$\sigma$ (km/s)')
p_MWs5=errorbar(r_ave_obs,sigma*0.655875368, yerr=dsigma*0.655,xerr=None,
                ecolor=None, elinewidth=None, capsize=3,
                barsabove=True, lolims=False, uplims=False,
                xlolims=False, xuplims=False)
s2=errorbar(rth1,vth1, yerr=dvth2,xerr=None,
                ecolor=None, elinewidth=None, capsize=3,
                barsabove=True, lolims=False, uplims=False,
                xlolims=False, xuplims=False)
savefig(outf1)
outf1.close
