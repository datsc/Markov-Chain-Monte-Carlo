###################################################################
# Calculating the surface brightness profile
###############################################################################
# Run with python sbpublic.py binary.00060 simb.dat simh.dat simb.png simh.png
###############################################################################
###########
# TIPSY   #
###########
#TIPSY FILES:1
#   header + #   gas + #   dark + #   stars
#FOR EACH PARTICLE:
#   mass
#   pos[0],pos[1],pos[2] (x,y,z)
#   vel[0],vel[1],vel[2] (V_x,V_y,V_z)
#   eps
# no  phi (potential) # no  density
#######################################################
# UNITS:
#   G = DOUBLE(1.)                      ;Gravitational constant
#   unit_t = DOUBLE(1.4907Gyrs)         ;Unit time (Giga-years)
#   unit_v = DOUBLE(1)                  ;Unit velocity (kpc/1.4907Gyrs)
#   unit_l = DOUBLE(1kpc)               ;Unit length 
#   unit_m = DOUBLE(1e5msun)            ;Unit mass 
#######################################################
##############################
# Importing python libraries #
##############################
import sys
from pytipsy import ifTipsy,TipsyHeader, TipsyDarkParticle,tipsypos, TipsyStarParticle
import pytipsy
import math
import numpy
import matplotlib
from matplotlib import pylab
from pylab import figure, plot, show,semilogx,xscale,yscale,errorbar, title,savefig, legend, xlabel, ylabel, close, xlim, ylim, semilogy,loglog,axis

######################################
# Assigning sun's coordinates        #                     
######################################
xsun=8.5;ysun=0.0 ;zsun=0.0; vxsun=0.0; vysun=200.0; vzsun=0.0; sol_m=1.0;
###############################
# Takes the arguments        #
###############################
if len(sys.argv) <=4:
    print 'Not enough input files'
filename=sys.argv[1]  #the final output of pkdgrav 
output1=sys.argv[2]   #the bulge in ascii
output2=sys.argv[3]   #the halo with ascii
outfig1=sys.argv[4]   #bulge figure
outfig2=sys.argv[5]   #halo figure
scllngdm=1480.8       #the outer radius of the outermost dispersion bin

###################################################
# Reading data from tipsy file and making arrays  #
###################################################
h   = TipsyHeader()                      
inp = ifTipsy(filename, "standard")      
inp.read(h)                              

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
m_di=sol_m*numpy.array(m_d_list);         m_si=sol_m*numpy.array(m_s_list);

inp.close

######################################
# Calculating the centre of mass     #
######################################
# INITIAL VALUES
####################
corenold=90000
coren=corenold
mtotoldcore=coren*m_si
xcoldcore=numpy.empty(coren)
ycoldcore=numpy.empty(coren)
zcoldcore=numpy.empty(coren)
vxcoldcore=numpy.empty(coren)
vycoldcore=numpy.empty(coren)
vzcoldcore=numpy.empty(coren)

for cc in xrange(corenold):
    if cc==0:
        xcoldcore[0]=x_s[0]
        ycoldcore[0]=y_s[0]
        zcoldcore[0]=z_s[0] 
        vxcoldcore[0]=vx_s[0]
        vycoldcore[0]=vy_s[0]
        vzcoldcore[0]=vz_s[0]
    if cc>=1:
        xcoldcore[cc]+=xcoldcore[cc-1]+x_s[cc]
        ycoldcore[cc]+=ycoldcore[cc-1]+y_s[cc]
        zcoldcore[cc]+=zcoldcore[cc-1]+z_s[cc]
        vxcoldcore[cc]+=vxcoldcore[cc-1]+vx_s[cc]
        vycoldcore[cc]+=vycoldcore[cc-1]+vy_s[cc]
        vzcoldcore[cc]+=vzcoldcore[cc-1]+vz_s[cc]    

mold=mtotoldcore 
vxcold=vxcoldcore[coren-1]/coren 
vycold=vycoldcore[coren-1]/coren 
vzcold=vzcoldcore[coren-1]/coren
xcold=xcoldcore[coren-1]/coren 
ycold=ycoldcore[coren-1]/coren
zcold=zcoldcore[coren-1]/coren

# 2- RELATIVE VALUES FOR INDIVIDUAL PARTICLES
##############################################
# STARS
xicold=x_s-xcold; yicold=y_s-ycold; zicold=z_s-zcold
ricold=(xicold**2+yicold**2+zicold**2)**0.5
# DM
xdicold=x_d-xcold; ydicold=y_d-ycold; zdicold=z_d-zcold
rdicold=(xdicold**2+ydicold**2+zdicold**2)**0.5


# 3- SORTING EVERYTHING WITH RESPECT TO DISTANCE FROM THE CENTRE
#################################################################
r_dinold=rdicold[rdicold.argsort()];  r_snold=ricold[ricold.argsort()]
m_dnold=m_di[rdicold.argsort()];      m_snold=m_si[ricold.argsort()]
x_dnold=x_d[rdicold.argsort()];    x_snold=x_s[ricold.argsort()]
y_dnold=y_d[rdicold.argsort()];    y_snold=y_s[ricold.argsort()]
z_dnold=z_d[rdicold.argsort()];    z_snold=z_s[ricold.argsort()]
vx_dnold=vx_d[rdicold.argsort()];    vx_snold=vx_s[ricold.argsort()]
vy_dnold=vy_d[rdicold.argsort()];    vy_snold=vy_s[ricold.argsort()]
vz_dnold=vz_d[rdicold.argsort()];    vz_snold=vz_s[ricold.argsort()]

# 4- LIMITS FOR THE FIRST ITERATION
#################################################
corenold=int(corenold)
rmas=r_snold[corenold]
corenn=int(0.98*corenold)
indls=corenn

n_sold=n_s
n_s=indls
mrlc=r_snold[corenn]
for mnm in xrange(n_d):
    if r_dinold[mnm]<=mrlc:
        mnm+=1
    elif r_dinold[mnm]>mrlc:
        indl=mnm
        break
n_dold=n_d
n_d=indl


# 5 - NEW ARRAYS FOR FEWER NUMBER OF PARTICLES
#################################################
print indls
x_ss=numpy.empty(indls);    x_ds=numpy.empty(indl)
y_ss=numpy.empty(indls);    y_ds=numpy.empty(indl)
z_ss=numpy.empty(indls);    z_ds=numpy.empty(indl)
vx_ss=numpy.empty(indls);   vx_ds=numpy.empty(indl)
vy_ss=numpy.empty(indls);   vy_ds=numpy.empty(indl)
vz_ss=numpy.empty(indls);   vz_ds=numpy.empty(indl)
m_ss=numpy.empty(indls);    m_ds=numpy.empty(indl)
r_ss=numpy.empty(indls);    r_ds=numpy.empty(indl)
ric=numpy.empty(indls);     r_din=numpy.empty(indl)

for fs in xrange(indls):
    x_ss[fs]=x_snold[fs]
    y_ss[fs]=y_snold[fs]
    z_ss[fs]=z_snold[fs]
    vx_ss[fs]=vx_snold[fs]
    vy_ss[fs]=vy_snold[fs]
    vz_ss[fs]=vz_snold[fs]
    m_ss[fs]=m_snold[fs]

for fn in xrange(indl):
    x_ds[fn]=x_dnold[fn]
    y_ds[fn]=y_dnold[fn]
    z_ds[fn]=z_dnold[fn]   
    vx_ds[fn]=vx_dnold[fn]
    vy_ds[fn]=vy_dnold[fn]
    vz_ds[fn]=vz_dnold[fn]
    m_ds[fn]=m_dnold[fn]

# 6 - ITERATION
##########################
myiter=0
for myiter in xrange(150): 
    m=sum(m_ss)#+sum(m_ds)
    xc = 0
    xc=(sum(x_ss)*m_si[0])/m
    yc=(sum(y_ss)*m_si[0])/m
    zc=(sum(z_ss)*m_si[0])/m
    vxc=(sum(vx_ss)*m_si[0])/m
    vyc=(sum(vy_ss)*m_si[0])/m
    vzc=(sum(vz_ss)*m_si[0])/m
    xic=x_ss-xc
    yic=y_ss-yc
    zic=z_ss-zc
    vxic=vx_ss-vxc
    vyic=vy_ss-vyc
    vzic=vz_ss-vzc
    ric=(xic**2+yic**2+zic**2)**0.5
    xdic=x_ds-xc
    ydic=y_ds-yc
    zdic=z_ds-zc
    vxdic=vx_ds-vxc
    vydic=vy_ds-vyc
    vzdic=vz_ds-vzc
    rdic=(xdic**2+ydic**2+zdic**2)**0.5
    r_din=rdic[rdic.argsort()];  r_sn=ric[ric.argsort()]
    m_dn=m_ds[rdic.argsort()];   m_sn=m_ss[ric.argsort()]
    x_dn=x_ds[rdic.argsort()];   x_sn=x_ss[ric.argsort()]
    y_dn=y_ds[rdic.argsort()];   y_sn=y_ss[ric.argsort()]
    z_dn=z_ds[rdic.argsort()];   z_sn=z_ss[ric.argsort()]
    vx_dn=vx_ds[rdic.argsort()];   vx_sn=vx_ss[ric.argsort()]
    vy_dn=vy_ds[rdic.argsort()];   vy_sn=vy_ss[ric.argsort()]
    vz_dn=vz_ds[rdic.argsort()];   vz_sn=vz_ss[ric.argsort()]

    dist=((xc-xcold)**2+(yc-ycold)**2+(zc-zcold)**2)**0.5
    if dist<0.005:
        break
    if dist>0.005:
        mold=m
        xcold=xc
        ycold=yc
        zcold=zc
        corenold=corenn
        corenn=int(0.98*corenold)
        indls=corenn
        n_s=indls
        mrlc=r_sn[corenn]
        for mnm in xrange(n_d):
            if r_din[mnm]<=mrlc:
                mnm+=1
            elif r_din[mnm]>mrlc:
                indl=mnm
                break
        n_s=indls
        n_d=indl
        myiter+=1
        print 'iteration',myiter
        x_ss=numpy.empty(indls)
        y_ss=numpy.empty(indls)
        z_ss=numpy.empty(indls)
        vx_ss=numpy.empty(indls)
        vy_ss=numpy.empty(indls)
        vz_ss=numpy.empty(indls)
        m_ss=numpy.empty(indls)
        x_ds=numpy.empty(indl)
        y_ds=numpy.empty(indl)
        z_ds=numpy.empty(indl)
        vx_ds=numpy.empty(indl)
        vy_ds=numpy.empty(indl)
        vz_ds=numpy.empty(indl)
        m_ds=numpy.empty(indl)
        ric=numpy.empty(indls);
        r_din=numpy.empty(indl); 
        for fs in xrange(indls):
            x_ss[fs]=x_sn[fs]
            y_ss[fs]=y_sn[fs]
            z_ss[fs]=z_sn[fs]
            vx_ss[fs]=vx_sn[fs]
            vy_ss[fs]=vy_sn[fs]
            vz_ss[fs]=vz_sn[fs]
            m_ss[fs]=m_sn[fs]
        for fn in xrange(indl):
            x_ds[fn]=x_dn[fn]
            y_ds[fn]=y_dn[fn]
            z_ds[fn]=z_dn[fn] 
            vx_ds[fn]=vx_dn[fn]
            vy_ds[fn]=vy_dn[fn]
            vz_ds[fn]=vz_dn[fn]
            m_ds[fn]=m_dn[fn]
        if myiter==150:
            print 'Warning! BAD ITERATION=',myiter
            raise StopIteration #MyError(myiter)
            break

################################################################################
#  Line of sight (c=center)
################################################################################
dlosc=((xc-xsun)**2+(yc-ysun)**2+(zc-zsun)**2)**0.5
xveclosc=(xc-xsun)/dlosc
yveclosc=(yc-ysun)/dlosc
zveclosc=(zc-zsun)/dlosc
                                         
bins=12
bind=100
n=numpy.zeros(bins)
corrected=numpy.empty(bind)
sumden_s=numpy.zeros(bins);                sumden_d=numpy.zeros(bind);
lrm=numpy.empty(bind);                     lrmd=numpy.empty(bind);
surf=numpy.empty(bins);                    dend=numpy.empty(bind);

# Relative and projected distances #
####################################
xic=x_s-xc
yic=y_s-yc
zic=z_s-zc
ric=(xic**2+yic**2+zic**2)**0.5
rc=(xc**2+yc**2+zc**2)**0.5

com=open("com_file.dat",'w')
print >>com,xc,yc,zc,vxc,vyc,vzc,dlosc,xveclosc,yveclosc,zveclosc
com.close

#print x_d,y_d,z_d
xdic=x_d-xc
ydic=y_d-yc
zdic=z_d-zc
rdic=(xdic**2+ydic**2+zdic**2)**0.5
r_dd=rdic[rdic.argsort()]
m_dd=m_di[rdic.argsort()]

r_si=(ric**2-((xic*xveclosc)+(yic*yveclosc)+(zic*zveclosc))**2)**0.5
r_di=rdic

########################
# Sorting the arrays:  #
########################
r_s=r_si[r_si.argsort()]
m_s=m_si[r_si.argsort()]
r_d=r_di[r_di.argsort()]
m_d=m_di[r_di.argsort()]
r3ds=ric[r_di.argsort()]

maxx=max(r_s) 
minn=min(r_s)
maxxd=max(r_d) 
minnd=min(r_d)
lrmax=math.log10(maxx)
lrmin=math.log10(minn)
dr=(lrmax-lrmin)/float(bind)
lrmaxd=math.log10(maxxd)
lrmind=math.log10(minnd)
drd=(lrmaxd-lrmind)/float(bind)


#Md inside rhalo
for ll in xrange(n_dold):
    if r_d[ll]<=scllngdm:
        nmd1dm=ll
    else:
        break
mmd1dm=m_d[0]*nmd1dm


#Ms inside r
for ll3 in xrange(n_sold):
    if r3ds[ll3]<=scllngdm:
        nms1dm=ll3
    else:
        break
mms1dm=m_s[0]*nms1dm
Msfinal=mms1dm*10.**5
Mdmfinal=mmd1dm*10.**5
totmfinal=(mms1dm+mmd1dm)*10.**5
mydensinrh=((mmd1dm+mms1dm)/((4.0*math.pi/3.0)*(scllngdm**3)))*1.989e35/(3.085e19)**3.0

##################################################################
# Calculating bin areas and the mean densities for each area     #
##################################################################
msun=1.989*10**30
kpc=3.085*10**19
arcmin=9.155397604e17/kpc

r_obs=numpy.empty(bins)

r_obs=[math.log10(0.0857038),math.log10(0.114725),math.log10(0.151126),math.log10(0.174134),math.log10(0.197232),math.log10(0.220339),math.log10(0.249179),math.log10(0.284005),math.log10(0.330016),math.log10(0.391652),math.log10(0.504713),math.log10(1.48085)] #rmax of the bins from matt's 2011data
r_ls=numpy.array(r_obs)
r_ls=10**r_ls

n=numpy.zeros(bins)
i=0
r_obd=numpy.empty(bind)
r_ld=numpy.array(r_obd)

rlbb=-2.0
for j2 in xrange(bind):
    r_ld[j2]=10**(rlbb)+(10**(rlbb))*j2*1.5

o=0
for k in xrange(bins):
    if k==0:
        n[k]=0
        print r_ls,r_s
        surf[k]=math.pi*((r_ls[k])**2.-(0)**2.)
        while 0<r_s[o]<=r_ls[k]: 
            sumden_s[k]+=m_s[o] 
            o+=1
            n[k]+=1
        print 'o,k and first',o,k,n[k] 
    elif k>=1:
        n[k]=0
        surf[k]=math.pi*((r_ls[k])**2.-(r_ls[k-1])**2.)
        print 'r_ls[k],k',r_ls[k],k            
        while r_ls[k-1]<r_s[o]<=r_ls[k]: 
            sumden_s[k]+=m_s[o] 
            o+=1
            n[k]+=1
            if o==n_sold:
                break
        print 'sumden_rest',sumden_s[k],n[k]
        if o==n_sold:
            break

g=0
nl=numpy.zeros(bins)
for g in xrange(bins):
    print 'n(g) inside',n[g],g
    if g==0:
      nmed=int(n[g]/2)
      r_obs[g]=r_s[nmed]
      nl[0]=n[0]
      print 'n(g[0]) inside',n[g],g,nmed,nl[g]
    elif g>=1:
        nl[g]=n[g]+nl[g-1]
        nmed=int(nl[g-1]+n[g]/2)
        r_obs[g]=r_s[nmed]
        print 'n nl[g] nl[g-1] nl[g-1]+n[g]/2 nmed r',n[g],nl[g],nl[g-1],nl[g-1]+n[g]/2,nmed,r_obs[g],g
        if nl[g]==n_sold:
            break
lrobs=numpy.log10(r_obs)

den_s=sumden_s/surf
dsigma=den_s/numpy.sqrt(n)
bulge=open(output1,'w')
print 'den_s',den_s
for l in xrange(bins):
    if den_s[l]>0:
        print >>bulge,den_s[l],dsigma[l],r_obs[l],lrobs[l],Msfinal,totmfinal
bulge.close()

#printed the median r of the bin which i used to calculate the density for the bin limited by the r_lim defined by the outer observed radii.
e=0
i=0

n=numpy.zeros(bind)
for e in xrange(bind):
    if e==0:
        n[e]=0
        print 'r_ld inside',r_ld
        dend[e]=4./3.*math.pi*((r_ld[e])**3.-(0)**3.)
        while  0<=r_d[i]<r_ld[e]: 
            sumden_d[e]+=m_d[i]
            i+=1
            n[e]+=1
    elif e>=1:
        n[e]=0
        dend[e]=4./3.*math.pi*((r_ld[e])**3.-(r_ld[e-1])**3.)
        while r_ld[e-1]<=r_d[i]<r_ld[e]: 
            sumden_d[e]+=m_d[i]  
            i+=1
            n[e]+=1
            if i==n_dold:
                break
        if i==n_dold:
            break

den_d=sumden_d/dend

g2=0
nl2=numpy.zeros(bind)
for g2 in xrange(bind):
    nl2[0]=0
    if g2==0:
      if n[g2]%2==0:
          nmed2=n[g2]/2
          r_obd[g2] = 0.5*(r_d[nmed2-1]+r_d[nmed2])
      else:
          nmed2=(n[g2]+1)/2
          r_obd[g2] = r_d[nmed2-1]
      r_obd[g2]=r_d[nmed2]
      nl2[0]=n[0]
    elif g2>=1:
        nl2[g2]=n[g2]+nl2[g2-1]
        if n[g2]%2==0:
            nmed2=nl2[g2-1]+n[g2]/2
            r_obd[g2]=0.5*(r_d[nmed2-1]+r_d[nmed2])
        else:
            nmed2=nl2[g2-1]+n[g2]/2
            r_obd[g2]=r_d[nmed2-1]
        if nl2[g2]==n_dold:
            break
lrobd=numpy.log10(r_obd)

dsumden=den_d/numpy.sqrt(n)

for p in xrange(bind):
    corrected[p]=den_d[p]-dsumden[p]

halo=open(output2,'w')
for p in xrange(bind):
     if den_d[p]>0:
         print >>halo,den_d[p],r_obd[p],lrobd[p],dsumden[p],Mdmfinal
halo.close()

######################################
# Including the observational data     #
######################################

OBSdata=open('obsbnoisy.dat','r') 
rth1_list=[]
rth2_list=[]
sbth1_list=[]
dsb1_list=[]

def file2listtheo(filename,listsb1,listdsb1,listrth1,listrth2):
    for line in filename:
        d=line.split()
        listsb1+=[float(d[0])]
        listdsb1+=[float(d[1])]
        listrth1+=[float(d[2])]
        listrth2+=[float(d[3])]
       
file2listtheo(OBSdata,sbth1_list,dsb1_list,rth1_list,rth2_list)  

rth1=numpy.array(rth1_list)
rth2=numpy.array(rth2_list)
sbth1=numpy.array(sbth1_list)
dsbth1=numpy.array(dsb1_list)

OBSdata.close()

lrobs=numpy.log10(r_obs)
lrobd=numpy.log10(r_obd)

####################################
#   PLOT surface brightness        #
####################################

outf1=open(outfig1,'w')
figure()
xlabel('r[kpc]')
ylabel(r'\Sigma')
loglog()
p_MWs=errorbar(r_ls,den_s, yerr=dsigma,xerr=None,
               ecolor=None, elinewidth=None, capsize=3,
               barsabove=True, lolims=False, uplims=False,
               xlolims=False, xuplims=False,label='Bulge')
p_obs=errorbar(rth1,sbth1, yerr=dsbth1,xerr=None,
                ecolor=None, elinewidth=None, capsize=3,
                barsabove=True, lolims=False, uplims=False,
                xlolims=False, xuplims=False,label='Stars-obs')
legend()
xlim(0.07,2.0)
savefig(outf1)
outf1.close()

outf2=open(outfig2,'w')
figure()
xlabel('r [kpc]')
p_MWs5=errorbar(r_obd,den_d, yerr=dsumden,xerr=None,
                ecolor=None, elinewidth=None, capsize=3,
                barsabove=True, lolims=False, uplims=False,
                xlolims=False, xuplims=False)
legend()
savefig(outf2)
outf2.close



