###################################################################
# PLOTTING THE SURFACE BRIGHTNESS PROFILE                         #
#
# using the tipsy input
# RUN with python denGITHUB.py  initial.bin bulge.dat halo.dat bulge.png halo.png
###########
#TIPSY FILES:
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
scllngdm=0.82 ##the radius of the outermost dispersion bin

######################################
# Assigning sun's coordinates        #                     
######################################
xsun=8.5;ysun=0.0 ;zsun=0.0; vxsun=0.0; vysun=200.0; vzsun=0.0; sol_m=1.0;#e5*1.989e30 
print 'msun is for pkdgrav would be sol_m=2.325e9*1.989e30 for galactics'

###################################################
# Reading data from tipsy file and making arrays  #
###################################################
h   = TipsyHeader()                      #create a location for the header and all the rest
inp = ifTipsy(filename, "standard")      #opens the file
inp.read(h)                              # read the header

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
inp.seekg(tipsypos(tipsypos.dark, 0))#jump directly to the first dm particle(0)
for i in xrange(n_d): #(h.h_nDark) interval from 0(in dark) to (h.h_nDark-1)
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
# 1- INITIAL VALUES
####################
corenold=99999
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

mold=mtotoldcore #sum(m_di)+sum(m_si)
vxcold=vxcoldcore[coren-1]/coren #(m_di[0]*sum(vx_d)+sum(vx_s)*m_si[0])/mold
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

# 4- CHOOSING THE LIMITS FOR THE FIRST ITERATION
#################################################
#rmas=max(r_snold)
rmas=r_snold[corenold]
print 'r_snold[0],r_snold[1],r_snold[2],r_snold[3],r_snold[4],rmas,n_s',r_snold[0],r_snold[1],r_snold[2],r_snold[3],r_snold[4],rmas,n_s
corenn=0.98*corenold
indls=corenn

indl=corenn
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

# 5 - NEW EMPTY ARRAYS FOR FEWER NUMBER OF PARTICLES
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

# 6 - ITERATION FOR THE CENTRE OF MASS
######################################
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
    # Distance criteria for aceptance
    dist=((xc-xcold)**2+(yc-ycold)**2+(zc-zcold)**2)**0.5
    if dist<0.002:
        break
    if dist>0.002:
        mold=m
        xcold=xc
        ycold=yc
        zcold=zc
        corenold=corenn
        corenn=0.98*corenold
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
        x_ss=numpy.empty(indls);        y_ss=numpy.empty(indls)
        z_ss=numpy.empty(indls);        vx_ss=numpy.empty(indls)
        vy_ss=numpy.empty(indls);        vz_ss=numpy.empty(indls)
        m_ss=numpy.empty(indls);        x_ds=numpy.empty(indl)
        y_ds=numpy.empty(indl);        z_ds=numpy.empty(indl);
        vx_ds=numpy.empty(indl);        vy_ds=numpy.empty(indl);
        vz_ds=numpy.empty(indl);        m_ds=numpy.empty(indl)
        ric=numpy.empty(indls);        r_din=numpy.empty(indl); 
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

##################################################################################
#  Calculating the distance and the unit vectors of the line of sight (c=center) #
##################################################################################
dlosc=((xc-xsun)**2+(yc-ysun)**2+(zc-zsun)**2)**0.5
xveclosc=(xc-xsun)/dlosc
yveclosc=(yc-ysun)/dlosc
zveclosc=(zc-zsun)/dlosc
                                         
bins=11
bind=100
n=numpy.zeros(bins)
corrected=numpy.empty(bind)
sumden_s=numpy.zeros(bins);                sumden_d=numpy.zeros(bind);
lrm=numpy.empty(bind);                     lrmd=numpy.empty(bind);
surf=numpy.empty(bins);                    dend=numpy.empty(bind);
#####################################################################################
#  Relative distance and velocity of each star from the center of Carina and        #  
#####################################################################################
# ic=star to center #
#####################
xic=x_s-xc
yic=y_s-yc
zic=z_s-zc
ric=(xic**2+yic**2+zic**2)**0.5
rc=(xc**2+yc**2+zc**2)**0.5

com=open("com_file.dat",'w')
print >>com,xc,yc,zc,vxc,vyc,vzc
com.close

# dic=star to center #
#####################
#print x_d,y_d,z_d
xdic=x_d-xc
ydic=y_d-yc
zdic=z_d-zc
rdic=(xdic**2+ydic**2+zdic**2)**0.5
r_dd=rdic[rdic.argsort()]
m_dd=m_di[rdic.argsort()]

print 'RDD ',r_dd

#############################
#  Projected distance(star) #
#############################
#r_si=(ric**2-((xic*xveclosc)**2+(yic*yveclosc)**2+(zic*zveclosc)**2))**0.5
r_si=(ric**2-((xic*xveclosc)+(yic*yveclosc)+(zic*zveclosc))**2)**0.5

#############################
#  Projected distance(dm)   #
#############################
r_di=rdic
#r_di=numpy.sqrt(xdic**2+ydic**2)
########################
# Sorting the arrays:  #
########################
r_s=r_si[r_si.argsort()]
m_s=m_si[r_si.argsort()]
r_d=r_di[r_di.argsort()]
m_d=m_di[r_di.argsort()]
r3ds=ric[r_di.argsort()]
####################################
# Calculating bin sizes            #
####################################
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


########################
# Calculating mass loss#
########################
#mloss=open(output3,'w')

####Md dark matter mass inside rhalo
for ll in xrange(n_dold):
    if r_d[ll]<=scllngdm:
        nmd1dm=ll
    else:
        break
mmd1dm=m_d[0]*nmd1dm


####Ms stellar mass inside rhalo    
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


#mloss=open(output3,'w')
#print >>mloss,mydensinrh,totmfinal
#mloss.close()


##################################################################
# Calculating bin areas and the mean densities for each area     #
##################################################################
#rm is the average r of each bin and lrm is log(rm)
msun=1.989*10**30#kg
kpc=3.085*10**19#m
arcmin=9.155397604e17/kpc #0.02967kpc

# defined the bin radii from Car_disp_obs.dat by chosing the limits half-way through the bins: 
# r_obs=[math.log10(1.449000*arcmin),math.log10(4.347000*arcmin),math.log10(7.245000*arcmin),math.log10(10.143000*arcmin),math.log10(14.490000*arcmin),math.log10(20.286000*arcmin),math.log10(26.082000*arcmin),math.log10(31.878000*arcmin)]
# lr_ave=numpy.array(r_obs)
# lr_ave=r_obs
# rl=[2.898*arcmin,5.796*arcmin,8.694*arcmin,12.3165*arcmin,17.388*arcmin,23.184*arcmin,28.98*arcmin,maxx]
# r_l=numpy.array(rl)

r_obs=numpy.empty(bins)
r_ls=numpy.empty(bins)

rlaa=-0.8286
for j in xrange(bins):
    r_ls[j]=10**(rlaa)+(10**(rlaa))*j
print 'here',r_ls[j]

n=numpy.zeros(bins)
i=0
r_obd=numpy.empty(bind)
r_ld=numpy.empty(bind)

rlbb=-2.0

for j2 in xrange(bind):
    r_ld[j2]=10**(rlbb)+(10**(rlbb))*j2*1.5

# SURFACE DENSITY PER STELLAR BIN (BINS)
o=0
for k in xrange(bins):
    if k==0:
        n[k]=0
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
for l in xrange(bins):
    if den_s[l]>0:
        print >>bulge,den_s[l],dsigma[l],r_obs[l],lrobs[l],Msfinal
bulge.close()


# DENSITY PER DARK MATTER BIN (BIND)
e=0
i=0

print 'Mass = ',m_d[30]

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

print 'totalmass',sum(sumden_d)+sum(sumden_s)
print 'if you use too many bins better to put a max limit for the radius of dm particles to take into account'

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

lrobs=numpy.log10(r_obs)
lrobd=numpy.log10(r_obd)


#############################
# Theoretical profiles      #
#############################
# NOT USED FOR THE CURRENT MODELS

#PLUMMER Bulge
#Rpb=1.0
#I_0=1.59202e+2 #given as 1.59202e+07 M_sun/kpc^3 un the Bulge.tab
#I_0=4.0/3.0*1.59202e+2 #given as 1.59202e+07 M_sun/kpc^3 un the Bulge.tab
#I=(I_0*Rpb**5)/(Rpb**2+r_obs**2)**2 


#PLUMMER Halo
#Rpd=1.0
#rhoh=1.59202e+3 #given as 1.59202e+08 M_sun/kpc^3 in the Halo.tab
#rhot=rhoh*Rpd**5/(Rpd**2+r_obd**2)**2.5

#theo=open('theoplumhalo.dat','w')
#q=0
#for q in xrange(bind):
 #   print >>theo,rhot[q],r_obd[q],lrobd[q]


#theoplum=open('theoplumbulge.dat','w')
#q=0
#for q in xrange(bins):
#    print >>theoplum,I[q],r_obs[q],lrobs[q]

######################################
# Including the observational data     #
######################################

OBSdata=open('WalkerCarina_SBP.dat','r')

rth1_list=[]
rth2_list=[]
sbth1_list=[]
dsb1_list=[]

def file2listtheo(filename,listrth1,listsb1,listdsb1):
    for line in filename:
        d=line.split()
        listrth1+=[float(d[0])]
        listsb1+=[float(d[1])]
        listdsb1+=[float(d[2])]

file2listtheo(OBSdata,rth1_list,sbth1_list,dsb1_list)  

rth1=numpy.array(rth1_list)
sbth1=numpy.array(sbth1_list)/1.0e5
dsbth1=numpy.array(dsb1_list)/1.0e5
OBSdata.close()

lrobs=numpy.log10(r_obs)
lrobd=numpy.log10(r_obd)


####################################
#   PLOT surface brightness        #
####################################

softening=r_obs-r_obs+0.05

from matplotlib import pylab
from pylab import figure, plot, show,semilogx,xscale,yscale,errorbar, title,savefig, legend, xlabel, ylabel, close, xlim, ylim, semilogy,loglog,axis


outf1=open(outfig1,'w')
figure()
xscale("log");yscale("log")
xlabel('r (kpc)')
ylabel('density ($10^5$M$_{\odot}$kpc$^{-2}$)')
p_MWs=errorbar(r_obs,den_s, yerr=dsigma,xerr=None,
                ecolor=None, elinewidth=None, capsize=3,
                barsabove=True, lolims=False, uplims=False,
                xlolims=False, xuplims=False,label='BULGE-simulation')
p_obs=errorbar(r_obs,sbth1, yerr=dsbth1,xerr=None,
                ecolor=None, elinewidth=None, capsize=3,
                barsabove=True, lolims=False, uplims=False,
                xlolims=False, xuplims=False,label='Bulge-observations')
peps=loglog(softening,den_s,'g--')
ylim(0.001,1000.0)
legend()
savefig(outf1)
outf1.close

softening=r_obd-r_obd+0.03

outf2=open(outfig2,'w')
figure()
xlabel('r (kpc)')
ylabel('density (M$_{\odot}$kpc$^{-3}$)')
p_MWs5=errorbar(r_obd,den_d, yerr=dsumden,xerr=None,
                ecolor=None, elinewidth=None, capsize=3,
                barsabove=True, lolims=False, uplims=False,
                xlolims=False, xuplims=False)#,label='HALO data')
peps2=loglog(softening,den_d,'g--')#,label='Softening length')
legend()
savefig(outf2)
outf2.close
