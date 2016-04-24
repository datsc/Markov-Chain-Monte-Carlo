###################################################################
#  Calculates the centre of mass
###################################################################
import sys
from pytipsy import ifTipsy,TipsyHeader, TipsyDarkParticle,tipsypos, TipsyStarParticle
import pytipsy
import math
import numpy

###############################
# Takes the arguments        #
###############################
if len(sys.argv) <=1:
    print 'Not enough input files'
filename=sys.argv[1]  #the final output of pkdgrav 
fileout=sys.argv[2]

xsun=8.5;ysun=0.0 ;zsun=0.0; vxsun=0.0; vysun=200.0; vzsun=0.0; sol_m=1.0

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

corenold=99990
mtotoldcore=corenold*m_si
xcoldcore=numpy.empty(corenold)
ycoldcore=numpy.empty(corenold)
zcoldcore=numpy.empty(corenold)
vxcoldcore=numpy.empty(corenold)
vycoldcore=numpy.empty(corenold)
vzcoldcore=numpy.empty(corenold)

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
vxcold=vxcoldcore[corenold-1]/corenold
vycold=vycoldcore[corenold-1]/corenold
vzcold=vzcoldcore[corenold-1]/corenold
xcold=xcoldcore[corenold-1]/corenold 
ycold=ycoldcore[corenold-1]/corenold 
zcold=zcoldcore[corenold-1]/corenold 

# RELATIVE VALUES FOR INDIVIDUAL PARTICLES
##############################################
xicold=x_s-xcold; yicold=y_s-ycold; zicold=z_s-zcold
ricold=(xicold**2+yicold**2+zicold**2)**0.5
xdicold=x_d-xcold; ydicold=y_d-ycold; zdicold=z_d-zcold
rdicold=(xdicold**2+ydicold**2+zdicold**2)**0.5


r_dinold=rdicold[rdicold.argsort()];  r_snold=ricold[ricold.argsort()]
m_dnold=m_di[rdicold.argsort()];      m_snold=m_si[ricold.argsort()]
x_dnold=x_d[rdicold.argsort()];    x_snold=x_s[ricold.argsort()]
y_dnold=y_d[rdicold.argsort()];    y_snold=y_s[ricold.argsort()]
z_dnold=z_d[rdicold.argsort()];    z_snold=z_s[ricold.argsort()]
vx_dnold=vx_d[rdicold.argsort()];    vx_snold=vx_s[ricold.argsort()]
vy_dnold=vy_d[rdicold.argsort()];    vy_snold=vy_s[ricold.argsort()]
vz_dnold=vz_d[rdicold.argsort()];    vz_snold=vz_s[ricold.argsort()]


rmas=r_snold[corenold]
print 'r_snold[0],r_snold[1],r_snold[2],r_snold[3],r_snold[4],rmas,n_s',r_snold[0],r_snold[1],r_snold[2],r_snold[3],r_snold[4],rmas,n_s
corenn=int(0.98*corenold)
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

# NEW ARRAYS FOR FEWER NUMBER OF PARTICLES
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
    if dist<0.002:
        break
    if dist>0.002:
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

r_c=(xc**2+yc**2+zc**2)**0.5
com=open(fileout,'w')
print >>com,xc,yc,zc,vxc*0.655875368,vyc*0.655875368,vzc*0.655875368,r_c
com.close
