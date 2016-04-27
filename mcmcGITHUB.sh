#PBS -l nodes=2:ppn=16
#PBS -l walltime=5:00:00:00
#PBS -l pvmem=1gb
#PBS -m bea
#PBS -M uu2@le.ac.uk
#-o  /scratch/$JOB_NAME.o$JOB_ID
#-e  /scratch/$JOB_NAME.e$JOB_ID
#PBS -v PYTHONPATH=/home/hannilib/
#PBS -S /bin/bash
#======================================
# Set up environment
#======================================
module load openmpi
module load gcc
module load python/2.7.3
module load swig
#======================================
cd /scratch/C2C_ART1_2013/C2C3

export dirnm=`date +%d%m%Y`
SIMDIR=/scratch/dp005/dc-ural1/C2C_ART1_2013/C2C3
PKDGRAV=/home/dc-ural1/mpkdgravopmpinopottry_2012/pkdgrav.d.gnu
cd $SIMDIR
echo 'The simulation is running in ' $SIMDIR
rm -f newnemscr.sh initial

# C1 Generating random numbers

l=1 #should be 1 less than scripts f$bla number
m=`echo "$l+1"|bc -l` #This is the value going with the new Valuestokeep, Accepted and inputs of the next chain.
INPARFILE=$SIMDIR/myinitialnoexpparam$l.dat  #To this one I write the values of the last model that started
LASTACCFILE=$SIMDIR/myfinalacceptedparam$l.dat #To this one I write the values of the last model that was accpeted including teh chi_sq
if [ -e $INPARFILE ];then
    echo 'great';
    else
    echo 'MODIFY THE INPUT';
    break
fi

###############################
# Read the initial parameters #
###############################
#The original target values are: 550000.0 2.474 0.279 140000000.0 629.814 1.0 0.125 0.11
# C2 Reading data from myinitialnoexpparam.dat
Ms=`echo |tail -n 1 $LASTACCFILE|awk '{print $1}'`    #M_S
mb1=`echo |tail -n 1 $LASTACCFILE|awk '{print $2}'`   #M_S/222288
rst=`echo |tail -n 1 $LASTACCFILE|awk '{print $3}'`   #r_s
Mhalo=`echo |tail -n 1 $LASTACCFILE|awk '{print $4}'` #M_DM
mdm1=`echo |tail -n 1 $LASTACCFILE|awk '{print $5}'`  #M_DM/222288
r1=`echo |tail -n 1 $LASTACCFILE|awk '{print $6}'`    #r_dm        
mualpha=`echo |tail -n 1 $LASTACCFILE|awk '{print $7}'` #mu_alpha*cos(delta) [So need to multiply Piatek values with cos(delta)]
mudelta=`echo |tail -n 1 $LASTACCFILE|awk '{print $8}'` #mu_delta
chi2_0=`echo |tail -n 1 $LASTACCFILE|awk '{print $9}'` #chi2
chi2sb_0=`echo |tail -n 1 $LASTACCFILE|awk '{print $10}'` #chi2_sb
chi2sig_0=`echo |tail -n 1 $LASTACCFILE|awk '{print $11}'` #chi2_sig

Msn=`echo |tail -n 1 $INPARFILE|awk '{print $1}'`    #M_S
mb1n=`echo |tail -n 1 $INPARFILE|awk '{print $2}'`   #M_S/222288
rstn=`echo |tail -n 1 $INPARFILE|awk '{print $3}'`   #r_s
Mhalon=`echo |tail -n 1 $INPARFILE|awk '{print $4}'` #M_DM
mdm1n=`echo |tail -n 1 $INPARFILE|awk '{print $5}'`  #M_DM/222288
r1n=`echo |tail -n 1 $INPARFILE|awk '{print $6}'`    #r_dm        
mualphan=`echo |tail -n 1 $INPARFILE|awk '{print $7}'` #mu_alpha*cos(delta) [So need to multiply Piatek values with cos(delta)]
mudeltan=`echo |tail -n 1 $INPARFILE|awk '{print $8}'` #mu_delta
o=`echo |tail -n 1 $INPARFILE|awk '{print $9}'` #model number

i=`echo "$o+1"|bc -l`
echo $o "and" $i

echo 'At the begining of rest2 the last used initial and bulge files are:'>>Checkingthefiles.txt
ls -lt bulge* >> Checkingfiles.txt
ls -lt BULGES/bulge* >> Checkingfiles.txt
echo 'The new file that I will write will be called the same?:bulgeo'$o>>Checkingthefiles.txt
rm -f newnemscr.sh initialeq$o.* initial$o.* bulge*$o comeq* B$o.dat H$o.dat B$o.ini H$o.ini halo*$o orbx$o.dat orbxpast$o.dat simb$o.dat simh$o.dat simdisp$o.dat simpa$o.dat simgrad$o.dat
ls -lt bulge* >> Checkingfiles.txt
echo 'The above is after I clean this directory'>>Checkingthefiles.txt

#C3 Setting the initial parameters as the old parameters
Msold=`echo "$Ms"|bc -l`
mb1old=`echo "$mb1"|bc -l`
rstold=`echo "$rst"|bc -l`
Mhaloold=`echo "$Mhalo"|bc -l`
mdm1old=`echo "$mdm1"|bc -l`
r1old=`echo "$r1"|bc -l`
mualphaold=`echo "$mualpha"|bc -l`
mudeltaold=`echo "$mudelta"|bc -l`
chi2old=`echo "$chi2_0"|bc -l`

#C4 Setting the first new parameters as the old parameters
Msnew=`echo "$Msn"|bc -l`
mb1new=`echo "$mb1n"|bc -l`
rstnew=`echo "$rstn"|bc -l` 
Mhalonew=`echo "$Mhalon"|bc -l`
mdm1new=`echo "$mdm1n"|bc -l`
r1new=`echo "$r1n"|bc -l`
mualphanew=`echo "$mualphan"|bc -l`
mudeltanew=`echo "$mudeltan"|bc -l`

#C5 Output to the values to keep
echo 'Initial values Msold mb1old rstold  Mhaloold mdm1old r1old mualhaold mudeltaold' $Msold $mb1old $rstold $Mhalo1old $mdm1old $r1old $mualphaold $mudeltaold >> Valuestokeep$m.txt
echo 'Initial values Ms mb1 rst Mhalo mdm1 r1 mualpha mudelta' $Msnew $mb1new $rstnew $Mhalonew $mdm1new $r1new $mualphanew $mudeltanew>> Valuestokeep$m.txt
echo 'o=' $o >> Valuestokeep$m.txt
#C6 output to the .o file
echo  'Initial values Msold mb1old rstold  Mhaloold mdm1old r1old mualhaold mudeltaold' $Msold $mb1old $rstold $Mhalo1old $mdm1old $r1old $mualphaold $mudeltaold
echo 'Initial values Ms mb1 rst Mhalo mdm1 r1 mualpha mudelta' $Msnew $mb1new $rstnew $Mhalonew $mdm1new $r1new $mualphanew $mudeltanew

#FINISHED 1##

##########################
# Calling tcsh for Nemo  #
##########################
# NEMO DETAILS:
#C7 mkhalo with 200000 particles
#C8  stellar component from fit inner=0.37 outer=7 eta=2.5 r= r_st  (0.3 for initial fit to data)
#C8NEWUgur stellar component from fit inner=0.515 outer=4.45 eta=3.745 r=0.3 
#C8NEWMark stellar component from fit inner=0.8711 outer=3.472 eta=6.046 r=0.3 
#C9 potential on the stellar component Hernquist r1 Mdm
#C8  dark matter component from fit inner=1 outer=4 eta=1 r= r_st 
#C9 potential on the dark matter stellar fit in, out, eta, rb Mb

##OLD NEMO PARAMS FOR SB
#mkhalo B0.ini 100000 inner=0.37 outer=7.0 r_s='$rstold 'M='$Msold 'eta=2.5 b=0.0 r_a=0.0 accpars=0.
#0,'$r1old','$mdm1old',1.0,4.0,1.0,0.0,0.0 accname=Halo tabfile=B1.tab r_t=10 WD_units=t

echo  '#!/bin/tcsh
if(! $?NEMO) then
setenv NEMO /data/dp005/dc-ural1/nemo_cvs
module load gcc
source $NEMO/nemo_start
endif
setenv FALCON $NEMO/usr/dehnen/falcON
source ${FALCON}/falcON_start
cd ' $SIMDIR '
mkhalo B'$o'.ini 100000 inner=0.515 outer=4.45 r_s='$rstnew 'M='$Msnew 'eta=3.745 b=0.0 r_a=0.0 accpars=0.0,'$r1new','$mdm1new',1.0,4.0,1.0,0.0,0.0 accname=Halo tabfile=B1.tab r_t=10 WD_units=t
mkhalo H'$o'.ini 100000 inner=1.0 outer=4.0 r_s='$r1new' M='$Mhalonew' eta=1.0 b=0.0 r_a=0.0 accpars=0.0,'$rstnew','$mb1new',0.515,4.45,3.745,0.0,0.0 accname=Halo tabfile=H1.tab r_t=10 WD_units=t 
s2a B'$o'.ini B'$o'.dat
s2a H'$o'.ini H'$o'.dat
source $FALCON/falcON_end
source $NEMO/nemo_end
' >> newnemscr.sh 

#C10 Running NEMO (newnemscr.sh) to generate the initial galaxy.
chmod +x newnemscr.sh
./newnemscr.sh

#C11 awk to write the bulgeo and haloo which will go to the unit conversion code (it reads the 7 columned part of the NEMO generated files.)

# C12 Old conversion factor for v_nemo-->v_pkd was 1.490932795
# C13 New version comes from The below conversion
# NEMO units     M=222288msun G=1 r=kpc t=gyrs v=r/t=977.792m/s
# PKDRAV units   M=10^5Msun G=1 r=kpc v=sqrt(GM/r)=655.875368m/s T=r/v=3.085677e19/v=1.49081988e9yrs
# Physical units M1=10^5msun m2=msun v=km/s r=kpc t=r/v=0.977792e9yrs
# Physical to PKDGRAV conversion vphys/vpkd=1.5246 tphys/tpkd=0.655875343
# Nemo to PKDGRAV conversion vnem/vpkd=1.49081988

awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' B$o.dat > bulgeo
awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' H$o.dat > haloo

#C14 Calculated the circular orbit (output is 290,73,0,35,-139,0 in km/s for pkdgrav)
echo >> orfcirc, "300.0 0.0 0.0 0.0 0.1475 0.0 0.0001 5000.0" 
python Ugur_orbit1.py orfcirc oroutcirc
echo 'initial galaxy generated'

echo '300.0 0.0 0.0 0.0 0.1475 0.0 0.00001 5000.0' > orfcirc

python Ugur_orbit1.py orfcirc oroutcirc
#C14  Bulgeo in bulge out for the circular orbit, the new coordinates using the new c.o.m. and converting the velocity from nemo(~km/s) to pkdgrav(~0.65km/s) 

python myunits5.py 200000 100000 100000 oroutcirc
#C15 Halo and bulge converted and centred on the far orbit.

./a2b 1.0 1.0
#C16 Printing the first binary file initial.bin from halo and bulge and tagging the halo and bulge files in the chains

mv haloo halooin$o
mv bulgeo bulgeoin$o
mv halo haloin$o
mv bulge bulgein$o

mpirun  /home/dc-ural1/mpkdgravopmpinopottry_2012/pkdgrav.d.gnu simeq.param
mv initial.bin initialeq$o.bin
#C17 the first 500Myrs run for the relaxation is made with the right version of the time step which is step=0.067077
#step=0.067077*1.490819e9=10^8years.

python centreofmass.py initialeq$o.bin comeq$o.dat
#C18 the centre of mass is calculated for the d=300kpc equilibrium run.. The output units are converte to kpc, and km/s instead of pkdgrav units... the vl printed out in [km/s] 

python checkinginitialbin2.py initialeq.00005 initialeq$o.00005.dat
#C19 outputs the binary file to an ascii file with m,x,y,z,vx,vy,yz,eps

head -n 100000 initialeq$o.00005.dat >> haloo
tail -n 100000 initialeq$o.00005.dat >> bulgeo
#C20 dividing the ascii file from the equilibrium run into the new haloo and bulgeo for this run

mv initialeq.00005 initialeq$o.00005
python centreofmass.py initialeq$o.00005 comeq5.dat
#C21  the centre of mass is calculated for the 6Gyrs orbit run

./adconvugur6gyrs 100.4029 -50.9661 100.9 $mualphaold $mudeltaold 223.1
#C22 this one outputs all the velocities[kpc/Myrs] with opposite signs which can be directly integrated backwards.
# The input to the orbit code comes as kpc and kpc/myrs and the output is kpc and -km/s so inverted again, ready to be integrated forward.

mv convvel.dat orbx$o.dat 
python Ugur_orbit1.py orbx$o.dat orbxpast$o.dat
#C23 From the final equilibrium file's haloo an bulgeo, it subtracts the c.o.m and it shifts it to orbit written in pkdgrav units in the halo and bulge file.[x=x-x_com_eq+x_orbit and vx=vx-vx_com_eq*1.5246+vx_orbit*1.5246]

python myunits6goingtozero.py 200000 100000 100000 orbxpast$o.dat  comeq5.dat
#Reads in the new (non-equil, real orbit) haloo, bulgeo and orbital parameters from 6gyrs ago(orbxpast0.dat) and from the equilibrium file(comeq5.dat)... Sorts particles by distances to the centre and shifts them to orbit in halo and bulge in pkdgrav units: so v[0.655km/s]

./a2b 1.0 1.0
#C24 Printing the binary file initial.bin for the first simulation.
mv haloo halooorb$o
mv halo halorb$o
mv bulgeo bulgeoorb$o
mv bulge bulgeorb$o
#C25 The first simulation on the orbit.
mpirun /home/dc-ural1/mpkdgravopmpinopottry_2012/pkdgrav.d.gnu sim1.param 
echo 'o'$o initialeq$o.00005
mv initial.bin initial$o.bin
mv final.00060 final$o.00060
#C26 Renaming the binary files of the first simulation.
echo 'Im haaaaapppy'
python dentry4coredart_2012.py final$o.00060 simb$o.dat simh$o.dat simb$o.png simh$o.png
# C27 Calculates the SB and mass of the first simulation.Need to change the scllngdm from 828pc to 1.5kpc as 0.82 is the average r of the last bin.

python dispGitHUB.p final$o.00060 simdisp$o.dat simdisp$o.png
#C28 calculates the dispersion profile of the first simulation... JUST MAKE SURE THAT you use recalculated gensimb and gensimdisp as these should be calculated for the right radii. Output is in kpc and km/s.

python checkinginitialbin2.py final$o.00060 final$o.00060.dat
tail -n 100000 final$o.00060.dat >> final$o.00060_stars.dat 
./posangle6_2013 100000 final$o.00060_stars.dat simpa$o.dat com_file.dat
python meanvelnew.py simpa$o.dat simgrad$o.dat 100000 
#C28_2 These codes write stars file and then uses position angle to find 2d positions on the sky and finally calculates the mean velocities along the major axis

python stattestGITHUBart.py obsdisp146noisy180.dat obsb146noisy141.dat simdisp$o.dat simb$o.dat simgrad$o.dat simchi$o.dat
#C29 Calculating the profiles and chi-sq or the first simulation without including the inner most data point



#C30 Chosing the step sizes for the parameters of the MCMC
#[ln(10^10)-ln(10^6)]/10=0.92
#[ln(0.744)-ln(0.0744)]/10=0.23
dalpha=`echo "0.078"|bc -l` 
ddelta=`echo "0.066"|bc -l` 
dMs=`echo "215000.0"|bc -l` 
dr1=`echo "0.23"|bc -l` 
dmhalo=`echo "0.92"|bc -l` 
drst=`echo "0.23"|bc -l` 


####################
#  LOOP            #
####################
#C31 Setting the chi^2 to zero so that the first model is always accepted.
finished=0

  seed=`echo "-235740136"|bc -l`
  pold=`echo "-0.1"|bc -l`

  echo 'BURASI NOKTA 1'
  nacc=`echo "0.0"|bc -l`
  nrej=`echo "0.0"|bc -l`
  nbatch=0
#C31 reading the chi^2 from the simchi file
chi2_0=`echo| awk '{print $1}' simchi$o.dat`
chi2sb_0=`echo| awk '{print $2}' simchi$o.dat`
chi2sig_0=`echo| awk '{print $3}' simchi$o.dat`
###echo 1 |awk '{srand('$seed');for(i=1;i<=5000;i++){printf("%lf\n",rand())}}' >>random.dat
#C32 printing out the random.dat that has the numbers to be compared to the likelihood ratios in the accept/reject point. 
ran_count=`echo "$o+1"|bc -l`
ranal_count=`echo "$o+1"|bc -l`
randel_count=`echo "$o+1"|bc -l`
ranms_count=`echo "$o+1"|bc -l`
ranr1_count=`echo "$o+1"|bc -l`
ranmdm_count=`echo "$o+1"|bc -l`
ranrst_count=`echo "$o+1"|bc -l`
chi2old=$chi2old;	
chi2_0=`echo| awk '{print $1}' simchi$o.dat`
echo 'First simulation, therefore, chi2old and chi2_0 are' $chi2old $chi2_0 

echo "IMPORTANT CHECKPOINT 1" $i "and" $o "and" $ran_count
ls bulge* final* initial*

################
#  COMPARISON  #
################
while [ $finished -lt 1 ]
    do    
echo "IMPORTANT CHECKPOINT 2" $i $o 
ls bulge* final* initial*
mv *png FIGURES/.
mv comeq5*dat CENTREofMASSES/.
mv orbx*dat ORBITS/.
mv bulge* BULGES/.
mv halo* HALOS/.
mv B* BULGES/.
mv H* HALOS/.
mv oldnem* GENSCRIPTS/.
    if [ $i> 1 ]; then
	#chi2old=$chi2_0;
	j=`echo "$i-1"|bc -l`
	chi2_0=`echo| awk '{print $1}' simchi$j.dat` 
	chi2sb_0=`echo| awk '{print $2}' simchi$j.dat` 
	chi2sig_0=`echo| awk '{print $3}' simchi$j.dat` 
	#C33 reading the new chi2 values for the simulations
	#C34 The part below is for when a simulation crashes and I do not have a chi^2 value. Then I just set chi^2 to avery high number and go to the next step.
	if [ -s  $chi2_0 ]  ; then
	    chi2_0=`echo "3000.0"|bc -l `
	    chi2sb_0=`echo "3000.0"|bc -l `
	    chi2sig_0=`echo "3000.0"|bc -l `
	    echo 'i is larger than 1, the chi^2 could not be calculated and was put back to 3000 '>> Valuestokeep$m.txt 
	    mv simb$j.dat CRASHEDMODELS/.
	    mv simh$j.dat CRASHEDMODELS/.
	    mv simdisp$j.dat CRASHEDMODELS/.
	    mv simchi$j.dat CRASHEDMODELS/.
	else
#C35 Setting the  chi^2s if they can be read properly
	chi2_0=`echo| awk '{print $1}' simchi$j.dat`
	chi2sb_0=`echo| awk '{print $2}' simchi$j.dat` 
	chi2sig_0=`echo| awk '{print $3}' simchi$j.dat` 
	echo 'i is larger than 1, therefore chi2old is '
	fi
    fi
echo 'chi2old and chi2_0= ' $chi2old $chi2_0 >> Valuestokeep$m.txt
echo 'chi2old and chi2_0= ' $chi2old $chi2_0  
echo 'Acceptance test with chi^2:'>> Valuestokeep$m.txt
echo 'Acceptance test with chi^2:'
#C36 Calculating the ratio of the chi^2 and then ratio of the likelihoods
chi12=`echo "$chi2old/$chi2_0"|bc -l`
likeli12=`echo "e(-($chi2_0-$chi2old)/2.0)"|bc -l`

#C37 Deciding to accept or reject the new model. First comparing the likelihood ratio to 1.
if [ $(echo "$likeli12>1"|bc) -eq 1 ]; then
    echo 'model' $i 'accepted'
    echo 'model' $i 'accepted' $chi2old $chi2_0 $chi2sbold $chi2sb_0 $chi2sigold $chi2sig_0 $chi12 $likeli12 >> Valuestokeep$m.txt
    echo 'model' $i 'y'>> Acceptancestats.txt
    echo 'model' $i 'y'>> Acceptancestats$m.txt #Just to check that it's consistent
    accept=1
else
    x=`head -n $ran_count random.dat | tail -n 1`
    echo 'My random number x_ran ='  $x    
    echo 'My random number x_ran ='  $x 'From line= '$ran_count>> Valuestokeep$m.txt
#C38  If this fails, comparing it to a random number between zero and 1.
    if [ $(echo "$likeli12 >($x)"|bc) -eq 1 ];then
    echo $i 'accepted anyway'$x $chi2_0 $chi12
    echo 'model' $i 'accepted anyway' $x $chi2old $chi2_0 $chi2sbold $chi2sb_0 $chi2sigold $chi2sig_0 $chi12 $likeli12 >> Valuestokeep$m.txt
    accept=1
    echo 'model' $i 'y'>> Acceptancestats$m.txt
    else
    echo 'model' $i 'rejected x= ' $x $chi2old $chi2_0 $chi2sbold $chi2sb_0 $chi2sigold $chi2sig_0 $chi12 $likeli12  >> Valuestokeep$m.txt
    accept=0
    fi
    ran_count=`echo "$ran_count+1"|bc -l`
fi 

#C39 If the new model is rejected, the param$old values are left as they were
if [ $accept == 0 ]; then
    chi2old=$chi2old
    chi2sigold=$chi2sigold;
    chi2sbold=$chi2sbold;
    Msold=$Msold; 
    mualphaold=$mualphaold; 
    mudeltaold=$mudeltaold; 
    r1old=$r1old;     
    rstold=$rstold; 
    Mhaloold=$Mhaloold; 
    mdm1old=$mdm1old;
    mb1old=$mb1old; 
    nrej=`echo "$nrej+1"|bc -l`;
#C40 reprinting the old coordinates in the parameter space, as we went back to this point.
    echo 'Model rejected so chi2old kept the same= '$chi2old >>   Valuestokeep$m.txt
    echo ' Msold='$Msold  'Mualphaold='$mualphaold  'Mudeltaold='$mudeltaold 'r1old='$r1old 'Mhaloold='$Mhaloold  'mdmold='$mdm1old 'mb1old='$mb1old 'rstold='$rstold 'chi2old='$chi2old'chi2sb='$chi2sbold 'chi2sig='$chi2sigold >> Acceptedmodels.dat
    echo ' Msold='$Msold  'Mualphaold='$mualphaold  'Mudeltaold='$mudeltaold 'r1old='$r1old 'Mhaloold='$Mhaloold  'mdmold='$mdm1old 'mb1old='$mb1old 'rstold='$rstold 'chi2old='$chi2old'chi2sb='$chi2sbold 'chi2sig='$chi2sigold >> Acceptedmodels$m.dat
    echo $Msold $mb1old $rstold $Mhaloold $mdm1old $r1old $mualphaold $mudeltaold $chi2old $chi2sbold $chi2sigold >> myfinalacceptedparam$m.dat
#C41 If the new model is accepted, the param$old values are set to the values of this model so we moved in the param space
elif  [ $accept == 1 ]; then
    chi2old=$chi2_0
    chi2sbold=$chi2sb_0
    chi2sigold=$chi2sig_0
    Msold=$Msnew;
    mualphaold=$mualphanew; 
    mudeltaold=$mudeltanew; 
    r1old=$r1new;
    rstold=$rstnew; 
    Mhaloold=$Mhalonew;
    mdm1old=$mdm1new; 
    mb1old=$mb1new; 
    nacc=`echo "$nacc+1"|bc -l`;
    echo 'Model accepted so chi2old will have the value of chi2_0= '$chi2old  >>   Valuestokeep$m.txt
    echo ' Ms='$Msold  'Mualpha='$mualphaold 'Mudelta='$mudeltaold 'r1='$r1old 'Mhalo='$Mhaloold  'mdm='$mdm1old 'mb1='$mb1old 'rst='$rstold 'chi2='$chi2old 'chi2sb='$chi2sbold 'chi2sig='$chi2sigold>> Acceptedmodels$m.dat
    echo ' Ms='$Msold  'Mualpha='$mualphaold 'Mudelta='$mudeltaold 'r1='$r1old 'Mhalo='$Mhaloold  'mdm='$mdm1old 'mb1='$mb1old 'rst='$rstold 'chi2='$chi2old 'chi2sb='$chi2sbold 'chi2sig='$chi2sigold>> Acceptedmodels.dat
    echo $Msold $mb1old $rstold $Mhaloold $mdm1old $r1old $mualphaold $mudeltaold $chi2old $chi2sbold $chi2sigold >> myfinalacceptedparam$m.dat
    #C42 printing the new coordinates in the parameter space.
fi
echo ' i='$i 'Msold='$Msold 'Mualphaold='$mualphaold 'Mudeltaold='$mudeltaold 'r1old='$r1old 'Mhaloold='$Mhaloold 'mdmold='$mdm1old 'mb1old='$mb1old 'rstold='$rstold >> Valuestokeep$m.txt
echo ' i='$i 'Msold='$Msold 'Mualphaold='$mualphaold 'Mudeltaold='$mudeltaold 'r1old='$r1old  'Mhaloold='$Mhaloold 'mdmold='$mdm1old  'mb1old='$mb1old 'rstold='$rstold

#########################
# NEW PARAMETERS        #
#########################
#C43 Taking the random step from the chosen new coordinates in the param space.
#C44 In the below, it reads in a random number for each param, multiplies it with the range of the param and adds this step to the previous value of the param. If the step is in a direction to make the parameter value to be outside the allowed range, it re-iterates with the next random number in the list.

rangeacc=0
    while [ $rangeacc -lt 1 ]    
    do 
	mualpharan=`head -n $ranal_count alpharan.dat | tail -n 1`
	mualpharan=`echo "scale=20;$mualpharan"|bc -l`
	mualphanew=`echo "$mualphaold+$mualpharan*$dalpha"|bc -l`
	echo 'a mualpharan' $mualpharan $dalpha $mualphanew $a
	if [ $(echo "$mualphanew <=0.61"|bc) -eq 1 ] && [ $(echo "$mualphanew >=-0.17"|bc) -eq 1 ]; then
	    rangeacc=1
	else 
	    rangeacc=0
	fi
	ranal_count=`echo "$ranal_count+1"|bc`
    done
#C45 the allowed ranges are:0.49 to -0.05 for mu_alpha so 22+-27
##CHECK THIS ONE!!new values from 2007 are:22+13 and 24+11 so it would be   #	if [ $(echo "$mualphanew <=0.61"|bc) -eq 1 ] && [ $(echo "$mualphanew >=-0.17"|bc) -eq 1 ]; then
    rangeacc=0;
    while [ $rangeacc -lt 1 ]
    do
	mudeltaran=`head -n $randel_count deltaran.dat | tail -n 1`
	mudeltaran=`echo "scale=20;$mudeltaran"|bc -l`
	mudeltanew=`echo "$mudeltaold+$mudeltaran*$ddelta"|bc -l`
	echo 'mudeltaran' $mudeltaran $ddelta $mudeltanew
	if [ $(echo "$mudeltanew >=-0.09"|bc) -eq 1 ] && [ $(echo "$mudeltanew <=0.57"|bc) -eq 1 ]; then
	    rangeacc=1
	else 
	    rangeacc=0
	fi
    randel_count=`echo "$randel_count+1"|bc`
    done
#C46 the allowed ranges are:-0.12 to -0.42 for mu_alpha so 15+-27
    rangeacc=0;
    while [ $rangeacc -lt 1 ]
    do
	Msran=`head -n $ranms_count msran.dat | tail -n 1`
	Msran=`echo "scale=20;$Msran"|bc -l`
	Msnew=`echo "$Msold+$Msran*$dMs"|bc -l`
	echo 'Msran' $Msran $dMs $Msnew	
	if [ $(echo "$Msnew <=2150000.0"|bc) -eq 1 ] && [ $(echo "$Msnew >=430000.0"|bc) -eq 1 ]; then
	    rangeacc=1
	else 
	    rangeacc=0
	fi
    ranms_count=`echo "$ranms_count+1"|bc`
    done
#C47 the allowed ranges are M* and 5M*
    rangeacc=0;
     while [ $rangeacc -lt 1 ]
    do
	rstran=`head -n $ranrst_count rstran.dat | tail -n 1`
	rstran=`echo "scale=20;$rstran"|bc -l`
	rstnew=`echo "e(l($rstold)+$rstran*$drst)"|bc -l`
	echo 'rstran' $ranrst_count $rstran $drst $rstnew
 	if [ $(echo "$rstnew <=0.744"|bc) -eq 1 ] && [ $(echo "$rstnew >=0.074"|bc) -eq 1 ]; then
	    rangeacc=1
	else 
	    rangeacc=0
	fi
    ranrst_count=`echo "$ranrst_count+1"|bc`
    done
    rangeacc=0;
#C48 the allowed ranges are 74pc to 740pc
    while [ $rangeacc -lt 1 ]
    do
	r1ran=`head -n $ranr1_count rdmran.dat | tail -n 1`
	r1ran=`echo "scale=20;$r1ran"|bc -l`
	r1new=`echo "e(l($r1old)+$r1ran*$dr1)"|bc -l`
	echo 'r1ran' $ranr1_count $r1ran $dr1 $r1new
	if [ $(echo "$r1new <=5.0"|bc) -eq 1 ] && [ $(echo "$r1new >=$rstnew"|bc) -eq 1 ]; then
	    rangeacc=1
	else 
	    rangeacc=0
	fi
    ranr1_count=`echo "$ranr1_count+1"|bc`
    done
#C49 the allowed ranges are rsmin and 5kpc
#CHECK THIS WITH MArk, should i go all the way down to 74 pc like in stars or should i stick to 2.5
    rangeacc=0;
   while [ $rangeacc -lt 1 ]
    do
	Mhaloran=`head -n $ranmdm_count mdmran.dat | tail -n 1`
	Mhaloran=`echo "scale=20;$Mhaloran"|bc -l`
	Mhalonew=`echo "e(l($Mhaloold)+$Mhaloran*$dmhalo)"|bc -l`
	echo 'Mhaloran' $ranmdm_count $Mhaloran $dmhalo $Mhalonew
	if [ $(echo "$Mhalonew <=10000000000.0"|bc) -eq 1 ] && [ $(echo "$Mhalonew >=1000000.0"|bc) -eq 1 ]; then
	    rangeacc=1
	else 
	    rangeacc=0
	fi
    ranmdm_count=`echo "$ranmdm_count+1"|bc`
#C50 the allowed ranges are 10^6 to 10^10  Check THIS!!!

    done
	mdm1new=`echo "$Mhalonew/222288.0"|bc -l`
	echo 'mdm1new '$mdm1new
	mb1new=`echo "$Msnew/222288.0"|bc -l`
	echo 'mb1new '$mb1new

echo 'mualpharan dalpha mualphanew' $mualpharan $dalpha $mualphanew >> Valuestokeep$m.txt
echo 'mudeltaran ddelta mudeltanew' $mudeltaran $ddelta $mudeltanew >> Valuestokeep$m.txt
echo 'Msran dMs Msnew' $Msran $dMs $Msnew  >> Valuestokeep$m.txt
echo 'r1ran dr1 r1new' $r1ran $dr1 $r1new  >> Valuestokeep$m.txt
echo 'Mhaloran dmhalo Mhalonew' $Mhaloran $dmhalo $Mhalonew  >> Valuestokeep$m.txt
echo 'rstran drst rstnew' $rstran $drst $rstnew  >> Valuestokeep$m.txt
#C51 wrote the new values to be tried

mv newnemscr.sh oldnemscr$i.sh
rm B0.*
rm H0.* 

#C52 Writing the new tcsh.sh file for nemo 

echo  '#!/bin/tcsh
if(! $?NEMO) then
setenv NEMO /data/dp005/dc-ural1/nemo_cvs
module load gcc
source $NEMO/nemo_start
endif
setenv FALCON $NEMO/usr/dehnen/falcON
source ${FALCON}/falcON_start
cd ' $SIMDIR '
mkhalo B'$i'.ini 100000 inner=0.515 outer=4.45 r_s='$rstnew 'M='$Msnew 'eta=3.745 b=0.0 r_a=0.0 accpars=0.0,'$r1new','$mdm1new',1.0,4.0,1.0,0.0,0.0 accname=Halo tabfile=B1.tab r_t=10 WD_units=t
mkhalo H'$i'.ini 100000 inner=1.0 outer=4.0 r_s='$r1new' M='$Mhalonew' eta=1.0 b=0.0 r_a=0.0 accpars=0.0,'$rstnew','$mb1new',0.515,4.45,3.745,0.0,0.0 accname=Halo tabfile=H1.tab r_t=10 WD_units=t 
s2a B'$i'.ini B'$i'.dat
s2a H'$i'.ini H'$i'.dat
source $FALCON/falcON_end
source $NEMO/nemo_end
' >> newnemscr.sh 

#C53 Running NEMO (newnemscr.sh) to generate the next model galaxy.
chmod +x newnemscr.sh
./newnemscr.sh

#C54 Cheking whether Nemo could find equilibrium models. The if  GOES WITH else at C61 and  fi AT C73
if [ -e B$i.dat ] && [ -e H$i.dat ]  ; then
awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' B$i.dat > bulgeo
awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' H$i.dat > haloo
#C55 awk to write the bulgeo and haloo ready for the conversion code (reads the 7 columned part of NEMO files.)


#C58 Running the unit conversion code reading the bulgeo and haloo and the orbital parameters for the circular orbit
# This code writes the new coordinates using the new c.o.m. and converting the vel. unit from nemo(~km/s) to pkdgrav(~0.65km/s) out to the files halo and bulge
python myunits5.py 200000 100000 100000 oroutcirc
#C59 Running a2b using the halo and bulge files  converted and centred on the far orbit.
./a2b 1.0 1.0
#C60 Printing the first binary file initial.bin and tagging the halo and bulge files in the chains
echo 'nemobitti'
else
#C61 If Nemo couldn't find equilibrium models, either for bulge or halo, it removes the generated component, and tries a new model with different parameters.
    newfile=0
    while [ $newfile -lt 1 ]
    do
    rm B$i.dat
    rm H$i.dat
    rm B$i.ini
    rm H$i.ini
    rm newnemscr.sh
    echo 'One of the files does not exist. Retrying.' $i 
    echo 'One of the files does not exist. Retrying.' $i >> Valuestokeep$m.txt
    ranal_count=`echo "$ranal_count+1"|bc -l`
    randel_count=`echo "$randel_count+1"|bc -l`
    ranr1_count=`echo "$ranr1_count+1"|bc -l`
    ranmdm_count=`echo "$ranmdm_count+1"|bc -l`
    ranms_count=`echo "$ranms_count+1"|bc -l`
    ranrst_count=`echo "$ranrst_count+1"|bc -l`
    rangeacc=0;
      	while [ $rangeacc -lt 1 ]
	do 
	    mualpharan=`head -n $ranal_count alpharan.dat | tail -n 1`
	    mualpharan=`echo "scale=20;($mualpharan)"|bc -l`
	    mualphanew=`echo "$mualphaold+$mualpharan*$dalpha"|bc -l`
	    echo 'a mualpharan' $mualpharan $dalpha $mualphanew $a
	    if [ $(echo "$mualphanew <=0.61"|bc) -eq 1 ] && [ $(echo "$mualphanew >=-0.17"|bc) -eq 1 ]; then
		rangeacc=1
	    else 
		rangeacc=0
	    fi
	    ranal_count=`echo "$ranal_count+1"|bc`
	done
#C62 the allowed ranges are:0.49 to -0.05 for mu_alpha so 22+-27
	rangeacc=0;
	while [ $rangeacc -lt 1 ]
	do
	    mudeltaran=`head -n $randel_count deltaran.dat | tail -n 1`
	    mudeltaran=`echo "scale=20;($mudeltaran)"|bc -l`
	    mudeltanew=`echo "$mudeltaold+$mudeltaran*$ddelta"|bc -l`
	    echo 'mudeltaran' $mudeltaran $ddelta $mudeltanew
	    if [ $(echo "$mudeltanew >=-0.09"|bc) -eq 1 ] && [ $(echo "$mudeltanew <=0.57"|bc) -eq 1 ]; then
		rangeacc=1
	    else 
		rangeacc=0
	    fi
	    randel_count=`echo "$randel_count+1"|bc`
	done
#C63 the allowed ranges are:-0.12 to -0.42 for mu_alpha so 15+-27
	rangeacc=0;
	while [ $rangeacc -lt 1 ]
	do
	    Msran=`head -n $ranms_count msran.dat | tail -n 1`
	    Msran=`echo "scale=20;$Msran"|bc -l`
	    Msnew=`echo "$Msold+$Msran*$dMs"|bc -l`
	    echo 'Msran' $Msran $dMs $Msnew	
	    if [ $(echo "$Msnew <=2150000.0"|bc) -eq 1 ] && [ $(echo "$Msnew >=430000.0"|bc) -eq 1 ]; then
		rangeacc=1
	    else 
		rangeacc=0
	    fi
	    ranms_count=`echo "$ranms_count+1"|bc`
	done
#C64 the allowed ranges are M* and 5M*
	rangeacc=0;
	while [ $rangeacc -lt 1 ]
	do
	    rstran=`head -n $ranrst_count rstran.dat | tail -n 1`
	    rstran=`echo "scale=20;$rstran"|bc -l`
	    rstnew=`echo "e(l($rstold)+$rstran*$drst)"|bc -l`
	    echo 'rstran' $ranrst_count $rstran $drst $rstnew
	    if [ $(echo "$rstnew <=0.744"|bc) -eq 1 ] && [ $(echo "$rstnew >=0.074"|bc) -eq 1 ]; then
		rangeacc=1
	    else 
		rangeacc=0
	    fi
	    ranrst_count=`echo "$ranrst_count+1"|bc`
	done
#C65 the allowed ranges are 74pc to 740pc
	rangeacc=0;	
	while [ $rangeacc -lt 1 ]
	do
	    r1ran=`head -n $ranr1_count rdmran.dat | tail -n 1`
	    r1ran=`echo "scale=20;$r1ran"|bc -l`
	    r1new=`echo "e(l($r1old)+$r1ran*$dr1)"|bc -l`
	    echo 'r1ran' $r1ran $dr1 $r1new
	    if [ $(echo "$r1new <=5.0"|bc) -eq 1 ] && [ $(echo "$r1new >=$rstnew"|bc) -eq 1 ]; then
		rangeacc=1
	    else 
		rangeacc=0
	    fi
	        ranr1_count=`echo "$ranr1_count+1"|bc`
	done
#C66 the allowed ranges are rsmin and 5kpc
#CHECK THIS WITH MArk, should i go all the way down to 74 pc like in stars or should i stick to 2.5
	rangeacc=0;
	while [ $rangeacc -lt 1 ]
	do
	    Mhaloran=`head -n $ranmdm_count mdmran.dat | tail -n 1`
	    Mhaloran=`echo "scale=20;$Mhaloran"|bc -l`
	    Mhalonew=`echo "e(l($Mhaloold)+$Mhaloran*$dmhalo)"|bc -l`
	    echo 'Mhaloran' $Mhaloran $dmhalo $Mhalonew
	    if [ $(echo "$Mhalonew <=10000000000.0"|bc) -eq 1 ] && [ $(echo "$Mhalonew >=1000000.0"|bc) -eq 1 ]; then
		rangeacc=1
	    else 
		rangeacc=0
	    fi
	    ranmdm_count=`echo "$ranmdm_count+1"|bc`
#C67 the allowed ranges are 10^6 to 10^10  Check THIS!!!
	done
	mdm1new=`echo "$Mhalonew/222288.0"|bc -l`
	echo 'mdm1new '$mdm1new
	mb1new=`echo "$Msnew/222288.0"|bc -l`

    echo 'mb1new '$mb1new
    echo 'mualpharan dalpha mualphanew' $mualpharan $dalpha $mualphanew >> Valuestokeep$m.txt
    echo 'mudeltaran ddelta mudeltanew' $mudeltaran $ddelta $mudeltanew >> Valuestokeep$m.txt
    echo 'Msran dMs Msnew' $Msran $dMs $Msnew  >> Valuestokeep$m.txt
    echo 'r1ran dr1 r1new' $r1ran $dr1 $r1new  >> Valuestokeep$m.txt
    echo 'Mhaloran dmhalo Mhalonew' $Mhaloran $dmhalo $Mhalonew  >> Valuestokeep$m.txt 
    echo 'rstran drst rstnew' $rstran $drst $rstnew  >> Valuestokeep$m.txt
#C68 wrote the new values to be tried

#C69 Writing the new tcsh.sh file for nemo    
echo  '#!/bin/tcsh
if(! $?NEMO) then
setenv NEMO /data/dp005/dc-ural1/nemo_cvs
module load gcc
source $NEMO/nemo_start
endif
setenv FALCON $NEMO/usr/dehnen/falcON
source ${FALCON}/falcON_start
cd ' $SIMDIR '
mkhalo B'$i'.ini 100000 inner=0.515 outer=4.45 r_s='$rstnew 'M='$Msnew 'eta=3.745 b=0.0 r_a=0.0 accpars=0.0,'$r1new','$mdm1new',1.0,4.0,1.0,0.0,0.0 accname=Halo tabfile=B1.tab r_t=10 WD_units=t
mkhalo H'$i'.ini 100000 inner=1.0 outer=4.0 r_s='$r1new' M='$Mhalonew' eta=1.0 b=0.0 r_a=0.0 accpars=0.0,'$rstnew','$mb1new',0.515,4.45,3.745,0.0,0.0 accname=Halo tabfile=H1.tab r_t=10 WD_units=t 
s2a B'$i'.ini B'$i'.dat
s2a H'$i'.ini H'$i'.dat
source $FALCON/falcON_end
source $NEMO/nemo_end
' >> newnemscr.sh 

#C70 Running NEMO (newnemscr.sh) to generate the next model galaxy.
    chmod +x newnemscr.sh
    ./newnemscr.sh
#C71 Cheking whether Nemo could find equilibrium models. If bulge or halo could both be generated it continues with writing files. 
    if [ -e B$i.dat ] && [ -e H$i.dat ]  ; then
    newfile=1
    else
    newfile=0
    rm B$i.dat 
    rm H$i.dat
    rm B$i.ini
    rm H$i.ini
    fi
    done
#C72 Set up newfile =0 or 1 so that we know if nemo worked or if we need to choose yet another model for equilibrium. 
    awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' B$i.dat > bulgeo
    awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' H$i.dat > haloo
    echo 'intcsh'
    python myunits5.py 200000 100000 100000 oroutcirc
    ./a2b 1.0 1.0
    echo 'nemobitti'
#C73 THE fi goes with if at C54 This stays at newfile=0 which follows the else in C61 unless newfile becomes 1 at C71
fi
#C73a printing the myinitialnoexpparam1.dat for the restarting of the chains.
echo $Msnew $mb1new $rstnew $Mhalonew $mdm1new $r1new $mualphanew $mudeltanew $i>>myinitialnoexpparam$m.dat
echo $Msnew $mb1new $rstnew $Mhalonew $mdm1new $r1new $mualphanew $mudeltanew $i $halooin$i>>textformyinitialparam$m.dat

#C74 Moving the generation files to files$i before running the simulation that is coming
mv haloo halooin$i
mv bulgeo bulgeoin$i
mv halo haloin$i
mv bulge bulgein$i
mpiexec  /home/dc-ural1/mpkdgravopmpinopottry_2012/pkdgrav.d.gnu simeq.param
mv initial.bin initialeq$i.bin
mv initialeq.00005 initialeq$i.00005
python centreofmass.py initialeq$i.00005 com_file.dat
mv com_file.dat comeq5$i.dat 

python checkinginitialbin2.py initialeq$i.00005 initialeq$i.00005.dat
head -n 100000 initialeq$i.00005.dat >> haloo
tail -n 100000 initialeq$i.00005.dat >> bulgeo
rm initialeq$i.00005.dat

./adconvugur6gyrs 100.4029 -50.9661 102.0 $mualphanew $mudeltanew 223.1
mv convvel.dat orbx$i.dat 
python Ugur_orbit1.py orbx$i.dat orbxpast$i.dat
python myunits6goingtozero.py 200000 100000 100000 orbxpast$i.dat comeq5$i.dat

./a2b 1.0 1.0
mv haloo halooorb$i
mv halo halorb$i
mv bulgeo bulgeoorb$i
mv bulge bulgeorb$i
mpiexec  /home/dc-ural1/mpkdgravopmpinopottry_2012/pkdgrav.d.gnu sim1.param 
echo 'i'$i initialeq$i.00005
mv initial.bin initial$i.bin
mv final.00060 final$i.00060
echo 'Im haaaaapppy'
python dentry4coredart_2012.py final$i.00060 simb$i.dat simh$i.dat simb$i.png simh$i.png
python dispGitHUB.py final$i.00060 simdisp$i.dat simdisp$i.png
python checkinginitialbin2.py final$i.00060 final$i.00060.dat
tail -n 100000 final$i.00060.dat >> final$i.00060_stars.dat 
./posangle6_2013 100000 final$i.00060_stars.dat simpa$i.dat com_file.dat
python meanvelnew.py simpa$i.dat  simgrad$i.dat 100000
python stattestGITHUBart.py obsdisp146noisy180.dat obsb146noisy141.dat simdisp$i.dat simb$i.dat simgrad$i.dat simchi$i.dat
#pythonfittingdatacored12.py simb$i.dat simh$i.dat ofit$i.dat ofit$i.png
rm final*dat
mv simb*dat PROFILES/.
mv simh*dat PROFILES/.
mv simdisp*dat PROFILES/.
mv simgrad*dat PROFILES/.
mv simpa*dat PROFILES/.
#C75 Doing the batchtest using the acceptance ratio, though I don't think it works!
batchtest=$[ $i % 50 ]
echo 'batchtest=' $batchtest 
#C76 cheks whether we already ran 50 models
if [ $i> 0 ] && [ $batchtest= 0 ]; then
    nbatch=`echo "$nbatch+1"|bc`
    ratioacc=`echo "scale=20;$nacc/($nacc+$nrej)"|bc -l`	
    echo "Acceptance ratio= "$ratioacc
    echo "nbatch" $nbatch  
    echo ' Step change factor' 
	dstep=`echo "scale=20;1/sqrt($nbatch)"|bc -l`
	dstepneg=`echo "scale=20;$dstep > 0.02"|bc -l`
	if [ $dstepneg> 0 ]; then
	    dstep=`echo "scale=20;0.01"|bc -l`
	fi
	checkacc=`echo "scale=20;$ratioacc>0.44"|bc -l`
	if  [ $checkacc< 1 ]; then
	    dstep=`echo "scale=20;$dstep*(-1.)"|bc -l`  
	    echo 'accepted too few,  decreasing the step size '$dstep
	else
	    dstep=`echo "scale=20;$dstep*(1.)"|bc -l`
	    echo 'accepted too many, increasing the step size '$dstep
	fi
	echo ' Acceptance ratio='$ratioacc >> Valuestokeep$m.txt
	echo 'nbatch='$nbatch >> Valuestokeep$m.txt
	echo 'Step change='$dstep  >> Valuestokeep$m.txt
	echo 'Old step sizes: dalpha='$dalpha 'ddelta='$ddelta 'dMs='$dMs 'dr1='$dr1 'dmhalo='$dmhalo>> Valuestokeep$m.txt
	dalpha=`echo "scale=10;e(((l($dalpha^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	ddelta=`echo "scale=10;e(((l($ddelta^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	dMs=`echo "scale=10;e(((l($dMs^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	dr1=`echo "scale=10;e(((l($dr1^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	drst=`echo "scale=10;e(((l($drst^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	dmhalo=`echo "scale=10;e(((l($dmhalo^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	echo 'The step sizes changed in this one. 
              New step sizes: dalpha='$dalpha 'ddelta='$ddelta 'dMs='$dMs 'dr1='$dr1 'drst='$drst 'dmhalo='$dmhalo>> Valuestokeep$m.txt
	echo 'The step sizes changed in this one'
fi
echo 'If the new step sizes arent quoted by now it means that they did not change...'

i=$(echo "$i+1"|bc)
    echo 'final i='$i
      if [ $i == 2000 ]; then
	  finished=1
	  echo 'finished'
      fi
done


