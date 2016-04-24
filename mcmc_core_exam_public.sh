#PBS -l nodes=2:ppn=16
#PBS -l walltime=5:00:00:00
#PBS -l pvmem=1gb
#-o  $JOB_NAME.o$JOB_ID
#-e  $JOB_NAME.e$JOB_ID
#PBS -v PYTHONPATH=/home/orblib/
#PBS -S /bin/bash
#======================================
# Set up environment
#======================================
module load openmpi
module load gcc
module load python/2.7.3
module load swig
#======================================

SIMDIR="ENTER YOUR DIRECTORY"
PKDGRAV=/home/pkdgrav
cd $SIMDIR

rm -f newnemscr.sh

# Generating random numbers
python myrangen.py   

#################################################
# Read and set the values for initial parameters#
#################################################
INPARFILE=$SIMDIR/initialparam.dat
Ms=`echo |awk '{print $1}' $INPARFILE `    #M_S
mb=`echo |awk '{print $2}' $INPARFILE `   #M_S/222288
rst=`echo |awk '{print $3}' $INPARFILE `   #r_s
Mhalo=`echo |awk '{print $4}' $INPARFILE ` #M_DM
mdm=`echo |awk '{print $5}' $INPARFILE `  #M_DM/222288
r1=`echo |awk '{print $6}' $INPARFILE `    #r_dm        
mualpha=`echo |awk '{print $7}' $INPARFILE ` #mu_alpha*cos(delta) 
mudelta=`echo |awk '{print $8}' $INPARFILE ` #mu_delta

#C Setting the initial old/new parameters
Msold=`echo "$Ms"|bc -l`
mbold=`echo "$mb"|bc -l`
rstold=`echo "$rst"|bc -l`
Mhaloold=`echo "$Mhalo"|bc -l`
mdmold=`echo "$mdm"|bc -l`
r1old=`echo "$r1"|bc -l`
mualphaold=`echo "$mualpha"|bc -l`
mudeltaold=`echo "$mudelta"|bc -l`
Msnew=`echo "$Ms"|bc -l`
mbnew=`echo "$mb"|bc -l`
rstnew=`echo "$rst"|bc -l` 
Mhalonew=`echo "$Mhalo"|bc -l`
mdmnew=`echo "$mdm"|bc -l`
r1new=`echo "$r1"|bc -l`
mualphanew=`echo "$mualpha"|bc -l`
mudeltanew=`echo "$mudelta"|bc -l`

echo 'Initial values Ms mb rst Mhalo mdm r1 mualpha mudelta' $Ms $mb $rst $Mhalo $mdm $r1 $mualpha $mudelta>> Valuestokeep.txt

###########################################################
# Calling tcsh for Nemo and running the first NEMO script #
###########################################################

echo  '#!/bin/tcsh
if(! $?NEMO) then
setenv NEMO /data/dp005/dc-ural1/nemo_cvs
module load gcc
source $NEMO/nemo_start
endif
setenv FALCON $NEMO/usr/dehnen/falcON
source ${FALCON}/falcON_start
cd ' $SIMDIR '
mkhalo B0.ini 100000 inner=0.515 outer=4.45 r_s='$rstold 'M='$Msold 'eta=3.745 b=0.0 r_a=0.0 accpars=0.0,'$r1','$mdmold',0.0,4.0,1.0,0.0,0.0 accname=Halo tabfile=B1.tab r_t=10 WD_units=t
mkhalo H0.ini 100000 inner=0.0 outer=4.0 r_s='$r1old' M='$Mhaloold' eta=1.0 b=0.0 r_a=0.0 accpars=0.0,'$rstold','$mbold',0.515,4.45,3.745,0.0,0.0 accname=Halo tabfile=H1.tab r_t=10 WD_units=t 
s2a B0.ini B0.dat
s2a H0.ini H0.dat
source $FALCON/falcON_end
 source $NEMO/nemo_end
' >> newnemscr.sh 

chmod +x newnemscr.sh
./newnemscr.sh

# NEMO units     M=222288msun G=1 r=kpc t=gyrs 
# PKDRAV units   M=10^5Msun G=1 r=kpc T=1.490*10^9yrs
awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' B0.dat > bulgeo
awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' H0.dat > haloo

echo 'initial galaxy generated'
echo '300.0 0.0 0.0 0.0 0.1475 0.0 0.00001 5000.0' > orfcirc

python Car_orbit.py orfcirc oroutcirc

python units1public.py 200000 100000 100000 oroutcirc
./ascii2tipsy 1.0 1.0

mv haloo halooin0
mv bulgeo bulgeoin0
mv halo haloin0
mv bulge bulgein0

mpirun PKDGRAV simeq.param 
mv initial.bin initialeq0.bin
python compublic.py initialeq0.bin comeq0.dat

python bin2asciipublic.py initialeq.00005 initialeq0.00005.dat
# outputs the binary file to an ascii file with m,x,y,z,vx,vy,yz,eps
head -n 100000 initialeq0.00005.dat >> haloo
tail -n 100000 initialeq0.00005.dat >> bulgeo

mv initialeq.00005 initialeq0.00005
python compublic.py initialeq0.00005 comeq5.dat
./adconv6gyrs 100.4029 -50.9661 100.9 $mualphaold $mudeltaold 223.1

mv convvel.dat orbx0.dat 
python Car_orbit.py orbx0.dat orbxpast0.dat
python units2public.py 200000 100000 100000 orbxpast0.dat comeq5.dat

./ascii2tipsy 1.0 1.0
mv haloo halooorb0
mv halo halorb0
mv bulgeo bulgeoorb0
mv bulge bulgeorb0
# Printing the binary file initial.bin for the first simulation and on orbit
mpirun PKDGRAV sim1.param 

mv initial.bin initial0.bin
mv final.00060 final0.00060

#ANALYSIS OF SIMULATION
python sbpublic.py final0.00060 simb0.dat simh0.dat simb0.png simh0.png
python disppublic.py final0.00060 simdisp0.dat simdisp0.png
python bin2asciipublic.py final0.00060 final0.00060.dat
tail -n 100000 final0.00060.dat >> final0.00060_stars.dat 
./papublic 100000 final0.00060_stars.dat simpa0.dat com_file.dat
python gradpublic.py simpa0.dat simgrad0.dat 100000 
python statpublic.py obsdispnoisy.dat obsbnoisy.dat simdisp0.dat simb0.dat simgrad0.dat simchi0.dat

# Chosing the step sizes for the parameters of the MCMC
dalpha=`echo "0.078"|bc -l` 
ddelta=`echo "0.066"|bc -l` 
dMs=`echo "215000.0"|bc -l` 
dr1=`echo "0.23"|bc -l` 
dmhalo=`echo "0.92"|bc -l`
drst=`echo "0.23"|bc -l` 

####################
#  LOOP            #
####################
# Setting the chi^2 to zero so that the first model is always accepted.
finished=0
  i=1
  seed=`echo "-6256326852"|bc -l`
  pold=`echo "-0.1"|bc -l`
  chi2old=`echo "3000.0"|bc -l`
  chi2sigold=`echo "3000.0"|bc -l`
  chi2sbold=`echo "3000.0"|bc -l`

  nacc=`echo "0.0"|bc -l`
  nrej=`echo "0.0"|bc -l`
  nbatch=0

chi2_0=`echo| awk '{print $1}' simchi0.dat`
chi2sb_0=`echo| awk '{print $2}' simchi0.dat`
chi2sig_0=`echo| awk '{print $3}' simchi0.dat`
echo 1 |awk '{srand('$seed');for(i=1;i<=5000;i++){printf("%lf\n",rand())}}' >>random.dat

ran_count=1
ranal_count=1
randel_count=1
ranms_count=1
ranr1_count=1
ranmdm_count=1
ranrst_count=1

################
#  COMPARISON  #
################
while [ $finished -lt 1 ]
    do    
mv *png FIGURES/.
mv comeq5*dat CENTREOFMASSES/.
mv orbx*dat ORBITS/.
mv bulge* BULGES/.
mv halo* HALOS/.
mv B* BULGES/.
mv H* HALOS/.
mv oldnem* GENSCRIPTS/.
mv initial*5 INFILES/.
mv initial*bin INFILES/.
mv final*60 FINALFILES/.

    if [ $i == 1 ]; then
	chi2old=$chi2old;	
	chi2sigold=$chi2sbold;	
	chi2sbold=$chi2sigold;	
	chi2_0=`echo| awk '{print $1}' simchi0.dat`
	chi2sb_0=`echo| awk '{print $2}' simchi0.dat`
	chi2sig_0=`echo| awk '{print $3}' simchi0.dat`
	echo 'i is 1, therefore chi2old and chi2_0 are' $chi2old $chi2_0 $chi2sbold $chi2sb $chi2sigold $chi2sig   
    fi
    if [ $i> 1 ]; then
	j=`echo "$i-1"|bc -l`
	chi2_0=`echo| awk '{print $1}' simchi$j.dat` 
	chi2sb_0=`echo| awk '{print $2}' simchi$j.dat` 
	chi2sig_0=`echo| awk '{print $3}' simchi$j.dat` 
	if [ -s  $chi2_0 ]  ; then
	    chi2_0=`echo "3000.0"|bc -l `
	    chi2sb_0=`echo "3000.0"|bc -l `
	    chi2sig_0=`echo "3000.0"|bc -l `
	    echo 'i is larger than 1, the chi^2 could not be calculated and was put back to 3000 '>> Valuestokeep.txt 
	    mv simb$j.dat CRASHEDMODELS/.
	    mv simh$j.dat CRASHEDMODELS/.
	    mv simdisp$j.dat CRASHEDMODELS/.
	    mv simchi$j.dat CRASHEDMODELS/.
	else
#C35 Setting the  chi^2s if they can be read properly
	chi2_0=`echo| awk '{print $1}' simchi$j.dat`
	chi2sb_0=`echo| awk '{print $2}' simchi$j.dat` 
	chi2sig_0=`echo| awk '{print $3}' simchi$j.dat` 
	echo 'i is larger than 1, and chi2old is calculated '
	fi
    fi
echo 'chi2old and chi2_0= ' $chi2old $chi2_0 >> Valuestokeep.txt
echo 'Acceptance test with chi^2:'>> Valuestokeep.txt

chi12=`echo "$chi2old/$chi2_0"|bc -l`
likeli12=`echo "e(-($chi2_0-$chi2old)/2.0)"|bc -l`

# Accept or reject the new model
if [ $(echo "$likeli12>1"|bc) -eq 1 ]; then
    echo 'model' $i 'accepted' $chi2old $chi2_0 $chi2sbold $chi2sb_0 $chi2sigold $chi2sig_0 $chi12 $likeli12 >> Valuestokeep.txt
    echo 'model' $i 'y'>> Acceptancestats.txt
    accept=1
else
    x=`head -n $ran_count random.dat | tail -n 1`
    echo 'My random number x_ran ='  $x 'From line= '$ran_count>> Valuestokeep.txt
# If this fails, comparing it to a random number between zero and 1.
    if [ $(echo "$likeli12 >($x)"|bc) -eq 1 ];then
	echo 'model' $i 'accepted anyway' $x $chi2old $chi2_0 $chi2sbold $chi2sb_0 $chi2sigold $chi2sig_0 $chi12 $likeli12 >> Valuestokeep.txt
    accept=1
    echo 'model' $i 'y'>> Acceptancestats.txt
    else
    echo 'model' $i 'rejected x= ' $x $chi2old $chi2_0 $chi2sbold $chi2sb_0 $chi2sigold $chi2sig_0 $chi12 $likeli12  >> Valuestokeep.txt
    accept=0
    fi
    ran_count=`echo "$ran_count+1"|bc -l`
fi 
# If the new model is rejected, the param$old values are left as they were
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
    mdmold=$mdmold;
    mbold=$mbold; 
    nrej=`echo "$nrej+1"|bc -l`;

    echo 'Model rejected so chi2old kept the same= '$chi2old >>   Valuestokeep.txt
    echo ' Msold='$Msold  'Mualphaold='$mualphaold  'Mudeltaold='$mudeltaold 'r1old='$r1old 'Mhaloold='$Mhaloold  'mdmold='$mdmold 'mbold='$mbold 'rstold='$rstold 'chi2old='$chi2old'chi2sb='$chi2sbold 'chi2sig='$chi2sigold >> Acceptedmodels.dat
    echo $Msold $mbold $rstold $Mhaloold $mdmold $r1old $mualphaold $mudeltaold $chi2old $chi2sbold $chisigold>> myfinalacceptedparam1.dat

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
    mdmold=$mdmnew; 
    mbold=$mbnew; 
    nacc=`echo "$nacc+1"|bc -l`;
    echo 'Model accepted so chi2old will have the value of chi2_0= '$chi2old  >>   Valuestokeep.txt
    echo ' Ms='$Msold  'Mualpha='$mualphaold 'Mudelta='$mudeltaold 'r1='$r1old 'Mhalo='$Mhaloold  'mdm='$mdmold 'mb='$mbold 'rst='$rstold 'chi2='$chi2old 'chi2sb='$chi2sbold 'chi2sig='$chi2sigold>> Acceptedmodels.dat
     echo $Msold $mbold $rstold $Mhaloold $mdmold $r1old $mualphaold $mudeltaold $chi2old $chi2sbold $chi2sigold>> myfinalacceptedparam1.dat
fi
echo ' i='$i 'Msold='$Msold 'Mualphaold='$mualphaold 'Mudeltaold='$mudeltaold 'r1old='$r1old 'Mhaloold='$Mhaloold 'mdmold='$mdmold 'mbold='$mbold 'rstold='$rstold >> Valuestokeep.txt
echo ' i='$i 'Msold='$Msold 'Mualphaold='$mualphaold 'Mudeltaold='$mudeltaold 'r1old='$r1old  'Mhaloold='$Mhaloold 'mdmold='$mdmold  'mbold='$mbold 'rstold='$rstold

#########################
# NEW PARAMETERS        #
#########################
#C43 Taking the random step from the chosen new coordinates in the param space.
#C44 In the below, I read in a random number for each param, multiply it with the range of the param and add this step to the previous value of the param. If the step is in a direction to make the parameter value to be outside the allowed range, I try agin with the next random number in the list.

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

    done
	mdmnew=`echo "$Mhalonew/222288.0"|bc -l`
	echo 'mdmnew '$mdmnew
	mbnew=`echo "$Msnew/222288.0"|bc -l`
	echo 'mbnew '$mbnew

echo 'mualpharan dalpha mualphanew' $mualpharan $dalpha $mualphanew >> Valuestokeep.txt
echo 'mudeltaran ddelta mudeltanew' $mudeltaran $ddelta $mudeltanew >> Valuestokeep.txt
echo 'Msran dMs Msnew' $Msran $dMs $Msnew  >> Valuestokeep.txt
echo 'r1ran dr1 r1new' $r1ran $dr1 $r1new  >> Valuestokeep.txt
echo 'Mhaloran dmhalo Mhalonew' $Mhaloran $dmhalo $Mhalonew  >> Valuestokeep.txt
echo 'rstran drst rstnew' $rstran $drst $rstnew  >> Valuestokeep.txt

mv newnemscr.sh oldnemscr$i.sh
rm B0.*
rm H0.*B0

# Writing the new tcsh.sh file for nemo 

echo  '#!/bin/tcsh
if(! $?NEMO) then
setenv NEMO /data/dp005/dc-ural1/nemo_cvs
module load gcc
source $NEMO/nemo_start
endif
setenv FALCON $NEMO/usr/dehnen/falcON
source ${FALCON}/falcON_start
cd ' $SIMDIR '
mkhalo B'$i'.ini 100000 inner=0.515 outer=4.45 r_s='$rstnew 'M='$Msnew 'eta=3.745 b=0.0 r_a=0.0 accpars=0.0,'$r1new','$mdmnew',0.0,4.0,1.0,0.0,0.0 accname=Halo tabfile=B1.tab r_t=10 WD_units=t
mkhalo H'$i'.ini 100000 inner=0.0 outer=4.0 r_s='$r1new' M='$Mhalonew' eta=1.0 b=0.0 r_a=0.0 accpars=0.0,'$rstnew','$mbnew',0.515,4.45,3.745,0.0,0.0 accname=Halo tabfile=H1.tab r_t=10 WD_units=t 
s2a B'$i'.ini B'$i'.dat
s2a H'$i'.ini H'$i'.dat
source $FALCON/falcON_end
source $NEMO/nemo_end
' >> newnemscr.sh 

chmod +x newnemscr.sh
./newnemscr.sh

if [ -e B$i.dat ] && [ -e H$i.dat ]  ; then
awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' B$i.dat > bulgeo
awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' H$i.dat > haloo
echo 'intcsh'

python units1public.py 200000 100000 100000 oroutcirc
./ascii2tipsy 1.0 1.0
echo 'nemo ended'
else
    newfile=0
    while [ $newfile -lt 1 ]
    do
    rm B$i.dat
    rm H$i.dat
    rm B$i.ini
    rm H$i.ini
    rm newnemscr.sh
    echo 'One of the files does not exist. Retrying.' $i 
    echo 'One of the files does not exist. Retrying.' $i >> Valuestokeep.txt
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

	done
	mdmnew=`echo "$Mhalonew/222288.0"|bc -l`
	echo 'mdmnew '$mdmnew
	mbnew=`echo "$Msnew/222288.0"|bc -l`

    echo 'mbnew '$mbnew
    echo 'mualpharan dalpha mualphanew' $mualpharan $dalpha $mualphanew >> Valuestokeep.txt
    echo 'mudeltaran ddelta mudeltanew' $mudeltaran $ddelta $mudeltanew >> Valuestokeep.txt
    echo 'Msran dMs Msnew' $Msran $dMs $Msnew  >> Valuestokeep.txt
    echo 'r1ran dr1 r1new' $r1ran $dr1 $r1new  >> Valuestokeep.txt
    echo 'Mhaloran dmhalo Mhalonew' $Mhaloran $dmhalo $Mhalonew  >> Valuestokeep.txt 
    echo 'rstran drst rstnew' $rstran $drst $rstnew  >> Valuestokeep.txt

#Writing the new tcsh.sh file for nemo    
echo  '#!/bin/tcsh
if(! $?NEMO) then
setenv NEMO /data/dp005/dc-ural1/nemo_cvs
module load gcc
source $NEMO/nemo_start
endif
setenv FALCON $NEMO/usr/dehnen/falcON
source ${FALCON}/falcON_start
cd ' $SIMDIR '
mkhalo B'$i'.ini 100000 inner=0.515 outer=4.45 r_s='$rstnew 'M='$Msnew 'eta=3.745 b=0.0 r_a=0.0 accpars=0.0,'$r1new','$mdmnew',0.0,4.0,1.0,0.0,0.0 accname=Halo tabfile=B1.tab r_t=10 WD_units=t
mkhalo H'$i'.ini 100000 inner=0.0 outer=4.0 r_s='$r1new' M='$Mhalonew' eta=1.0 b=0.0 r_a=0.0 accpars=0.0,'$rstnew','$mbnew',0.515,4.45,3.745,0.0,0.0 accname=Halo tabfile=H1.tab r_t=10 WD_units=t 
s2a B'$i'.ini B'$i'.dat
s2a H'$i'.ini H'$i'.dat
source $FALCON/falcON_end
source $NEMO/nemo_end
' >> newnemscr.sh 

    chmod +x newnemscr.sh
    ./newnemscr.sh
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

    awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' B$i.dat > bulgeo
    awk '(NF==7 && $1!="#"){printf("%E %E %E %E %E %E %E\n",$1*2.22288,$2,$3,$4,$5*1.49081988,$6*1.49081988,$7*1.49081988)}' H$i.dat > haloo
    echo 'intcsh'
    python units1public.py 200000 100000 100000 oroutcirc
    ./ascii2tipsy 1.0 1.0
    echo 'nemo ended'
fi
echo $Msnew $mbnew $rstnew $Mhalonew $mdmnew $r1new $mualphanew $mudeltanew $i>>myinitialnoexpparam1.dat
echo $Msnew $mbnew $rstnew $Mhalonew $mdmnew $r1new $mualphanew $mudeltanew $i $halooin$i>>textformyinitialparam.dat
mv haloo halooin$i
mv bulgeo bulgeoin$i
mv halo haloin$i
mv bulge bulgein$i
mpiexec PKDGRAV simeq.param
mv initial.bin initialeq$i.bin
mv initialeq.00005 initialeq$i.00005
python compublic.py initialeq$i.00005 com_file.dat
mv com_file.dat comeq5$i.dat 

python bin2asciipublic.py initialeq$i.00005 initialeq$i.00005.dat
head -n 100000 initialeq$i.00005.dat >> haloo
tail -n 100000 initialeq$i.00005.dat >> bulgeo
rm initialeq$i.00005.dat

./adconv6gyrs 100.4029 -50.9661 102.0 $mualphanew $mudeltanew 223.1
mv convvel.dat orbx$i.dat 
python Car_orbit.py orbx$i.dat orbxpast$i.dat
python units2public.py 200000 100000 100000 orbxpast$i.dat comeq5$i.dat

./ascii2tipsy 1.0 1.0
mv haloo halooorb$i
mv halo halorb$i
mv bulgeo bulgeoorb$i
mv bulge bulgeorb$i
mpiexec  PKDGRAV sim1.param
echo 'i'$i initialeq$i.00005
mv initial.bin initial$i.bin
mv final.00060 final$i.00060
echo 'Im haaaaapppy'
python sbpublic.py final$i.00060 simb$i.dat simh$i.dat simb$i.png simh$i.png
python disppublic.py final$i.00060 simdisp$i.dat simdisp$i.png
python bin2asciipublic.py final$i.00060 final$i.00060.dat
tail -n 100000 final$i.00060.dat >> final$i.00060_stars.dat 
./papublic 100000 final$i.00060_stars.dat simpa$i.dat com_file.dat
python gradpublic.py simpa$i.dat  simgrad$i.dat 100000
python statpublic.py obsdispnoisy.dat obsbnoisy.dat simdisp$i.dat simb$i.dat simgrad$i.dat simchi$i.dat

rm final*dat
mv simb*dat PROFILES/.
mv simh*dat PROFILES/.
mv simdisp*dat PROFILES/.
mv simgrad*dat PROFILES/.
mv simpa*dat PROFILES/.
mv initial*5 INFILES/.
mv initial*bin INFILES/.
mv final*60 FINALFILES/.

batchtest=$[ $i % 50 ]
echo 'batchtest=' $batchtest 
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
	echo ' Acceptance ratio='$ratioacc >> Valuestokeep.txt
	echo 'nbatch='$nbatch >> Valuestokeep.txt
	echo 'Step change='$dstep  >> Valuestokeep.txt
	echo 'Old step sizes: dalpha='$dalpha 'ddelta='$ddelta 'dMs='$dMs 'dr1='$dr1 'dmhalo='$dmhalo>> Valuestokeep.txt
	dalpha=`echo "scale=10;e(((l($dalpha^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	ddelta=`echo "scale=10;e(((l($ddelta^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	dMs=`echo "scale=10;e(((l($dMs^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	dr1=`echo "scale=10;e(((l($dr1^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	drst=`echo "scale=10;e(((l($drst^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	dmhalo=`echo "scale=10;e(((l($dmhalo^2)/l(10)+$dstep)/2.0)*l(10.0))"|bc -l` 
	echo 'The step sizes changed in this one. 
              New step sizes: dalpha='$dalpha 'ddelta='$ddelta 'dMs='$dMs 'dr1='$dr1 'drst='$drst 'dmhalo='$dmhalo>> Valuestokeep.txt
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


