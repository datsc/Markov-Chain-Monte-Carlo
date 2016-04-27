# this script is run in the final analysis directory. It combines the results
# of the 8 Markov Chains run for the Cuspy Dark Matter Halo Models.
# From each chain it first subtracts the first 100 simulations as these
# are necessary during the burn in of the Markov Chain. It uses the velocity
#files (vc), the final masses within the observed radius (finalmasses), orbits
#(orbpar). There were about 10000 simulations results to be organised in files
#that can be easily analysed statistically.

#!/bin/bash

i=1;
for((i==1;i<=9;i+=1));
do
cp Core2C$i/Acceptedmodels.dat AcceptedmodelsArtCusp_"$i".dat 
cp Core2C$i/chicombArtCusp"$i".dat  .
cp Core2C$i/textformyinitialparamArtCusp_$i.txt .
cat AcceptedmodelsArtCusp_"$i".dat >> AcceptedmodelsArtCuspComb.dat
if [ $i == 1 ]; then t=`echo 978|bc -l`  tn=`echo "${t}-100"|bc -l`
elif [ $i == 2 ]; then t=`echo 937|bc -l` tn=`echo "${t}-100"|bc -l`
elif [ $i == 3 ]; then t=`echo 1034|bc -l`  tn=`echo "${t}-100"|bc -l`
elif [ $i == 4 ]; then t=`echo 1133|bc -l`  tn=`echo "${t}-100"|bc -l`
elif [ $i == 5 ]; then t=`echo 800|bc -l` tn=`echo "${t}-100"|bc -l`
elif [ $i == 6 ]; then t=`echo 1390|bc -l`  tn=`echo "${t}-100"|bc -l`
elif [ $i == 7 ]; then t=`echo 1034|bc -l` tn=`echo "${t}-100"|bc -l`
elif [ $i == 8 ]; then t=`echo 874|bc -l` tn=`echo "${t}-100"|bc -l`
elif [ $i == 9 ]; then t=`echo 775|bc -l`  tn=`echo "${t}-100"|bc -l`
fi
head -n $t vc_artcusp${i}.dat|tail -n $tn >> vcArtCusp${i}_corr.dat
head -n $t vc_inartcusp${i}.dat|tail -n $tn >> vcinArtCusp${i}_corr.dat
head -n $t finalmasses_artcusp${i}.dat|tail -n $tn >> finalmassesArtCusp${i}_corr.dat
head -n $t vcnew_artcusp${i}.dat|tail -n $tn >> vcArtCuspnew${i}_corr.dat
head -n $t orbparAC${i}.dat|tail -n $tn >> orbparAC${i}_corr.dat
done


k=1
for((k==1;k<=9;k+=1));
do
cat vcArtCuspnew${k}_corr.dat >> vcArtCuspnew.dat
cat orbparAC${k}_corr.dat >> orbparArtCuspshort.dat
cat finalmassesArtCusp${k}_corr.dat >> finalmassesArtCusp.dat
cat vcArtCusp${k}_corr.dat >> vcArtCusp.dat
cat vcinArtCusp${k}_corr.dat >> vcinArtCusp.dat
done
