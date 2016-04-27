## This script read the orbit output of one simulation in the Markov chain
## and runs a different version of the orbit code that outputs every 100Myrs
## for 6Gyrs. It was used to visualise the orbits of the good models after the
## chains were completed.


#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:59:59
#PBS -l pvmem=1gb
#PBS -q devel
#PBS -m bea
#PBS -M uu2@le.ac.uk
#-o  /scratch/MODEL146/$JOB_NAME.o$JOB_ID
#-e  /scratch/MODEL146/$JOB_NAME.e$JOB_ID
#PBS -v PYTHONPATH=/home/dc-ural1/hannilib/
#PBS -S /bin/bash

#======================================
# Set up environment
#======================================
module load openmpi
module load gcc
module load python/2.7.3
module load swig
#======================================
cd /scratch/MODEL146

#!/bin/bash
cp orbx0.dat orbitin1.dat
echo 'hi'
for ((i==1;i<=100;i+=1));
do
echo $i
new=$i"00"
eval 'sed 's/6000/$new/g' orbitin1.dat > orf$i'
done

for ((j==1;j<=100;j+=1));
do
fileout='orout'$j
echo $fileout
python Ugur_orbit1.py orf$j $fileout
done

for ((k==1;k<=100;k+=1));
do 
filein='orout'$k
head $filein >>allorb.dat 
done
