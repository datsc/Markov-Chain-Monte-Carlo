This directory will soon have all the codes used in the Dwarfed Masses project to simulate and estimate the mass of a dwarf galaxy interacting with the Milky Way. The list is given below.

For the moment it just has the main bash script for MCMC, the statistics code calculating chi-squared and the velocity dispersion code.

The data can be accessed here: http://vo.aip.de/dwarfedmasses/

The publication of this project and scientific details can be found here:
http://www.nature.com/ncomms/2015/150702/ncomms8599/full/ncomms8599.html


## ANALYSIS CODES
####################
To only use the analysis codes on the provided binary files, you need to install the tipsylib within the orblib and have the following codes and the observed profile files for either the Carina or the Mock chains.

bin2asciipublic.py -- reading the binary file and writing out ascii.
sbpublic.py -- calculates the surface brightness and prints out the c.o.m used in the velocity codes.
disppublic.py   -- calculates the velocity dispersion.
papublic.py --  calculates position along major axis to be used in the gradient code.
gradpublic.py  -- calculates the velocity gradient.
statpublic.py -- calculates the chi squared.

An example for how to run them is:

python sbpublic.py fin0.00060 simb0.dat simh0.dat simb0.png simh0.png
python disppublic.py fin0.00060 simdisp0.dat simdisp0.png
python bin2asciipublic.py fin0.00060 fin0.00060.dat
tail -n 100000 fin0.00060.dat >> fin0.00060_stars.dat 
./papublic 100000 fin0.00060_stars.dat simpa0.dat com_file.dat
python gradpublic.py simpa0.dat simgrad0.dat 100000 
python statpublic.py obsdispnoisy.dat obsbnoisy.dat simdisp0.dat simb0.dat simgrad0.dat simchi0.dat

## Running MCMC
####################
The main script is mcmc_core_exam_public.sh

After getting all the software ready to run (orblib and tipsylib within,pkdgrav,NEMO and falcon within), you need to provide this script with your simulation directory name (the SIMDIR variable), pkdgrav directory name (PKDGRAV variable) and the directory for the library of the orbit code (PYTHONPATH=hannilib). The parameter "seed" need to be provided each time.

The example file initialparam.dat provides the MCMC script with the values of the parameters the first simulation of the chain will use. It has total stellar mass and the scalelength, total halo mass and scalelength and the proper motion. 

The codes below all need to be in the simulation directory, as well as the observed profiles obsdisp.dat and obsb.dat 

disppublic.py  
gradpublic.py
sbpublic.py
statpublic.py
rangenpublic.py
units1public.py
units2public.py
compublic.py
Car_orbit.py
bin2asciipublic.py
ascii2tipsy.cc (compiled to ascii2tipsy)
Adconv.c (compiled to adconv6gyrs)
papublic.c (compiled to papublic)
xdrfuncs.h
