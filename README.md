# README file for HEXTOR (Habitable EBM for eXoplaneT ObseRvations)
# Release 1.2.1

1) Edit ./model/Makefile and set the values of the Fortran compiler FC and location of this model WDIR

2) Edit runEBM.sh set the value of "wdir"

3) Check that ncl is installed on your system (optional, needed for plots only)

4) Adjust any model parameters by editing the file ./input.nml. 

5) Run the script by executing ./runEBM.sh

6) Model output is located in ./model/out/, and optional plots are *.eps files in ./plots/

