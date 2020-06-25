# HEXTOR (Habitable EBM for eXoplaneT ObseRvations)

README file for release 3.5.0

1) Edit ./model/Makefile and set the values of the Fortran compiler FC and location of this model WDIR

2) Edit runEBM.sh set the value of "wdir"

3) Check that ncl is installed on your system (optional, needed for plots only)

4) Adjust any model parameters by editing the file ./input.nml. 

5) Adjust the length of the model integration by setting the value of the "niter" parameter in ./model/driver.f.

6) Run the script by executing ./runEBM.sh

7) Model output is located in ./model/out/, and optional plots are *.eps files in ./plots/

-- The file "tempseries.out" lists columns as: 
    time, temperature, pg0, pco2, co2condensation, co2liquid, fh2, diffusion

-- CO2 condensation is indicated by a 1 (condensing) or 0 (not condensing) in the file above.

