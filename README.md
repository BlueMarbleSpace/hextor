# HEXTOR (Habitable EBM for eXoplaneT ObseRvations)

README file for release 4.0.2

## Obtaining the Code
Use the command below to clone the repository with submodules. This is required to obtain the ffhash module, which is required for the radiaitive transfer lookup table.

`git clone --recurse-submodules https://github.com/BlueMarbleSpace/hextor.git`

## Running the Model

1. Edit ./model/Makefile and set the values of the Fortran compiler FC and location of this model WDIR

2. Edit runEBM.sh set the value of "wdir"

3. Copy a desired namelist template file from ./namelists/ or ./fillet/ to the . directory. Adjust any model parameters by editing the file ./input.nml

4. Adjust the length of the model integration by setting the value of the "niter" parameter in ./model/driver.f

5. Adjust the star by uncommenting the appropriate line in radiation/radiation.f90

6. Run the script by executing ./runEBM.sh

7. Model output is located in ./model/out/

- The file "tempseries.out" lists columns as: 
    time, max temp, min temp, pg0, pco2, co2condensation, fh2, diffusion

- CO2 condensation is indicated by a 1 (condensing) or 0 (not condensing) in the file above.

- H2 is not implemented in the current version of the lookup table

