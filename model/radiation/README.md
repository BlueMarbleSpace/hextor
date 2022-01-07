# Radiaitive Transfer Hash Table for HEXTOR

README file for the radiaitive transfer lookup table module to accompany HEXTOR release 4.0.0

## Using this Module with HEXTOR

If you are using this code to run the HEXTOR climate model, then you may need to modify the header of radiation.f90 to select a different host star. You can uncomment/comment the appropriate .h5 database to toggle between different host stars. 

Current options are the Sun and a 2600 K BT-Settl star.

## Using this Module as an API

You can also use this radiaitve transfer module in standalone mode to calculate radiaitive transfer properties with the lookup table. You may also choose to couple this module to your own climate model.

The getOLR() and getPALB() routines serve as the API for interfacing with the lookup table.

The main.f90 program provides a way to call the lookup table in standalone mode. This file can serve as a template for using the lookup table with other models. A Makefile is provided for main.f90.

Note that you must edit the header of radiation.f90 to run the model in standalone model. See the comments in the header for details.

