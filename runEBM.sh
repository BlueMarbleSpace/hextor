#!/bin/csh
# Runscript for EBM (02 August 2023 JHM)
#---------------------------------------------------------------------------------------------------
#set echo

#---------------------------------------------------------------------------------------------------
# Edit this value of wdir to match your system
set wdir 	= /models/hextor

# Update the symbolic link in the config folder to match your system configuration
$wdir/config/machine.sh

#---------------------------------------------------------------------------------------------------
# Do not edit below this line

set modeldir 	= $wdir/model
set plotdir 	= $wdir/plots
set outdir	= $modeldir/out
set outfile	= $outdir/model.out
set namelist    = input.nml
set sep         = "----------------------------------------"

#---------------------------------------------------------------------------------------------------

cp $wdir/$namelist $modeldir
cd $modeldir

make -f Makefile

#---------------------------------------------------------------------------------------------------

echo $sep
./driver

head -n 12 $outfile && tail -n 30 $outfile | head -n 5
head -n 35 $outfile | tail -n 18 > $outdir/zonal.out
tail -n 18 $outfile > $outdir/geog.out
