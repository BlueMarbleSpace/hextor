#!/bin/csh
# Runscript for Stochastic EBM (14 Feb 2014 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

#---------------------------------------------------------------------------------------------------
# Edit this value of wdir to match your system

set wdir 	= /models/hextor
#---------------------------------------------------------------------------------------------------

set modeldir 	= $wdir/model
set plotdir 	= $wdir/plots
set outdir	= $modeldir/out
set outfile	= $outdir/model.out
set namelist    = input.nml
set sep         = "----------------------------------------"

#---------------------------------------------------------------------------------------------------

cp $wdir/$namelist $modeldir
cd $modeldir

# This was needed to fix an old error with ifort, no longer necessary
# echo null > \(
# echo null > AS_NEEDED

make -f Makefile

#---------------------------------------------------------------------------------------------------

echo $sep
./driver

head -n 12 $outfile && tail -n 30 $outfile | head -n 5
head -n 33 $outfile | tail -n 18 > $outdir/zonal.out
tail -n 18 $outfile > $outdir/geog.out


#---------------------------------------------------------------------------------------------------

echo $sep
cd $plotdir

# Uncomment for NCL plots
# ncl -Q plotZonal.ncl
# ncl -Q plotTempSeries.ncl
# ncl -Q plotDailyAvgTemp.ncl
# ncl -Q plotDailyLatTemp.ncl
# ncl -Q plotPowerSpectrum.ncl
