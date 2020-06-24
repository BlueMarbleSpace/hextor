#!/bin/csh
# Runscript for Darren Williams' EBM (2009.06.01 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set wdir 	= /models/hextor
set modeldir 	= $wdir/model
set plotdir 	= $wdir/plots
set outdir	= $modeldir/out
set outfile	= $outdir/model.out
set sep         = "----------------------------------------"

#---------------------------------------------------------------------------------------------------

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
# Uncomment for NCL plots
#
# echo $sep
# cd $plotdir
#
# ncl plotZonal.ncl
# ps2pdf plotZonal.ps && rm plotZonal.ps
