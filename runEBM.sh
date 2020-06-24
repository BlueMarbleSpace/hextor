#!/bin/csh
# Runscript for Darren Williams' EBM (2009.06.01 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set wdir 	= /neptune/s1/jhaqqmis/williams_ebm
set modeldir 	= $wdir/model
set outdir	= $modeldir/out
set outfile	= $outdir/model.out

cd $modeldir

echo null > \(
echo null > AS_NEEDED

make -f Makefile
./driver

#---------------------------------------------------------------------------------------------------

head -n 12 $outfile && tail -n 30 $outfile | head -n 5
