#!/bin/csh
# Calculate bifurcation curves for HEXTOR (06 October 2020 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set wdir 	= /models/hextor
set modeldir 	= $wdir/model
set plotdir 	= $wdir/plots
set outdir	= $modeldir/out
set outfile	= $outdir/model.out
set namelist    = input.nml
set namelistB   = $wdir/input.nml.bifurcation
set icelines    = $outdir/icelines.out
set bifurcCold  = $outdir/bifurcationCold.out
set bifurcWarm  = $outdir/bifurcationWarm.out
set bifurcHot   = $outdir/bifurcationHot.out
set bifurcBack  = $outdir/bifurcationBack.out
set solmin      = 0.6
set solmax      = 1.4
set solintC     = 0.02
set solintW     = -0.02
#set solintC     = 0.01
#set solintW     = -0.01
set sep         = "----------------------------------------"

cd $modeldir
make -f Makefile

#---------------------------------------------------------------------------------------------------
# First start from a cold state and increase the solar constant
#---------------------------------------------------------------------------------------------------

cat /dev/null > $bifurcCold

foreach solcon ( `seq $solmin $solintC $solmax` )

  cp $namelistB $namelist
  sed -i "s/XXX/${solcon}/" $namelist
  sed -i "s/RRR/0/" $namelist
  sed -i "s/TTT/233.0/" $namelist

  echo $sep
  ./driver

  cat $icelines >> $bifurcCold
  head -n 5 $outfile | tail -n 1

end


#---------------------------------------------------------------------------------------------------
# Then start from a present Earth state and increase the solar constant
#---------------------------------------------------------------------------------------------------

cat /dev/null > $bifurcBack

foreach solcon ( `seq $solmin $solintC $solmax` )

  cp $namelistB $namelist
  sed -i "s/XXX/${solcon}/" $namelist
  sed -i "s/RRR/0/" $namelist
  sed -i "s/TTT/273.0/" $namelist

  echo $sep
  ./driver

  cat $icelines >> $bifurcBack
  head -n 5 $outfile | tail -n 1 

end

#---------------------------------------------------------------------------------------------------
# Next start from a warm present-Earth state and decrease the solar constant
#---------------------------------------------------------------------------------------------------

cat /dev/null > $bifurcWarm

foreach solcon ( `seq $solmax $solintW $solmin` )

  cp $namelistB $namelist
  sed -i "s/XXX/${solcon}/" $namelist
  sed -i "s/RRR/0/" $namelist
  sed -i "s/TTT/273.0/" $namelist

  echo $sep
  ./driver

  cat $icelines >> $bifurcWarm
  head -n 5 $outfile | tail -n 1 

end


#---------------------------------------------------------------------------------------------------
# Finally, start from a hot state and decrease the solar constant
#---------------------------------------------------------------------------------------------------

cat /dev/null > $bifurcHot

foreach solcon ( `seq $solmax $solintW $solmin` )

  cp $namelistB $namelist
  sed -i "s/XXX/${solcon}/" $namelist
  sed -i "s/RRR/0/" $namelist
  sed -i "s/TTT/300.0/" $namelist

  echo $sep
  ./driver

  cat $icelines >> $bifurcHot
  head -n 5 $outfile | tail -n 1 

end


#---------------------------------------------------------------------------------------------------

head -n 33 $outfile | tail -n 18 > $outdir/zonal.out
tail -n 18 $outfile > $outdir/geog.out

echo $sep
cd $plotdir
ncl plotBistability.ncl
