#!/bin/csh
# Calculate FILLET Experiment 3 cases with HEXTOR (21 June 2024 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set expname     = "exp3"
set wdir 	= /models/hextor
set modeldir 	= $wdir/model
set filletdir   = $wdir/fillet/${expname}/
set outdir	= $modeldir/out
set outfile	= $outdir/fillet.out
set outglobal	= $outdir/fillet_global.out
set filletout   = $filletdir/global_output_HEXTOR_${expname}.dat
set namelist    = input.nml
set namelistB   = $filletdir/input.nml.${expname}
set solmin      = 0.8
set solmax      = 1.5
set solint      = 0.0125
set solintW     = -0.0125
set Tcold       = 233.0
set Twarm       = 288.0
set count       = 0

cd $modeldir
make -f Makefile

cat $filletdir/global_header.txt > $filletout

#---------------------------------------------------------------------------------------------------
# Start from a cold state and increase the solar constant
#---------------------------------------------------------------------------------------------------

foreach solcon ( `seq $solmin $solint $solmax` )

  cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat
  sed -i "s/XXX/${solcon}/" $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat
  sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat

  cp $namelistB $namelist
  sed -i "s/XXX/${solcon}/" $namelist
  sed -i "s/TTT/${Tcold}/" $namelist

  echo ${count}: " Solcon =" ${solcon} "Cold Start"
  ./driver > /dev/null

  tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
  cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat

  @ count++
end

#---------------------------------------------------------------------------------------------------
# Start from a warm state and decrease the solar constant
#---------------------------------------------------------------------------------------------------

foreach solcon ( `seq $solmax $solintW $solmin` )

  cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat
  sed -i "s/XXX/${solcon}/" $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat
  sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat

  cp $namelistB $namelist
  sed -i "s/XXX/${solcon}/" $namelist
  sed -i "s/TTT/${Twarm}/" $namelist

  echo ${count}: " Solcon =" ${solcon} "Warm Start"
  ./driver > /dev/null

  tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
  cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat

  @ count++
end
