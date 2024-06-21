#!/bin/csh
# Calculate FILLET Experiment 4 cases with HEXTOR (21 June 2024 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set expname     = "exp4"
set wdir 	= /models/hextor
set modeldir 	= $wdir/model
set filletdir   = $wdir/fillet/${expname}/
set outdir	= $modeldir/out
set outfile	= $outdir/fillet.out
set outglobal	= $outdir/fillet_global.out
set filletout   = $filletdir/global_output_HEXTOR_${expname}.dat
set namelist    = input.nml
set namelistB   = $filletdir/input.nml.${expname}
set co2min      = 0       # min value is 1, adjusted below
set co2max      = 100000
set co2int      = 100
set co2intW     = -100
set Tcold       = 233.0
set Twarm       = 288.0
set count       = 0

cd $modeldir
make -f Makefile

cat $filletdir/global_header.txt > $filletout

#---------------------------------------------------------------------------------------------------
# Start from a cold state and increase the solar constant
#---------------------------------------------------------------------------------------------------

foreach co2 ( `seq $co2min $co2int $co2max` )

  if ( ${co2} == 0 ) then
    @ co2++
  endif

  cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat
  sed -i "s/XXX/${co2}/" $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat
  sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat

  set fco2 = `python -c "print( '{:.2E}'.format( ${co2} / 1e6 ) )"`

  cp $namelistB $namelist
  sed -i "s/XXX/${fco2}/" $namelist
  sed -i "s/TTT/${Tcold}/" $namelist

  echo ${count}: " co2 =" ${co2} "Cold Start"
  ./driver > /dev/null

  tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
  cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat

  @ count++
end

#---------------------------------------------------------------------------------------------------
# Start from a warm state and decrease the solar constant
#---------------------------------------------------------------------------------------------------

foreach co2 ( `seq $co2max $co2intW $co2min` )

  if ( ${co2} == 0 ) then
    @ co2++
  endif

  cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat
  sed -i "s/XXX/${co2}/" $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat
  sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat

  set fco2 = `python -c "print( '{:.2E}'.format( ${co2} / 1e6 ) )"`

  cp $namelistB $namelist
  sed -i "s/XXX/${fco2}/" $namelist
  sed -i "s/TTT/${Twarm}/" $namelist

  echo ${count}: " co2 =" ${co2} "Warm Start"
  ./driver > /dev/null

  tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
  cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat

  @ count++
end
