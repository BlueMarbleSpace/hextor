#!/bin/csh
# Calculate FILLET Experiment 4 cases with HEXTOR (22 June 2024 JDH)
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
set Tcold       = 233.0
set Twarm       = 288.0
set count       = 0

cd $modeldir
make -f Makefile

cat $filletdir/global_header.txt > $filletout

#---------------------------------------------------------------------------------------------------
# Start from a cold state and increase co2
#---------------------------------------------------------------------------------------------------

foreach fco2 ( 1.00000000e-06 1.26485522e-06 1.59985872e-06 2.02358965e-06 2.55954792e-06 3.23745754e-06 4.09491506e-06 5.17947468e-06 6.55128557e-06 8.28642773e-06 1.04811313e-05 1.32571137e-05 1.67683294e-05 2.12095089e-05 2.68269580e-05 3.39322177e-05 4.29193426e-05 5.42867544e-05 6.86648845e-05 8.68511374e-05 1.09854114e-04 1.38949549e-04 1.75751062e-04 2.22299648e-04 2.81176870e-04 3.55648031e-04 4.49843267e-04 5.68986603e-04 7.19685673e-04 9.10298178e-04 1.15139540e-03 1.45634848e-03 1.84206997e-03 2.32995181e-03 2.94705170e-03 3.72759372e-03 4.71486636e-03 5.96362332e-03 7.54312006e-03 9.54095476e-03 1.20679264e-02 1.52641797e-02 1.93069773e-02 2.44205309e-02 3.08884360e-02 3.90693994e-02 4.94171336e-02 6.25055193e-02 7.90604321e-02 1.00000000e-01 )

  set co2 = `python -c "print( '{:.1F}'.format( ${fco2} * 1e6 ) )"`

  cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat
  sed -i "s/XXX/${co2}/" $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat
  sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_cold_$count.dat

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
# Start from a warm state and decrease co2
#---------------------------------------------------------------------------------------------------

foreach fco2 ( 1.00000000e-01 7.90604321e-02 6.25055193e-02 4.94171336e-02 3.90693994e-02 3.08884360e-02 2.44205309e-02 1.93069773e-02 1.52641797e-02 1.20679264e-02 9.54095476e-03 7.54312006e-03 5.96362332e-03 4.71486636e-03 3.72759372e-03 2.94705170e-03 2.32995181e-03 1.84206997e-03 1.45634848e-03 1.15139540e-03 9.10298178e-04 7.19685673e-04 5.68986603e-04 4.49843267e-04 3.55648031e-04 2.81176870e-04 2.22299648e-04 1.75751062e-04 1.38949549e-04 1.09854114e-04 8.68511374e-05 6.86648845e-05 5.42867544e-05 4.29193426e-05 3.39322177e-05 2.68269580e-05 2.12095089e-05 1.67683294e-05 1.32571137e-05 1.04811313e-05 8.28642773e-06 6.55128557e-06 5.17947468e-06 4.09491506e-06 3.23745754e-06 2.55954792e-06 2.02358965e-06 1.59985872e-06 1.26485522e-06 1.00000000e-06 )

  set co2 = `python -c "print( '{:.1F}'.format( ${fco2} * 1e6 ) )"`

  cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat
  sed -i "s/XXX/${co2}/" $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat
  sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat

  cp $namelistB $namelist
  sed -i "s/XXX/${fco2}/" $namelist
  sed -i "s/TTT/${Twarm}/" $namelist

  echo ${count}: " co2 =" ${co2} "Warm Start"
  ./driver > /dev/null

  tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
  cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}_warm_$count.dat

  @ count++
end
