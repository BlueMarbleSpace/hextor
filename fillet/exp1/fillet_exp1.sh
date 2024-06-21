#!/bin/csh
# Calculate FILLET Experiment 1 cases with HEXTOR (20 June 2024 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set expname     = "exp1"
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
set solmax      = 1.25
set solint      = 0.025
set oblmin      = 0.0
set oblmax      = 90.0
set oblint      = 10.0
set count       = 0

cd $modeldir
make -f Makefile

#---------------------------------------------------------------------------------------------------
# Start from a warm state with zero obliquity, and increase the solar constant and obliquity
#---------------------------------------------------------------------------------------------------

cat $filletdir/global_header.txt > $filletout

foreach solcon ( `seq $solmin $solint $solmax` )
  foreach obliquity ( `seq $oblmin $oblint $oblmax` )

    cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_$count.dat
    sed -i "s/XXX/${solcon}/" $filletdir/lat_output_HEXTOR_${expname}_$count.dat
    sed -i "s/QQQ/${obliquity}/" $filletdir/lat_output_HEXTOR_${expname}_$count.dat
    sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_$count.dat

    cp $namelistB $namelist
    sed -i "s/XXX/${solcon}/" $namelist
    sed -i "s/QQQ/${obliquity}/" $namelist

    echo ${count}: " Solcon =" ${solcon} "  Obliquity =" ${obliquity}
    ./driver > /dev/null

    tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
    cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}_$count.dat

    @ count++
  end
end
