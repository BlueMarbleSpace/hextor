#!/bin/csh
# Calculate FILLET Experiment 1a cases with HEXTOR (20 June 2024 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set expname     = "exp2a"
set wdir 	= /models/hextor
set modeldir 	= $wdir/model
set filletdir   = $wdir/fillet/${expname}/
set outdir	= $modeldir/out
set outfile	= $outdir/fillet.out
set outglobal	= $outdir/fillet_global.out
set filletout   = $filletdir/global_output_HEXTOR_${expname}.dat
set namelist    = input.nml
set namelistB   = $filletdir/input.nml.${expname}
set axismin      = 0.875
set axismax      = 1.1
set axisint      = -0.0125
set oblmin      = 0.0
set oblmax      = 90.0
set oblint      = 10.0
set count       = 0

cd $modeldir
make -f Makefile

#---------------------------------------------------------------------------------------------------
# Start from a warm state with zero obliquity, and vary the orbital distance and obliquity
#---------------------------------------------------------------------------------------------------

cat $filletdir/global_header.txt > $filletout

foreach axis ( `seq $axismax $axisint $axismin` )
  foreach obliquity ( `seq $oblmin $oblint $oblmax` )

    set axis_cm = `python -c "print( '{:.9E}'.format( ${axis} * 1.495978707E13 ) )"`
    set solflux = `python -c "print( '{:.0F}'.format( 1361 / ${axis}**2 ) )"`
    set solcon  = `python -c "print( '{:.3F}'.format( 1 / ${axis}**2 ) )"`

    cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}_$count.dat
    sed -i "s/XXX/${solcon}/" $filletdir/lat_output_HEXTOR_${expname}_$count.dat
    sed -i "s/QQQ/${obliquity}/" $filletdir/lat_output_HEXTOR_${expname}_$count.dat
    sed -i "s/CCC/${count}/" $filletdir/lat_output_HEXTOR_${expname}_$count.dat

    cp $namelistB $namelist
    sed -i "s/XXX/${solflux}/" $namelist
    sed -i "s/ZZZ/${axis_cm}/" $namelist
    sed -i "s/QQQ/${obliquity}/" $namelist

    echo ${count}: " axis =" ${axis} "  Obliquity =" ${obliquity}
    ./driver > /dev/null

    tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
    cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}_$count.dat

    @ count++
  end
end
