#!/bin/csh
# Calculate FILLET Benchmark 1 cases with HEXTOR (21 June 2024 JDH)
#---------------------------------------------------------------------------------------------------
#set echo

set expname     = "ben1"
set wdir 	= /models/hextor
set modeldir 	= $wdir/model
set filletdir   = $wdir/fillet/${expname}/
set outdir	= $modeldir/out
set outfile	= $outdir/fillet.out
set outglobal	= $outdir/fillet_global.out
set filletout   = $filletdir/global_output_HEXTOR_${expname}.dat
set namelist    = input.nml
set namelistB   = $filletdir/input.nml.${expname}
set count       = 0

cd $modeldir
make -f Makefile

#---------------------------------------------------------------------------------------------------

cat $filletdir/global_header.txt > $filletout
cp $filletdir/lat_header.txt $filletdir/lat_output_HEXTOR_${expname}.dat
cp $namelistB $namelist
echo ${count}: ${expname}

./driver > /dev/null

tail -n 1 $outglobal | sed '0,/0/{s/0/'${count}'/}' >> $filletout
cat $outfile >> $filletdir/lat_output_HEXTOR_${expname}.dat

