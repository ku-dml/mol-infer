#!/bin/bash

for f in ./*.sdf
do
  outf=$(echo $f | cut -d'_' -f 1)
  outfile="${outf}_fv4+.csv"
  echo $outf
  echo $outfile
  ./fv4 $f $outfile
done
