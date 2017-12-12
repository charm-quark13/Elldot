#! /bin/bash

function text1()
{

cat<<EOF
title "C2H4-C2F4 Dimer Optimization"
start canon-dimer
echo

geometry units bohr
  C        0.0000000000      1.2651715494      5.6691779657
  C        0.0000000000     -1.2651715494      5.6691779657
  H       -1.6146914887      2.3323319404      5.6691779657
  H        1.6146914887      2.3323319404      5.6691779657
  H       -1.6146914887     -2.3323319404      5.6691779657
  H        1.6146914887     -2.3323319404      5.6691779657
  C        0.0000000000     -1.2387153855     -1.8897259886
  C        0.0000000000      1.2387153855     -1.8897259886
  F        2.0713286561     -2.6253963159     -1.8897259886
  F       -2.0713286561     -2.6253963159     -1.8897259886
  F       -2.0713286561      2.6253963159     -1.8897259886
  F        2.0713286561      2.6253963159     -1.8897259886
end

basis
  * library 6-31g*
end

dft
  xc xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.00
     cam 0.33 cam_alpha 0.19 cam_beta 0.46
  direct
  vectors output gs-dimer.movecs
end

dplot
  title Dimer Localization
  vectors dimer.movecs
  LimitXYZ
  -6.0 6.0 100
  -6.0 6.0 100
  -4.0 8.0 100
  gaussian
  orbitals view
    1
    $1
  output mo-$1.cube
end

task dft
task dplot

EOF
}

x=$1
y=$2

for i in $(seq $x $y);
do
echo "$(text1 $i)" > dimer
nwchem dimer
done

rm dimer

outfile="$3.cube"

echo "$(cat mo-$x.cube)" > "$outfile"

OLDIFS=$IFS
IFS=' '
read -ra zeta <<< $(sed "3q;d" mo-$x.cube)
skip=$((zeta+7))
#echo $skip
#echo $zeta
IFS=$OLDIFS

for i in $(seq $(($x+1)) $y)
do
#echo "### $i" >> "$outfile"
echo "$(tail --lines=+$skip mo-$i.cube)" >> "$outfile"
done

rm mo-*.cube

echo "Don't forget to add number of mo's to the first line!"
echo "Insert list of mo's included in the cube file after final atom coordinates!"
echo "MO list should look like: # of MO's : 1st-MO : 2nd-MO : ... : 9th-MO /n"

