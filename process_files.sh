#!/bin/bash

#I should have optional argument for the base name, rather than doing the more risky '*.freqs'

mv assemble.txt assemble.bak

for point_dir in `ls -d r_*` 
do
cd $point_dir
point=${point_dir#r_} #strip just to get the step point
#delh needs both E(elec) (actual) from relaxed scan, and ZPE from the freq calculation
eelec=`cat *.ediff`
zpe=`awk '/Non-thermal \(ZPE\) correction/ { print $4 * 2625.499638 }' *.out`
delh=`echo "$eelec + $zpe" | bc`
echo "ctst  'TS-${point//.}'  ${delh}     ${point}               !!!!!!!!!!!    <= loose TS" >> ../assemble.txt
echo "FORMULA" >> ../assemble.txt
echo "1. VTST" >> ../assemble.txt
echo "2. (blank comment line)" >> ../assemble.txt
echo "3. (blank comment line)" >> ../assemble.txt
echo "1   1   1" >> ../assemble.txt
echo "0.0   2" >> ../assemble.txt
cat *.freqs >> ../assemble.txt 
num_vibs=`awk 'NR==1{print $1}' *.freqs`
awk -v n_vibs=$num_vibs 'BEGIN { FS=":"; } /Bk/ { print n_vibs+1 "      kro    "$2"   1       1     !     K-rotor"  }' *.rotors >> ../assemble.txt
awk -v n_vibs=$num_vibs 'BEGIN { FS=":"; } /B2d/ { print n_vibs+2 "      jro    "$2"   1       2     !     2D-rotor"  }' *.rotors >> ../assemble.txt
echo "" >> ../assemble.txt
cd ..
done
