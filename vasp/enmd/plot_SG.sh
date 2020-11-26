#!/bin/bash

# output setting
SG_data="SG.dat"

potim=$(grep "POTIM" INCAR | awk '{print $3}')
increm=$(grep "INCREM" INCAR | awk '{print $3}')

echo "Get STEP v.s. LAMBDA in REPORT"
echo "# step lambda POTIM=${potim} INCREM=${increm}" > $SG_data

#grep 'b_m' REPORT | cat -n | awk '{print $1*'$increm'"   "$6/$4}' >> $SG_data
grep 'b_m' REPORT | cat -n | awk '{print $1"   "$3}' >> $SG_data

gnuplot -persist << EOF

set term png
set output "SG.png"

set title "Slow-Growth Molecular Dynamics"

set ylabel "Free Energy Gradient / eV/Å"

set xlabel "ksi (reaction coordinate)"

plot "$SG_data" u 1:2 w l t "Free Energy Gradient", \
    0 w l t "Zero"
     

EOF
