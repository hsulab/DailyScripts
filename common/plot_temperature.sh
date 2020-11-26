#!/bin/bash

dat=$1
echo $dat

gnuplot -persist << EOF

set term png
set output "md_temp.png"

set title "Temperature Change in Molecular Dynamics"

set ylabel "Temperature (K)"

set xlabel "Time Step"

plot "${dat}" u 1:2 w l axis x1y1 t "Temperature", \
     300 w l lw 2 axis x1y1 t 'Equilibrium Temperature'

EOF
