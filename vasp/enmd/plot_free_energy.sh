#!/bin/bash

gnuplot -persist << EOF

set term png
set output "fe.png"

set title "Slow-Growth Molecular Dynamics"

set ylabel "Free Energy / eV"

set xlabel "ksi (reaction coordinate)"

plot "fe.dat" u 1:2 w l t "Energy Gradient", \
     "fe.dat" u 1:3 w l t "Free Energy"

EOF
