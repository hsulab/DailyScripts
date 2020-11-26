#!/bin/bash
MD_data=MD.dat

echo "# Time  Temp  Ener" > $MD_data
grep 'E0' print-out | awk '{print "  "$1"  "$3"  "$9}' >> $MD_data

echo "Write MD Data to $MD_data"

gnuplot -persist << EOF

set term png
set output "MD.png"

set title "Molecular Dynamics"

set ylabel "Temperature / K"
set y2label "Potential Energy / eV"

set y2tics
set ytics nomirror

set xlabel "Time / fs"

plot "$MD_data" u 1:2 w l axis x1y1 t "Temperature", \
     "$MD_data" u 1:3 w l axis x1y2 t "Potential Energy", \
     300 w l axis x1y1 t "Equilibrium Temperature"

EOF

echo "Plot Energy/Temperature v.s. Timestep to MD.png"
