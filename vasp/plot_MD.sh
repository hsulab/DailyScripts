#!/bin/bash
echo "# Time  Temp  Ener" > MD.dat
grep 'E0' OSZICAR | awk '{print "  "$1"  "$3"  "$9}' >> MD.dat

echo "Write MD Data to MD_data"

gnuplot -persist << EOF

set term png
set output "MD.png"

set title "Molecular Dynamics"

set ylabel "Temperature / K"
set y2label "Potential Energy / eV"

set y2tics
set ytics nomirror

set xlabel "Time / fs"

plot "MD.dat" u 1:2 w l axis x1y1 t "Temperature", \
     "MD.dat" u 1:3 w l axis x1y2 t "Potential Energy", \
     300 w l axis x1y1 t "Equilibrium Temperature"

EOF

echo "Plot Energy/Temperature v.s. Timestep to MD.png"
