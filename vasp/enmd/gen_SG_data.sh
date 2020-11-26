#!/bin/bash
echo "Get STEP v.s. LAMBDA in REPORT"
echo "# step lambda" > SG.dat

increm=0.01
#grep 'b_m' REPORT | cat -n | awk '{print $1*'$increm'"   "$6/$4}' >> SG_data
grep 'b_m' REPORT | cat -n | awk '{print $1"   "$6/$4}' >> SG.dat
