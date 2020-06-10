#!/bin/bash

# -----
# Author: Jiayan Xu, Queen's University Belfast
# Version: 0.0.1
# Description:
#   check SCF convergence of VASP
# -----


#
if [ $# -eq 1 ]
then 
    echo "Set nscf tp $1"
    echo ""
    nscf=$1
else
    echo "Set nscf default 120."
    echo ""
    nscf=120
fi

# get all scf steps
grep -B 1 "E0" print-out > tmp_scf

# get wrong steps
wrong_scf=`awk 'NR%3==1 {if ($2 >= '$nscf') print int(NR/3)+1}' tmp_scf`

# print-out some info
for scf in $wrong_scf
do 
    echo "SCF STEP $scf MAYBE WRONG!!!"
    grep -B 3 " $scf T=.*" print-out
    echo ""
done

rm tmp_scf
