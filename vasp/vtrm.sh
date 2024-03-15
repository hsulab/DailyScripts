#!/bin/bash
# ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
# Author: Jiayan Xu
# Description:
#   Remove redundant output files.
# ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

# default saved files
saved_fnames=('INCAR' 'POSCAR' 'KPOINTS' 'POTCAR' 'CONTCAR' \
              'WAVECAR' 'CHG' 'CHGCAR' 'ICONST' 'PENALTYPOT' \
              'vdw_kernel.bindat' \
              'POSCAR_origin' 'RCOCAR' 'PCVCAR' \
              'plumed.inp' 'plumed.dat' 'plumed.out' \
              'vasp.pbs' 'vasp.slurm' 'vasp.sh' 'vasp.script' \
              'fort.188' 'ase-sort.dat')

# add saved files from arguments
saved_fnames=(${saved_fnames[@]} $@)

# out some infomation
echo "Remove redundant files, leaving directories unaffected."
echo ${saved_fnames[@]}

# remove files
for file in `ls`
do 
  # check if file in saved list
  removed_flag=0
  for fname in ${saved_fnames[@]}
  do 
    if [ $file == $fname ]
    then 
      removed_flag=0
      break
    else
      if [ ! -d $file ]
      then
         removed_flag=1
      fi
    fi
  done
  # remove marked files
  if [ $removed_flag == 1 ]
  then
    echo "remove $file"
    rm $file
  fi
done
