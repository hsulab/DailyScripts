#!/bin/bash 

# few variables
dimension=3 # optimize 2D (x,y) or 3D (x,y,z)

dat_file='bonds.dat'
geo_file='POSCAR'

ia=1
pz=0.364156360055
lz=18.87256432

# run
echo "Sample Selected Bond at a Given Interval" > print-out
if [ -d CONTCARs ] 
then 
    rm -r CONTCARs
    mkdir CONTCARs
    echo "# Length Energy" > ${dat_file}
else
    mkdir CONTCARs
    echo "# Length Energy" > ${dat_file}
fi

cp ${geo_file} ${geo_file}.org

# BM optimization 
for s in `seq 1.4 0.1 2.8`
do 
  cp ${geo_file} ${geo_file}.bak

  echo "Bond Length ${s}" >> print-out
  awk '{ if(NR=='$((ia + 9))') {printf "  %16.12f  %16.12f  %16.12f   %s   %s   %s\n" $1 $2 '${s}'/'${lz}'+'$pz' $4 $5 $6} else {print $0} }' ${geo_file}.bak > ${geo_file}

  cur_dis=$(python3 ./dist.py)

  echo $cur_dis

  # aprun

  break

  # get energy and back up
  energy=$(grep "sigma->" OUTCAR | tail -n 1)
  energy=${energy##*=}

  #suffix=$((s * 100)) # Warning! Change original variable.
  cp CONTCAR CONTCARs/CONTCAR_${s}
  printf "%8.4f  %12.8f\n" $s $energy >> ${dat_file}

  # break

done

cp ${geo_file}.org ${geo_file}

