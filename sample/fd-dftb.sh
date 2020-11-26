# DFTB+

# run
dimension=3 # optimize 2D (x,y) or 3D (x,y,z)
fd_step=0.01 # finite difference step

# MPI SETTING
export KMP_AFFINITY=disabled # for intel compiled code

# Set the number of threads to 12
export OMP_NUM_THREADS=12 # the maximum for a single NUMA region

EXEC=/work/e89/e89/jyxu/apps/dftbplus/19.1/intel-2017/_install/bin/dftb+

echo "Calculate Force along Bond Distance with $fd_step Finite Difference" > print-out

echo '# Finite Difference Data' > fd.dat
echo '# File For_Dis Bak_Dis For_En Bak_En Force' >> fd.dat

CPREFIX='geo_in.gen'
cp $CPREFIX $CPREFIX.org

# BM optimization 
for car in `ls ./GENs/`
do 
  # read current distance
  cp ./GENs/$car ${CPREFIX}
  cp ${CPREFIX} ${CPREFIX}.bak
  cur_dis=$(python3 ./dist-geo.py)
  printf "cur %2.12f\n" $cur_dis
  for_dis=$(echo $cur_dis | awk '{print $1+'$fd_step'}')
  printf "for %2.12f\n" $for_dis
  bak_dis=$(echo $cur_dis | awk '{print $1-'$fd_step'}')
  printf "bak %2.12f\n" $bak_dis

  # break

  # forward
  scale=$(echo $for_dis $cur_dis | awk '{print $1/$2}')
  scale=$(printf "%2.12f\n" $scale)
  awk '{if(NR>=6 && NR<='$((dimension + 5))') {{for(i=1;i<=3;i++) printf "%14.8f",$i*'$scale'} {printf "\n"}} else print $0}' ${CPREFIX}.bak > ${CPREFIX}
  for_dis=$(python3 ./dist-geo.py)
  echo "$car - Forward D" >> print-out
  aprun -cc none -n 2 -N 2 -S 1 -d 12 $EXEC 2>&1 >> print-out

  for_en=$(grep "Extrapolated to 0" detailed.out | awk '{print $6}')
  echo " " >> print-out

  # backward
  scale=$(echo $bak_dis $cur_dis | awk '{print $1/$2}')
  scale=$(printf "%2.12f\n" $scale)
  awk '{if(NR>=6 && NR<='$((dimension + 5))') {{for(i=1;i<=3;i++) printf "%14.8f",$i*'$scale'} {printf "\n"}} else print $0}' ${CPREFIX}.bak > ${CPREFIX}
  bak_dis=$(python3 ./dist-geo.py)
  echo "$car - Backward D" >> print-out
  aprun -cc none -n 2 -N 2 -S 1 -d 12 $EXEC 2>&1 >> print-out

  bak_en=$(grep "Extrapolated to 0" detailed.out | awk '{print $6}')
  echo " " >> print-out

  rep_for=$(echo $for_dis $bak_dis $for_en $bak_en | awk '{print -($3-$4)/($1-$2)}')
  printf "%s  %2.12f  %2.12f  %2.12f  %2.12f  %2.12f\n" $car $for_dis $bak_dis $for_en $bak_en $rep_for >> fd.dat

  break

done

cp ${CPREFIX}.org ${CPREFIX}
