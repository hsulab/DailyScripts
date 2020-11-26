# run
dimension=3 # optimize 2D (x,y) or 3D (x,y,z)
fd_step=0.01 # finite difference step
EXEC=/home/e89/e89/peijun/bin/vasp.5.4.1-TS/vasp_std

echo "Calculate Force along Bond Distance with $fd_step Finite Difference" > print-out

echo '# Finite Difference Data' > fd.dat
echo '# File For_Dis Bak_Dis For_En Bak_En Force' >> fd.dat

cp POSCAR POSCAR.bak

# BM optimization 
for car in `ls CONTCARs/`
do 
  # read current distance
  cp ./CONTCARs/$car POSCAR
  cur_dis=$(python3 ./dist.py)
  #printf "cur %2.12f\n" $cur_dis
  for_dis=$(echo $cur_dis | awk '{print $1+'$fd_step'}')
  #printf "for %2.12f\n" $for_dis
  bak_dis=$(echo $cur_dis | awk '{print $1-'$fd_step'}')
  #printf "bak %2.12f\n" $bak_dis

  # forward
  scale=$(echo $for_dis $cur_dis | awk '{print $1/$2}')
  scale=$(printf "%2.12f\n" $scale)
  awk '{if(NR>=3 && NR<='$((dimension + 2))') {{for(i=1;i<=3;i++) printf "%14.8f",$i*'$scale'} {printf "\n"}} else print $0}' POSCAR.bak > POSCAR
  for_dis=$(python3 ./dist.py)
  echo "$car - Forward D" >> print-out
  aprun -n 24 $EXEC 2>&1 >> print-out

  for_en=$(grep "sigma->" OUTCAR | tail -n 1)
  for_en=${for_en##*=}
  echo " " >> print-out

  # backward
  scale=$(echo $bak_dis $cur_dis | awk '{print $1/$2}')
  scale=$(printf "%2.12f\n" $scale)
  awk '{if(NR>=3 && NR<='$((dimension + 2))') {{for(i=1;i<=3;i++) printf "%14.8f",$i*'$scale'} {printf "\n"}} else print $0}' POSCAR.bak > POSCAR
  bak_dis=$(python3 ./dist.py)
  echo "$car - Backward D" >> print-out
  aprun -n 24 $EXEC 2>&1 >> print-out

  bak_en=$(grep "sigma->" OUTCAR | tail -n 1)
  bak_en=${bak_en##*=}
  echo " " >> print-out

  rep_for=$(echo $for_dis $bak_dis $for_en $bak_en | awk '{print -($3-$4)/($1-$2)}')
  printf "%s  %2.12f  %2.12f  %2.12f  %2.12f  %2.12f\n" $car $for_dis $bak_dis $for_en $bak_en $rep_for >> fd.dat

  # break

done

cp POSCAR.bak POSCAR
