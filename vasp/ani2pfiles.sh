#!/bin/bash

# ====== ====== ====== ====== ====== 
# Filename: ani2pfiles.sh 
# Author: jyxu 

# Format of OUT.ANI looks like 
# 1 natoms
# 2 STEP = 1
# 3 symbol x y z 
# ......
# N-1 natoms 
# N STEP = N 
# N+1 symbol x y z 

# ====== ====== ====== ====== ====== 

# default settings 
ani="./OUT.ANI" # OUT.ANI to be read 
oripos="./POSCAR" # original poscar, POSCAR for getting structure info 
pdir="../pfiles" # directory for POSCARs 
# be careful with the different vasp format 
nlines=8 # poscar foramt, vasp4 format uses 8 while vasp5 uses 9

if [ $# -eq 0 ]
then 
    echo "ani2pfiles.sh [OUT.ANI] [POSCAR] [PFILES]"
elif [ $# -eq 1 ]
then 
    echo "Write POSCARs to a customized directory $1."
    pdir=$1
fi 

# set directory path for storing poscars 
if [ ! -f $ani ] 
then 
    echo "$ani does not exist"
    exit 
fi 

if [ ! -f $oripos ] 
then 
    echo "$oripos does not exist"
    exit 
fi 

if [ -d $pdir ] 
then 
    echo "Directory for POSCARs exists -> $pdir"
    echo "Better check $pdir before running this script"
else
    echo "Create a new directory for POSCARs -> $pdir"
    mkdir $pdir
fi

# get number of steps 
nsteps=`grep "STEP" $ani | tail -n 1 | wc -l`

# get number of atoms 
natoms=`head -n 1 $ani`
natoms=${natoms##* }

# get lattice info 
strinfo=`head -n $nlines $oripos`
a1=`sed -n '3p' $oripos | awk '{printf}'`

# read POSCAR to get fixes (T or F)
oldifs=$IFS
IFS=$'\n'
cors=$(`tail -n +$((nlines+1)) $oripos | head -n +$natoms`)
tfs=[]
count=1
for cor in ${cors[@]}
do 
    cor=${cor##*[0,1,2,3,4,5,6,7,8,9]}
    tfs[$count]=$cor
    (( count++ ))
done

# read OUT.ANI and write to POSCARs
for i in `seq 1 $nsteps`
do 
    # read OUT.ANI and get coordinates 
    (( firstline=(i-1)*(natoms+2)+3 ))
    (( lastline=i*(natoms+2) ))
    cors=`tail -n +$firstline $ani | head -n +$natoms`
    cors=($cors) 
    xyzs=[]
    count=1
    for cor in ${cors[@]} 
    do 
        cor=${cor:5}
        xyzs[$count]=$cor 
        (( count++ ))
    done 
    
    # write to poscar 
    sufx=$( printf "%04d" "$i" )
    posfile=$pdir/"p${sufx}"
    if [ -f $posfile ] 
    then 
        echo "$posfile exists, remove and write a new one"
        rm $posfile
    else 
        echo "$posfile is written"
    fi

    echo "$strinfo" >> $posfile 
    for j in `seq 1 $natoms` 
    do 
        echo "     ${xyzs[$j]}    ${tfs[$j]}" >> $posfile
    done
    if [ $i -eq 10 ]
    then 
        exit # for test 
    fi
done
IFS=$oldifs
