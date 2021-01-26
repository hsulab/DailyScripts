#!/bin/bash -l

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request time of wallclocl time (format hours:minutes:seconds).
#$ -l h_rt=3:00:00

# 3. Request 1 gigabyte of RAM per process (nust be an integer).
#$ -l mem=1G

# 4. Request 10 gigabyte of TMPDIR space per node (default is 10 GB).
#$ -l tmpfs=10G

# 5. Set the name of the job.
#$ -N CdS_cell_PBE

# 6. Select the MPI parallel environment and number of processes.
#$ -pe mpi 96

# 7. Set the working directory.
#$ -cwd

#add to QUB_C
#$ -P Gold
#$ -A QUB_chem

# 8. Run the MPI job, gerun is a wrapper on Thomas.
dimension=3 # optimize 2D (x,y) or 3D (x,y,z)

echo "${dimension}D Lattice Optimization using BM Equation" > print-out
if [ -d CONTCARs ] 
then 
    rm -r CONTCARs
else
    mkdir CONTCARs
    echo "${dimension}D Lattice Optimization using BM Equation" > fitting_data
fi

cp POSCAR POSCAR_origin

# BM optimization 
for s in `seq 0.9 0.01 1.1`
do 
    # change lattice in POSCAR
    scale=$(printf "%2.12f\n" $s)
    awk '{if(NR>=3 && NR<='$((dimension + 2))') {{for(i=1;i<=3;i++) printf "%14.8f",$i*'$scale'} {printf "\n"}} else print $0}' POSCAR_origin > POSCAR
    echo "Ions Optimization with Scaled $scale Lattice" >> print-out

    gerun $VASPPATH 2>&1 >> print-out

    # get energy and back up
    energy=$(grep "sigma->" OUTCAR | tail -n 1)
    energy=${energy##*=}

    #suffix=$((s * 100)) # Warning! Change original variable.
    suffix=$s # Warning! Change original variable.
    cp CONTCAR CONTCARs/CONTCAR_$suffix
    printf "%12f    %12f\n" $scale $energy >> fitting_data
done

cp POSCAR_origin POSCAR

BM_FIT.py >> print-out

echo `date "+%Y-%m-%d %H:%M:%S"` `pwd` >> $HOME/finished
