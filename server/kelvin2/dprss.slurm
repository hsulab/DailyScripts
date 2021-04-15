#!/bin/bash
#SBATCH --partition=k2-medpri        # queue
#SBATCH --job-name=dprss-iter-001    # Job name
#SBATCH --nodes=1-1
#SBATCH --ntasks=1                   # Number of cores
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=slurm-%j.txt        # Standard output and error log
#SBATCH --error=slurm-%j.err         # Standard output and error log

export OMP_NUM_THREADS=4

source activate tf2 

echo `date "+%Y-%m-%d %H:%M:%S"` `pwd` >> $HOME/submitted
dprss dyn -i ./random_structures.xyz -w ./dyn 2>&1 > dprss.out
echo `date "+%Y-%m-%d %H:%M:%S"` `pwd` >> $HOME/finished
