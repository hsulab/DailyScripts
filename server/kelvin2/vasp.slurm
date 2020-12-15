#!/bin/bash
#SBATCH --partition=medpri           # queue
#SBATCH --job-name=pt-fcc-md-2200    # Job name
#SBATCH --ntasks=20                  # Number of cores
#SBATCH --time=12:00:00              # Time limit hrs:min:sec
#SBATCH --output=output.txt          # Standard output and error log
#SBATCH --error=error.txt            # Standard output and error log

module purge
module load services/s3cmd
module load compilers/intel/2016u1 # not necessary
module load libs/intel/2016u1
module load libs/intel-mkl/2016u1/bin
module load mpi/intel-mpi/2016u1/bin

echo `date "+%Y-%m-%d %H:%M:%S"` `pwd` >> $HOME/submitted
mpirun -np 20 $HOME/apps/vasp/bin/5.4.1-TS/vasp_std 2>&1 > print-out
echo `date "+%Y-%m-%d %H:%M:%S"` `pwd` >> $HOME/finished
