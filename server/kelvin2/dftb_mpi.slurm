#!/bin/bash
#SBATCH --partition=hipri           # queue
#SBATCH --job-name=DFTBplus_mpi_Kelvin2 # Job name
#SBATCH --ntasks=20                  # Number of cores
#SBATCH --time=0:30:00               # Time limit hrs:min:sec
#SBATCH --output=output.txt          # Standard output and error log
#SBATCH --error=error.txt            # Standard output and error log

# remove all modules first and load all intel-related ones
module purge
module load services/s3cmd
module load compilers/intel/2016u1 # not necessary
module load libs/intel/2016u1
module load libs/intel-mkl/2016u1/bin
module load mpi/intel-mpi/2016u1/bin

# change VASP version here
EXEC="/mnt/scratch/chemistry-apps/jxu/dftbplus/installed/20.1/intel-2016/parallel/bin/dftb+"

echo `date "+%Y-%m-%d %H:%M:%S"` `pwd` >> $HOME/submitted
mpirun -np $EXEC 2>&1 > dftb.out
echo `date "+%Y-%m-%d %H:%M:%S"` `pwd` >> $HOME/finished
