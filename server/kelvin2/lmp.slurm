#!/bin/bash -l 
#SBATCH --partition=k2-hipri         # queue
#SBATCH --job-name=bPtO2-nvt         # Job name
#SBATCH --time=00:30:00              # Time limit hrs:min:sec
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=5                   # Number of cores
#SBATCH --cpus-per-task=4            # Number of cores per MPI task 
#SBATCH --mem=10G                    # Number of cores per MPI task 
#SBATCH --output=slurm.o%j           # Standard output and error log
#SBATCH --error=slurm.e%j            # Standard output and error log

module purge
module load services/s3cmd
module load libs/intel/2016u1
module load libs/intel-mkl/2016u1/bin
module load mpi/intel-mpi/2016u1/bin

conda activate dpdev

srun --ntasks=1 --cpus-per-task=4 cd 100 && mpirun -n 4 lmp -in in.lammps 2>&1 > lmp.out && cd - &
srun --ntasks=1 --cpus-per-task=4 cd 200 && mpirun -n 4 lmp -in in.lammps 2>&1 > lmp.out && cd - &
srun --ntasks=1 --cpus-per-task=4 cd 300 && mpirun -n 4 lmp -in in.lammps 2>&1 > lmp.out && cd - &
srun --ntasks=1 --cpus-per-task=4 cd 400 && mpirun -n 4 lmp -in in.lammps 2>&1 > lmp.out && cd - &
srun --ntasks=1 --cpus-per-task=4 cd 500 && mpirun -n 4 lmp -in in.lammps 2>&1 > lmp.out && cd - &
wait
