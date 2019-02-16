#!/bin/bash
#SBATCH -J mrbayes                               # Job name 
#SBATCH -o mrbayes.%a.%A.out                     # File to which stdout will be written
#SBATCH -n 8                                    # number of ranks
#SBATCH -c 1                                    #number of threads
#SBATCH --contiguous
#SBATCH -t 7-00:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem-per-cpu 3000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared                     # Partition general, serial_requeue, unrestricted, interact
#SBATCH --export=level #may be useful for species level analyses in the future
#SBATCH --mail-type=ALL
#SBATCH --mail-user=souzademedeiros@fas.harvard.edu

#echo level $level

source new-modules.sh

module purge
#module load gcc/4.8.2-fasrc01 openmpi/1.10.0-fasrc01 mrbayes/3.2.5-fasrc01
#module load gcc/4.8.2-fasrc01 openmpi/1.10.0-fasrc01 mrbayes/3.2.6-fasrc01
module load intel/15.0.0-fasrc01 openmpi/1.10.0-fasrc01 mrbayes/3.2.6-fasrc01
export PATH
export LD_LIBRARY_PATH 


infile=$(find . -maxdepth 1 -name '*.nex' | sort | sed -n ${SLURM_ARRAY_TASK_ID}p)

srun -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK --mpi=pmi2 --export=ALL mb $infile
#mpirun -np 8 mb mrbayes_input.nex

