#!/bin/bash
#SBATCH -J splitfreq                              # Job name 
#SBATCH -o splitfreq.%A.out                     # File to which stdout will be written
#SBATCH -n 32                               # Number of cores/cpus
#SBATCH -N 1                               # Number of cores/cpus
#SBATCH -t 0-08:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 120000                               # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p test                     # Partition general, serial_requeue, unrestricted, interact

module purge
module load python gcc/7.1.0-fasrc01 openmpi/3.1.3-fasrc01

source activate egg_treestats

ipcluster start --cluster-id=IP$SLURM_JOB_ID --n=$SLURM_NTASKS --ip=* &

sleep 60

python get_split_freqs.py --profile default --id IP$SLURM_JOB_ID 

ipcluster stop
