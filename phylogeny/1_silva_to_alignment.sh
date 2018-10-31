#!/bin/bash

#SBATCH -J slv_algn                               # Job name 
#SBATCH -o silva_to_alignment.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 7-00:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 500                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p general                       # Partition general, serial_requeue, unrestricted, interact
#SBATCH --export gene,email

#this script call silva_to_alignment.py, to generate a fasta file with sequences to be used for adding to reference alignment and then phylogenetic placement in raxml
#environment variable gene has to be specified as 18S or 28S before calling slurm script

module purge
module load python/2.7.8-fasrc01
source activate PYTHON
module load java/1.8.0_45-fasrc01

echo processing $gene
python silva_to_alignment.py $gene -e $email
