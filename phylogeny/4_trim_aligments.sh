#!/bin/bash
#SBATCH -J trim                               # Job name
#SBATCH -o trim.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 0-10:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 80000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p general,interact        
#SBATCH --export gene

# This script calls a python script to trim alignments
# The alignment to be trimmed should be given as a command-line argument


export alignment=$gene/${gene}_aligned_final.fasta
#now trim
echo TRIMMING ALIGNMENTS FOR $gene


./trim_alignments.py -a $alignment -m 100 -c 0.2 -s 10 
