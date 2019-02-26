#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-08:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test   # Partition to submit to
#SBATCH --mem=10000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o conv_table.%A.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -J conv_table

module load gcc/7.1.0-fasrc01 R/3.5.0-fasrc02 eigen/3.3.4-fasrc03

Rscript --no-restore get_convergence_table.R  
