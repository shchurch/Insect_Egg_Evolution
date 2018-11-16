#!/bin/bash
#SBATCH -J build                               # Job name
#SBATCH -o build.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 0-10:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 20000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p general,interact        
#SBATCH --export gene
#SBATCH --mail-type=END
#SBATCH --mail-user=souzademedeiros@fas.harvard.edu

#the following makes script stop on first error
set -e

# This script calls a python script to trim alignments
# The alignment to be trimmed should be given as a command-line argument


echo CONCATENATING AND BUILDING RAXML INPUT
python -u ./get_concat_alignment_for_phenotypes.py


echo BUILDING MRBAYES INPUTS FOR GENERA
python -u build_mrbayes.py

echo DONE
