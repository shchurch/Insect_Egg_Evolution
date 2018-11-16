#!/bin/bash
#SBATCH -J stitch                               # Job name 
#SBATCH -o stitch.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # number of nodes
#SBATCH -n 1                               # Number of cores/cpus
#SBATCH -t 3-00:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 50000                               # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p interact,general,unrestricted                     # Partition general, serial_requeue, unrestricted, interact
#SBATCH --export=level    #species or genus

module purge
module load python/2.7.8-fasrc01
source activate PYTHON

cd ./mrbayes_smalltrees

cd rainford
echo STITCHING RAINFORD
python -u ../../stitch_mrbayes.py --threads $SLURM_NTASKS --aligndata ../../dataset_genus_overlap.csv --n-posterior 5 --backbone ../rainford_backbone.nexml --trees *.t


cd ../misof
#echo STITCHING MISOF
python -u ../../stitch_mrbayes.py --threads $SLURM_NTASKS  --aligndata ../../dataset_genus_overlap.csv --n-posterior 5 --backbone ../misof_backbone.nexml --trees *.t
