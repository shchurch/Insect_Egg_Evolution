#!/bin/bash
#SBATCH -J stitch                               # Job name 
#SBATCH -o stitch.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # number of nodes
#SBATCH -n 1                               # Number of cores/cpus
#SBATCH -t 3-00:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 100000                               # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p serial_requeue                     # Partition general, serial_requeue, unrestricted, interact
#SBATCH --export=level    #species or genus

module purge
module load python/2.7.14-fasrc02
source activate general



cd mrbayes_smalltrees/treelinks_misof
echo STITCHING MISOF
python -u ../../stitch_mrbayes.py --threads $SLURM_NTASKS  --aligndata ../../dataset_genus_overlap.csv --n-posterior 100 --backbone ../misof_backbone.nexml --trees *.t



cd ../treelinks_rainford
echo STITCHING RAINFORD
python -u ../../stitch_mrbayes.py --threads $SLURM_NTASKS --aligndata ../../dataset_genus_overlap.csv --n-posterior 100 --backbone ../rainford_backbone.nexml --trees *.t
