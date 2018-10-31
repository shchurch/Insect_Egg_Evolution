#!/bin/bash
#SBATCH -J egg_epa                               # Job name 
#SBATCH -o raxml_epa.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 32                                    # Number of cores/cpus
#SBATCH -t 3-00:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 20000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p general,unrestricted,interact                      # Partition general, serial_requeue, unrestricted, interact
#SBATCH --export gene


#environment variable gene has to be specified as 18S or 28S before calling slurm script
echo processing $gene

module purge
module load gcc/6.2.0-fasrc01 python/2.7.8-fasrc01 mafft/7.245-fasrc01
source activate PYTHON 

cd raxml_epa

#align query to reference sequences
mafft --nuc --anysymbol --keeplength --thread $SLURM_NTASKS --addlong ../${gene}/query_unaligned.fasta --auto ${gene}_reference.fasta > ${gene}_epa_alignment.fasta

#create a phylip file suitable for raxml
python convert_fasta_to_phylip.py $gene

#now run raxml
rm RAxML*$gene*
raxmlHPC-PTHREADS-SSE3 --epa-keep-placements=1000 --epa-prob-threshold=0.00001 -T $SLURM_NTASKS -f v -s ${gene}_epa_alignment.phy -t misof_reference.tre -m GTRGAMMA -n $gene
