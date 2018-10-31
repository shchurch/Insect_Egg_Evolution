#!/bin/bash
#SBATCH -J fast_align                               # Job name 
#SBATCH -o fast_align_final.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 32                                    # Number of cores/cpus
#SBATCH -t 3-00:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 80000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p general,unrestricted,interact                    # Partition general, serial_requeue, unrestricted, interact
#SBATCH --export gene


#environment variable gene has to be specified as 18S or 28S before calling slurm script
echo processing $gene using $SLURM_NTASKS cores

source deactivate
module purge
module load python/2.7.8-fasrc01 
source activate PYTHON 

#get sequences to align
python get_final_dataset.py $gene

#align sequences
cd $gene

#rm *shortnames*
#rm pasta*
rm -r tmp_upp
rm *shortnames*
#first, simplify names for upp
cat ${gene}_fullseqs_final.fasta | cut -d '|' -f 1 | cut -d '.' -f 1 > ${gene}_fullseqs_final_shortnames.fasta

#run upp
run_upp.py -B 700 -l 5 -x $SLURM_NTASKS -T 0.20 -M 2000 -s ${gene}_fullseqs_final_shortnames.fasta -o ${gene}_shortnames -p $PWD/tmp_upp/

#use code commented below if pasta was succesfully computed 
#run_upp.py -l 5 -x $SLURM_NTASKS -t pasta.fasttree -a pasta.fasta -s ${gene}_fullseqs_final_shortnames.fasta -o ${gene}_shortnames -p $PWD/tmp_upp/
#-B 500 the default is 1000, but sometimes this causes problems with pasta and a smaller number has more chances of success
#-l 5 branches longer than 5 times meadian branch not in backbone
#-x number of processors availab;e
#-T 0.40 sequences from 0.8 to 1.2 times the median size will be included in the backbone 
#-M 2000 median size of sequences used as backbone
#-s input
#-o output
#-p temporary directory

#recover long sequence names
rm ${gene}_aligned_final.fasta

cat ${gene}_shortnames_alignment.fasta | while read line; do
    if [[ $line == ">"* ]]; then
        grep "$line" ${gene}_fullseqs_final.fasta >> ${gene}_aligned_final.fasta
    else
        echo $line >> ${gene}_aligned_final.fasta
    fi
done

