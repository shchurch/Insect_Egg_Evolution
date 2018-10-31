#!/bin/bash

#SBATCH -J build_rrna                               # Job name 
#SBATCH -o assemble_rrna.%a.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 2                                    # Number of cores/cpus
#SBATCH -t 2-00:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 12000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p general                       # Partition general, serial_requeue, unrestricted, interact


#this script builds rRNA sequences from the Misof et al 2015 dataset
#accessions.csv was made based on their supplement

###load modules
module purge
module load bowtie2/2.2.9-fasrc01
module load trinity/2.2.0-fasrc01
module load legacy
module load centos6/Trimmomatic-0.32

#get file info
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.csv | cut -d , -f 1)
NAME=${NAME// /_} #replace spaces by underlines in name
SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.csv | cut -d , -f 3)


#fetch SRA fastq
echo processing $NAME $SRA
fastq-dump --split-files /n/home08/souzademedeiros/regaldir/ncbi-temp/sra/$SRA.sra

#use bowtie to separate rRNA sequences, ignoring bowtie stdout
bowtie2 --quiet -k 1 -t -p 2 --al-conc ${SRA}_28S.fq -x ./bowtie/28S -1 ${SRA}_1.fastq -2 ${SRA}_2.fastq --un-conc ${SRA}_not_rRNA.fq -S /dev/null
bowtie2 --quiet -k 1 -t -p 2 --al-conc ${SRA}_18S.fq -x ./bowtie/18S -1 ${SRA}_1.fastq -2 ${SRA}_2.fastq --un-conc ${SRA}_not_rRNA.fq -S /dev/null

#use trimmomatic to trim poor quality reads
java -jar $TRIMMOMATIC/trimmomatic-0.32.jar \
PE -threads 2 -phred33 ${SRA}_28S.1.fq ${SRA}_28S.2.fq \
${SRA}_28S_R1.pair.fastq ${SRA}_28S_R1.single.fastq \
${SRA}_28S_R2.pair.fastq ${SRA}_28S_R2.single.fastq \
ILLUMINACLIP:$TRIMMOMATIC/adapters/illuminaClipping_main.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25

java -jar $TRIMMOMATIC/trimmomatic-0.32.jar \
PE -threads 2 -phred33 ${SRA}_18S.1.fq ${SRA}_18S.2.fq \
${SRA}_18S_R1.pair.fastq ${SRA}_18S_R1.single.fastq \
${SRA}_18S_R2.pair.fastq ${SRA}_18S_R2.single.fastq \
ILLUMINACLIP:$TRIMMOMATIC/adapters/illuminaClipping_main.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25


#use Trinity to make assemblies
Trinity --max_memory 11G --CPU 2 --output Trinity_${SRA}_28S_$NAME --left ${SRA}_28S_R1.pair.fastq \
--right ${SRA}_28S_R2.pair.fastq --seqType fq --full_cleanup

Trinity --max_memory 11G --CPU 2 --output Trinity_${SRA}_18S_$NAME --left ${SRA}_18S_R1.pair.fastq \
--right ${SRA}_18S_R2.pair.fastq --seqType fq --full_cleanup

#delete all temporary files
rm ${SRA}*
