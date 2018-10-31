#!/bin/bash

#SBATCH -J choose_seq                               # Job name 
#SBATCH -o choose_seq.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 2-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 1000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p serial_requeue                       # Partition general, serial_requeue, unrestricted, interact


#this script uses CD-HIT to remove redundant sequences for each sample
#and biopython to choose the longest contig
#it saves all contigs in a fasta file with accession number, and species name
#accessions.csv was made based on Misof et al supplement

###load modules
module purge
module load cd-hit/4.6.4-fasrc02
module load python/2.7.8-fasrc01
source activate PYTHON 


###create folder for all Trinity sequences and cd-hit output
mkdir -p trinity
mkdir -p cdhit
mv Trinity*.fasta trinity/


###run cd-hit on all sequences and run python one-liner to get the longest sequence
for i in `seq 2 114`
    do
    #get file info
    NAME=$(sed -n "${i}p" accessions.csv | cut -d , -f 1)
    UNAME=${NAME// /_} #replace spaces by underlines in name
    SRA=$(sed -n "${i}p" accessions.csv | cut -d , -f 3)

    cd-hit -c 0.98 -M 800 -T 1 -i trinity/Trinity_${SRA}_18S_${UNAME}.Trinity.fasta -o cdhit/unique_18S_${SRA}_${UNAME}.fasta
    python -c "from Bio import SeqIO; dnaseqs = [seq for seq in SeqIO.parse('./cdhit/unique_18S_${SRA}_${UNAME}.fasta','fasta')]; lens = [len(record.seq) for record in dnaseqs]; i = lens.index(max(lens)); dnaseqs[i].description=''; dnaseqs[i].id='accession:' + '$SRA' + '|name:' + '$NAME'; SeqIO.write(dnaseqs[i], open('18S_long.fasta', 'a'), 'fasta') if max(lens) >= 1200 else SeqIO.write(dnaseqs[i], open('18S_short.fasta', 'a'), 'fasta')"
    cd-hit -c 0.98 -M 800 -T 1 -i trinity/Trinity_${SRA}_28S_${UNAME}.Trinity.fasta -o cdhit/unique_28S_${SRA}_${UNAME}.fasta
    python -c "from Bio import SeqIO; dnaseqs = [seq for seq in SeqIO.parse('./cdhit/unique_28S_${SRA}_${UNAME}.fasta','fasta')]; lens = [len(record.seq) for record in dnaseqs]; i = lens.index(max(lens)); dnaseqs[i].description=''; dnaseqs[i].id='accession:' + '$SRA' + '|name:' + '$NAME'; SeqIO.write(dnaseqs[i], open('28S_long.fasta', 'a'), 'fasta') if max(lens) >= 1900 else SeqIO.write(dnaseqs[i], open('28S_short.fasta', 'a'), 'fasta')"
done

