#!/bin/bash

#SBATCH -J download_rrna                               # Job name 
#SBATCH -o download_rrna.%A.out                     # File to which stdout will be written
#SBATCH -N 1                                    # Ensure that all cores are on one machine
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -t 0-01:00                             # Runtime in DD-HH:MM
#SBATCH --mem 500                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p serial_requeue                       # Partition general, serial_requeue, unrestricted, interact


#this script downloads fastq files from the Misof et al 2015 dataset
#accessions.csv was made based on their supplement

###load modules
module purge
#module load sratoolkit/2.7.0-fasrc01
#module load legacy
#module load hpc/aspera

for i in `seq 2 114`
    do
    #get file info
    NAME=$(sed -n "${i}p" accessions.csv | cut -d , -f 1)
    NAME=${NAME// /_} #replace spaces by underlines in name
    SRA=$(sed -n "${i}p" accessions.csv | cut -d , -f 3)

    #fetch SRA - sometimes gets error for no reason, loop until it fetches
    #echo processing $NAME $SRA
    #until prefetch $SRA
    #do 
    #    sleep 10 
    #done

    #prefetch still with problems, trying wget
    if [ -e ../ncbi-temp/sra/$SRA.sra ] 
        then
        echo $SRA downloaded already, skipping
    else
        echo Downloading $SRA \($NAME\) using wget 
        wget --no-verbose -P ../ncbi-temp/sra/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/${SRA:0:3}/${SRA:0:6}/$SRA/$SRA.sra &
    fi
done

wait
