#!/bin/bash
#SBATCH -J mrbayes                               # Job name 
#SBATCH -o mrbayes.%a.%A.out                     # File to which stdout will be written
#SBATCH -n 16                                    # number of ranks
#SBATCH -c 1                                    #number of threads
#SBATCH -t 4-05:00:00                             # Runtime in DD-HH:MM
#SBATCH --mem-per-cpu 3000                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH -p shared                    # Partition general, serial_requeue, unrestricted, interact
#SBATCH --export=level #may be useful for species level analyses in the future
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=souzademedeiros@fas.harvard.edu

# This script runs mrbayes continuously in the cluster
# After starting analysis, mcmcp is set to append
# The analysis is run for 4 days
# After that, the script is called again
# Before resuming, output files are checked to make sure they are compatible
# Sometimes, there is an I/O error in some of the files, and ckp ends up being ahead
# We first make sure the ckp generation is behind all output files, then relaunch mrbayes

#if this is resuming, archive previous output file
if [ "$last_out" ]; then
        echo Zipping output from last run: $last_out
	mv $last_out old_outs/
	cd old_outs
	gzip -9 $last_out
	cd ..
        export last_out=""
fi

module purge
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 mrbayes/3.2.6-fasrc02

export PATH
export LD_LIBRARY_PATH 

infile=$(find . -maxdepth 1 -name '*.nex' | sort | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo $infile

root=$(basename $infile .nex)


if [ -f ${root}.ckp ] #if resuming, make sure the checkpoint is in the correct position
then 
        echo Run resuming, checking if previous output files are OK

	#get number of significant digists
	ndigits=$(tail -n 2 ${root}.run*t | grep -v end | grep -v == | cut -d " " -f 5 | cut -d "." -f 2 | grep -Ev "^$" | wc -L)

	#finds the last tree in tree file with smallest number of trees
	mintrees=$(tail -n 2 ${root}.run*t | grep -v end | grep -v == | cut -d " " -f 5 | cut -d "." -f 2 | grep -Ev "^$" | xargs -I {} printf "%0"$ndigits"d\n" {} | sort | head -n 1)

	#finds the last parameter line in parameter file with smallest number of parameters
	minp=$(tail -n 2 ${root}.run*p ${root}.mcmc | grep -v == | cut -f 1 | grep -Ev "^$" | xargs -I {} printf "%0"$ndigits"d\n" {} | sort | head -n 1)

	#finds last generation in which everything worked and remove leading zeroes
	minall=$(echo -e "$mintrees\n$minp" | sort | head -n 1 | sed 's/^0*//')

	#finds the last multilple of 1,000,000 (i. e. the last checkpoint generation when everything worked)
	lastckp=$(($minall / 1000000 * 1000000))

	#replaces checkpoint generation in the ckp file
	sed -Ei "s/generation: [0-9]+/generation: $lastckp/g" ${root}.ckp 


	srun --time=4-00:00:00 -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK --mpi=pmi2 --export=ALL mb $infile
        

else
	srun --time=4-00:00:00 -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK --mpi=pmi2 --export=ALL mb $infile
fi

#now check if job ran until timeout signal:
#first let's wait a little bit to make sure sacct will be updated
sleep 6000

#exit_code=$(sacct --job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.0 -o ExitCode | tail -n 1 | cut -d ":" -f 2)
#echo Exit code:
#sacct --job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.0 -o ExitCode

#sometimes the srun job id is .1, not .0
#if [ -z "$exit_code" ]
#then
#      echo "\$exit_code is empty, trying again with ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.1"
#      exit_code=$(sacct --job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.1 -o ExitCode | tail -n 1 | cut -d ":" -f 2)
#      sacct --job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.1 -o ExitCode
#fi


#if [ $exit_code == "125" -o $exit_code == "15" ]
echo sacct --job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} -o ExitCode | grep -q ":15"
sacct --job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} -o ExitCode | grep -q ":15"
if [ $? != 0 ] #if different from zero, job did not terminate with a timeout signal
then
	#if did not timeout, do not restart
        echo Mrbayes ended before timeout. Check.
else
	#if concluded with success, call script again
	echo Mrbayes ran until timeout
	echo Turning on append for next run: sed -i "s/\[mcmcp append=yes;\]/mcmcp append=yes;/g" $infile
	sed -i "s/\[mcmcp append=yes;\]/mcmcp append=yes;/g" $infile
	echo Calling next run:
	export last_out=mrbayes.${SLURM_ARRAY_TASK_ID}.${SLURM_ARRAY_JOB_ID}.out 
	sbatch --export last_out --array=$SLURM_ARRAY_TASK_ID -J ${root} ../get_mrbayes_trees_4runs.sh

fi
exit