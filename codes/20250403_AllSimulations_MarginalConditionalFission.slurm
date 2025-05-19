#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J MarginalVsConditionalFission# A MODIFIER.
# walltime (hh:mm::ss)
#SBATCH -t 00:05:00 

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH --tasks-per-node=1 
#SBATCH --array 1-1000


# change working directory
# SBATCH --chdir=.

### In filenames, %j=jobid, %a=index in job array
#SBATCH -o res_bash/test_%j_%a.out # standard out goes to this file
#SBATCH -e res_bash/test_%j_%a.out # standard err goes to this file

# fin des directives PBS
#############################


# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

#############################


module load R/4.2.0 # charger R
Rscript codes/20250403_AllSimulations_MarginalConditionalFission.R # A MODIFIER