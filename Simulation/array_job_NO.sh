#!/bin/bash
#$ -q short.q
#$ -cwd
#$ -N array_job_NO
## -M mattias.wetscher@uibk.ac.at
## -m e
## comment out next line for nobs %in% c(500, 1000, 5000, 10000)
## -pe openmp 1 
#$ -t 1:100
#$ -l h_vmem=4G
#$ -l h_rt=10:00:00

## PARSING INPUT ARGUMENTS.
## HARDCODE SCALE!!!
if [ $# -ne 3 ] ; then
   printf "Sorry, wrong usage of this script.\n"
   printf "Try again.\n"
   printf "Proper Usage:\n"
   printf "   ${0} <nobs> <nnoise> <rho>\n\n"
   exit 666
else
   nobs=$1
   nnoise=$2
   rho=$3
fi

## Determine the runtime based on the value of 'what'
##if [ "$nobs" -le 5000 ]; then
##   queue="short.q"
##   h_rt="5:00:00"
##   h_vmem="8G"
##elif [ "$nobs" -eq 10000 ]; then
##   queue="short.q"
##   h_rt="10:00:00"
##   h_vmem="12G"
##else
##   queue="short.q"
##   h_rt="10:00:00"
##   h_vmem="20G"
##fi


## Just for testing
printf "Got: nrep = %d, nobs = %d, nnoise = %d, rho = %d\n\n" "${SGE_TASK_ID}" "${nobs}" "${nnoise}" "${rho}"

##set working directory
cd /scratch/c4031051/stagewise/simulation/

## setup R within Anaconda
module load Anaconda3/2023.03/r-4.2-conda-2023.03
source $UIBK_CONDA_PROFILE
conda activate r-base
 
## Starting the R script
/usr/site/hpc/bin/sysconfcpus -n $NSLOTS Rscript --no-save array_job_NO.R --nobs ${nobs} --nnoise ${nnoise} --rho ${rho}

##close anaconda
conda deactivate
