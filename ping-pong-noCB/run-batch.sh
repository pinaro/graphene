#!/bin/bash
# loop over number of jobs
# print the number (so we have some visual progress indicator)
# then submit the  jobs to SLURM
#
FILE="/home/pinar/2018-paper-graphene/ping-pong-noCB/single-job.sh"

for COUNT in {1..50}; do
 echo ${COUNT}
 sbatch -o ${FILE}.${COUNT}.stdout.txt -e ${FILE}.${COUNT}.stderr.txt ${FILE}
 sleep 1 # pause to be kind to the scheduler
done
