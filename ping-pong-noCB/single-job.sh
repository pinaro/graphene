#!/bin/sh
#SBATCH -J noCB-ping
#SBATCH -n 1
#SBATCH --partition=longq

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

srun python /home/pinar/2018-paper-graphene/ping-pong-noCB/graphene.py
