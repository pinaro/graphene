#!/bin/sh
#SBATCH -J pingpong
#SBATCH -n 1
#SBATCH --partition=longq

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

srun python /home/pinar/2018-paper-graphene/notes_pinar/basic_sim/reboot/graphene.py
