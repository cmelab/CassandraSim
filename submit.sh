#!/bin/bash -l
#SBATCH -p batch 
#SBATCH -J MosCas
#SBATCH -o job.%j.o
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00

conda activate mc

# PUT YOUR CODE IN HERE!
python -u MosCas.py

