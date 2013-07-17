#!/bin/bash
#PBS -N himeno 
#PBS -l mppwidth=64
#PBS -l mppnppn=4
#PBS -l mppdepth=4
#PBS -l walltime=1:00:00
#PBS -j oe
cd $PBS_O_WORKDIR
date
export OMP_NUM_THREADS=4
export CRAY_CUDA_PROXY=1
aprun -n 64  -N 4 -d 4 ./a.out+pat

