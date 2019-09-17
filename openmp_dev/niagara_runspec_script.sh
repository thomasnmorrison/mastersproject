#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=1:00:00
#SBATCH --job-name=runspec_job
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

module load intel/2018.2
module load fftw/3.3.7

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./runspec
