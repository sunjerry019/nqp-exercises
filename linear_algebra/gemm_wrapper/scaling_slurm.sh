#!/bin/bash
#
#SBATCH --job-name=GEMMScaling
#SBATCH --comment="Check GEMMScaling"
#SBATCH --ntasks=10
#SBATCH --partition=cip
#SBATCH --mem-per-cpu=2048
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Yudong.Sun@physik.uni-muenchen.de
#SBATCH --chdir=/home/y/Yudong.Sun/0_modules/NQP/nqp-exercises/linear_algebra/gemm_wrapper
#SBATCH --output=/home/y/Yudong.Sun/0_modules/NQP/nqp-exercises/linear_algebra/gemm_wrapper/slurm/slurm.%j.%N.out
#SBATCH --error=/home/y/Yudong.Sun/0_modules/NQP/nqp-exercises/linear_algebra/gemm_wrapper/slurm/slurm.%j.%N.err.out

# source /etc/profile.d/modules.sh
source /home/y/Yudong.Sun/envvars.sh
source $MODULESHOME/init/bash
source ../../init_modules.sh

mpiexec -n $SLURM_NTASKS python3 scaling.py