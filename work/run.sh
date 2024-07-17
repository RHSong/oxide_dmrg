#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --job-name=o2
#SBATCH --mem=500000

source /home/yangjunjie/.bashrc
module purge
module load gcc/9.2.0
module load binutils/2.26
module load cmake-3.6.2

export OMP_NUM_THREADS=28;
export MKL_NUM_THREADS=28;
export PYSCF_MAX_MEMORY=$SLURM_MEM_PER_NODE;
echo $PYSCF_MAX_MEMORY

conda activate py310-block2
export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH

export TMPDIR=/scratch/global/yangjunjie/$SLURM_JOB_NAME/$SLURM_JOB_ID/
export PYSCF_TMPDIR=$TMPDIR
echo TMPDIR       = $TMPDIR
echo PYSCF_TMPDIR = $PYSCF_TMPDIR
mkdir -p $TMPDIR

python main.py

