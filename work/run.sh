#!/bin/bash
#SBATCH --partition=serial
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --job-name=o2
#SBATCH --mem=200000

module purge
module load gcc/9.2.0
module load binutils/2.26
module load cmake-3.6.2

export OMP_NUM_THREADS=28;
export MKL_NUM_THREADS=28;
export PYSCF_MAX_MEMORY=$SLURM_MEM_PER_NODE;
echo $PYSCF_MAX_MEMORY

source /home/yangjunjie/intel/oneapi/setvars.sh --force;
export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH

export TMPDIR=/scratch/global/yangjunjie/$SLURM_JOB_NAME/$SLURM_JOB_ID/
export PYSCF_TMPDIR=$TMPDIR
echo TMPDIR       = $TMPDIR
echo PYSCF_TMPDIR = $PYSCF_TMPDIR
mkdir -p $TMPDIR

export PATH=/home/yangjunjie/packages/block2/block2-preview-p0.5.2rc10/pyblock2/driver/:$PATH
export PYTHONPATH=/home/yangjunjie/packages/pyscf/pyscf-main/;
export PYTHONPATH=/home/yangjunjie/packages/block2/block2-preview-p0.5.2rc10/build-py-ext/:$PYTHONPATH;
export PYSCF_EXT_PATH=/home/yangjunjie/intel/oneapi/intelpython/python3.9/lib/python3.9/site-packages/;
alias block2main="/home/yangjunjie/packages/block2/block2-preview-p0.5.2rc10/pyblock2/driver/block2main"

python main.py
