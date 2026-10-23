#!/usr/bin/env bash
#SBATCH --job-name=camr_eos_gpu
#SBATCH --account=nn12019k
#SBATCH --partition=accel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00              # L0_200 all four backends took 32 min on a derated-FP64 card (GERG 28 min); H100 should be far less
#SBATCH --output=camr_eos_gpu_%j.out
# GPU half of the four-backend comparison (Exec/CO2_PipeBreak/bench/OLIVIA.md).
# Runs from the aarch64 (GH200) checkout, which holds the CUDA executables:
#   CAMR2d.gnu.TPROF.MPI.CUDA.PS.<eos>.ex, built on an accel node with CUDA_ARCH=90.
set -euo pipefail
# Lmod here is hierarchical: the NRIS/<stack> module exposes the neoverse_v2 software, so never `module purge`.
# The accel nodes' default init loads the GPU stack; if `module load GCC/13.3.0` fails here, add that stack module
# (see `module list` in an accel allocation) in front of it, as the CPU script does with NRIS/CPU.
module load GCC/13.3.0 OpenMPI/5.0.3-GCC-13.3.0 CUDA/12.6.0
export AMREX_HOME=${AMREX_HOME:-$HOME/src/amrex}
cd "${CAMR_GPU_TREE:-$HOME/src/CAMR-gh200}/Exec/CO2_PipeBreak"
echo "tree $(pwd)  commit $(git rev-parse --short HEAD)  node $(hostname)"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
# both tags, GPU only; CORES/GPU_CORE empty = no taskset, SLURM binds
TAG=L0_200    MAX_LEVEL=0 STEPS=200                              PLATFORMS=gpu EOS="PR PRTab GERG GERGTab" CORES="" GPU_CORE="" bench/backend_bench.sh
TAG=L0big_200 MAX_LEVEL=0 STEPS=200 NCELL="1024 512" MAXGRID=128 PLATFORMS=gpu EOS="PR PRTab"             CORES="" GPU_CORE="" bench/backend_bench.sh
