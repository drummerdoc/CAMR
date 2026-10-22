#!/usr/bin/env bash
#SBATCH --job-name=camr_eos_gpu
#SBATCH --account=ACCOUNT            # edit
#SBATCH --partition=PARTITION        # the GPU (GH200 / H100) partition; edit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --time=02:00:00              # L0_200 all four backends: GERG on a derated-FP64 card took 28 min; H100 should be far less
#SBATCH --output=camr_eos_gpu_%j.out
# GPU half of the four-backend comparison (Exec/CO2_PipeBreak/bench/OLIVIA.md).
set -euo pipefail
module purge
module load GCC/13 OpenMPI CUDA/12 Python/3     # edit to the site's names
source "$HOME/camr_py/bin/activate"
export AMREX_HOME=${AMREX_HOME:-$HOME/src/amrex}
cd "$HOME/src/CAMR/Exec/CO2_PipeBreak"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
# both tags, GPU only; CORES/GPU_CORE empty = no taskset, SLURM binds
TAG=L0_200    MAX_LEVEL=0 STEPS=200                           PLATFORMS=gpu EOS="PR PRTab GERG GERGTab" CORES="" GPU_CORE="" bench/backend_bench.sh
TAG=L0big_200 MAX_LEVEL=0 STEPS=200 NCELL="1024 512" MAXGRID=128 PLATFORMS=gpu EOS="PR PRTab"             CORES="" GPU_CORE="" bench/backend_bench.sh
