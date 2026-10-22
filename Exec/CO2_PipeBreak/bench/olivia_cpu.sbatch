#!/usr/bin/env bash
#SBATCH --job-name=camr_eos_cpu
#SBATCH --account=ACCOUNT            # edit
#SBATCH --partition=PARTITION        # the CPU (x86) partition, or the GH200 partition for Grace-core timings; edit
#SBATCH --nodes=1
#SBATCH --ntasks=8                   # override with RANKS=$(nproc) ... for the full-node number
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00              # L0_200: PR 26 s, PRTab 21 s, GERG 342 s, GERGTab 51 s on 8 EPYC ranks; L0big_200 PR/PRTab ~16x
#SBATCH --output=camr_eos_cpu_%j.out
# CPU half of the four-backend comparison (Exec/CO2_PipeBreak/bench/OLIVIA.md).
#   sbatch bench/olivia_cpu.sbatch                                   # 8 ranks
#   RANKS=$(nproc) TAGSUFFIX=_fullnode sbatch --ntasks=$(nproc) bench/olivia_cpu.sbatch
set -euo pipefail
module purge
module load GCC/13 OpenMPI Python/3             # edit to the site's names
source "$HOME/camr_py/bin/activate"
export AMREX_HOME=${AMREX_HOME:-$HOME/src/amrex}
cd "$HOME/src/CAMR/Exec/CO2_PipeBreak"
RANKS=${RANKS:-${SLURM_NTASKS:-8}}; TAGSUFFIX=${TAGSUFFIX:-}
lscpu | grep -E 'Model name|^CPU\(s\)' | head -2
TAG=L0_200$TAGSUFFIX    MAX_LEVEL=0 STEPS=200                           PLATFORMS=cpu EOS="PR PRTab GERG GERGTab" RANKS=$RANKS CORES="" bench/backend_bench.sh
TAG=L0big_200$TAGSUFFIX MAX_LEVEL=0 STEPS=200 NCELL="1024 512" MAXGRID=128 PLATFORMS=cpu EOS="PR PRTab"             RANKS=$RANKS CORES="" bench/backend_bench.sh
