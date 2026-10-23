#!/usr/bin/env bash
#SBATCH --job-name=camr_eos_cpu
#SBATCH --account=nn12019k
#SBATCH --partition=small           # the x86 (zen5) CPU partition used for the CPU build; change if the site names it differently
#SBATCH --nodes=1
#SBATCH --ntasks=8                   # 8 ranks to compare with peart; for the full-node number: sbatch --ntasks=<cores> bench/olivia_cpu.sbatch
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00              # peart, 8 EPYC ranks: L0_200 PR 26 s, PRTab 21 s, GERG 342 s, GERGTab 51 s; L0big_200 PR 261 s, PRTab 220 s
#SBATCH --output=camr_eos_cpu_%j.out
# CPU half of the four-backend comparison (Exec/CO2_PipeBreak/bench/OLIVIA.md).
# Runs from the x86 (zen5) checkout, which holds the host executables
#   CAMR2d.gnu.TPROF.MPI.PS.<eos>.ex.  RANKS follows --ntasks.
set -euo pipefail
# Lmod here is hierarchical: NRIS/CPU exposes the zen5 stack, so never `module purge` (it hides GCC/13.3.0).
export LMOD_PAGER=none PAGER=cat GIT_PAGER=cat   # no pager: the module stack's libtinfo breaks the system `more`/`less`
module load NRIS/CPU GCC/13.3.0 OpenMPI/5.0.3-GCC-13.3.0
# zen5 GCCcore/13.3.0 ships a stale include-fixed/pthread.h (glibc mismatch): put the node's own
# pthread.h first.  Set up once with:  mkdir -p ~/fixinc && ln -sf /usr/include/pthread.h ~/fixinc/pthread.h
[ -d "$HOME/fixinc" ] && export CPATH=$HOME/fixinc${CPATH:+:$CPATH}
export AMREX_HOME=${AMREX_HOME:-$HOME/src/amrex}
cd "${CAMR_CPU_TREE:-$HOME/src/CAMR}/Exec/CO2_PipeBreak"
RANKS=${SLURM_NTASKS:-8}; TAGSUFFIX=${TAGSUFFIX:-}
echo "tree $(pwd)  commit $(git rev-parse --short HEAD)  node $(hostname)  ranks $RANKS"
lscpu | grep -E 'Model name|^CPU\(s\)' | head -2
TAG=L0_200$TAGSUFFIX    MAX_LEVEL=0 STEPS=200                              PLATFORMS=cpu EOS="PR PRTab GERG GERGTab" RANKS=$RANKS CORES="" bench/backend_bench.sh
TAG=L0big_200$TAGSUFFIX MAX_LEVEL=0 STEPS=200 NCELL="1024 512" MAXGRID=128 PLATFORMS=cpu EOS="PR PRTab"             RANKS=$RANKS CORES="" bench/backend_bench.sh
