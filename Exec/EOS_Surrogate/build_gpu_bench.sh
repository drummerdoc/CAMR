#!/usr/bin/env bash
# Build gpu_bench against a device-callable copy of CAMR's PR headers.
# The HEM_NO_AMREX shim in Source/EOS/PR/hem_saturation_amrex.H defines the
# AMReX qualifier macros as empty and its amrex::max/min/abs helpers as plain
# host inlines; for nvcc both must be __host__ __device__.  The copies live in
# cuda_shim/ (git-ignored) and are regenerated on every build.
set -euo pipefail
cd "$(dirname "$0")"
PR=../../Source/EOS/PR; TAB=../../Source/EOS/PRTab
mkdir -p cuda_shim
sed -E -e 's/^(\s*)#define AMREX_GPU_HOST_DEVICE\s*$/\1#define AMREX_GPU_HOST_DEVICE __host__ __device__/' \
       -e 's/^(\s*)#define AMREX_GPU_DEVICE\s*$/\1#define AMREX_GPU_DEVICE __device__/' \
       -e 's/^(\s*)#define AMREX_FORCE_INLINE inline\s*$/\1#define AMREX_FORCE_INLINE __forceinline__/' \
       -e 's/(namespace Math \{ )inline double abs/\1__host__ __device__ inline double abs/' \
       -e 's/^(\s*template <typename T> )inline T (max|min)/\1__host__ __device__ inline T \2/' \
       $PR/hem_saturation_amrex.H > cuda_shim/hem_saturation_amrex.H
cp $PR/hem_pr_state.H cuda_shim/hem_pr_state.H
ARCH=${ARCH:-sm_120}
nvcc -O3 -std=c++17 -arch=$ARCH -DHEM_NO_AMREX -DAMREX_USE_GPU -Icuda_shim -I$TAB -Xptxas -v \
     gpu_bench.cu $TAB/prtab_table_data.cpp -o gpu_bench "$@"
