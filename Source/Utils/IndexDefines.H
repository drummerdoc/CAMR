#ifndef INDEX_DEFINES_H
#define INDEX_DEFINES_H

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>
#include <AMReX_GpuMemory.H>

#define Interior   0
#define Inflow     1
#define Outflow    2
#define Symmetry   3
#define SlipWall   4
#define NoSlipWall 5

#define NUM_SPECIES 3

#define URHO 0
#define UMX 1
#define UMY 2

#define QRHO 0
#define QU 1
#define QV 2

#define GDU 1
#define GDV 2

#if (AMREX_SPACEDIM == 2)
#define UEDEN 3
#define UEINT 4
#define UTEMP 5
#define UFA 6

#define QGAME 3
#define QPRES 4
#define QREINT 5
#define QTEMP 6
#define QFA 7

#define GDPRES 3
#define GDGAME 4
#define NGDNV 5

#define NTHERM 6
#define QTHERM 7

#elif (AMREX_SPACEDIM == 3)
#define UMZ 3
#define UEDEN 4
#define UEINT 5
#define UTEMP 6
#define UFA 7

#define QW 3
#define QGAME 4
#define QPRES 5
#define QREINT 6
#define QTEMP 7
#define QFA 8

#define GDW 3
#define GDPRES 4
#define GDGAME 5
#define NGDNV 6

#define NTHERM 7
#define QTHERM 8
#endif

#define NUM_ADV 0

#define QGAMC 0
#define QC 1
#define QCSML 2
#define QDPDR 3
#define QDPDE 4
#define NQAUX 5

#ifndef CAMR_USE_AUX
#define NUM_AUX 0
#endif

// ---------------------------------------------------------------------
//  Pelanti-Shyue six-equation state extension (co2-eos-cfd Phase 4).
//
//  When USE_PS_HYDRO is defined (per-case, via -DUSE_PS_HYDRO in the
//  Exec GNUmakefile), the conservative array is extended by five
//  components — volume fraction α₁, per-phase densities α_k ρ_k, and
//  per-phase total energy densities α_k ρ_k E_k.  These extra slots
//  are declared as passive scalars so that Godunov / MOL harmlessly
//  advect them along the contact velocity; the PelantiShyue hydro
//  module (Phase 4c) will replace that treatment with the correct
//  A±ΔQ wave-propagation update.
//
//  When USE_PS_HYDRO is NOT defined (default) the layout is identical
//  to upstream CAMR — no additional slots, no changes to NVAR or the
//  passive-scalar map.
// ---------------------------------------------------------------------
#ifdef USE_PS_HYDRO
#define NUM_PS 5
#else
#define NUM_PS 0
#endif

#define UFS (UFA + NUM_ADV)
#define UFX (UFS + NUM_SPECIES)
#define UPS (UFX + NUM_AUX)
#define NVAR (NTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_PS)

#define QFS (QFA + NUM_ADV)
#define QFX (QFS + NUM_SPECIES)
#define QPS (QFX + NUM_AUX)
#define QVAR (QTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_PS)
#define NQ QVAR

#define NPASSIVE (NUM_SPECIES + NUM_ADV + NUM_AUX + NUM_PS)

#ifdef USE_PS_HYDRO
// PS conservative slots
#define UALPHA1 (UPS + 0)   // volume fraction α₁                (dimensionless)
#define UM1RHO1 (UPS + 1)   // α₁ ρ₁                             [kg/m³]
#define UM2RHO2 (UPS + 2)   // α₂ ρ₂                             [kg/m³]
#define UE1     (UPS + 3)   // α₁ ρ₁ E₁ (per-phase total energy) [J/m³]
#define UE2     (UPS + 4)   // α₂ ρ₂ E₂                          [J/m³]

// Berger-LeVeque fluctuation register width (task #22): the phase-energy
// WP-vs-Godunov defect, one flux-form component per phase.  Its two
// components map to state slots UE1, UE2 (contiguous, so a single
// destcomp=UE1 offset applies the whole register).
// Per-face scratch (fcorr) width for the defect fluctuation register:
//   [0] = UE1 defect,  [1] = UE2 defect   (one-sided, flux-form; #22 P3).
// (The non-conservative alpha C-F correction is handled by the
// capacity-form reflux in CAMR::reflux() — task #36 — which moves alpha1
// with its own conservatively-refluxed mass and needs no fcorr slots.)
#define NUM_PS_FLUCT 2

// PS primitive slots (populated by PS_ctoprim in Phase 4c)
#define QALPHA1 (QPS + 0)   // α₁
#define QRHO1   (QPS + 1)   // ρ₁ per-phase density
#define QRHO2   (QPS + 2)   // ρ₂
#define QP1     (QPS + 3)   // P₁ per-phase pressure
#define QP2     (QPS + 4)   // P₂
#endif

struct PassMap
{
  amrex::GpuArray<int, NPASSIVE> upassMap = {{0}};
  amrex::GpuArray<int, NPASSIVE> qpassMap = {{0}};
};

void init_pass_map(PassMap* pmap);

#endif
