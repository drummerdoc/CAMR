#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// =====================================================================
//  Exec/CO2_RiemannSuite  —  generalized CO2 Riemann problem
//
//  A single executable that can run ANY case from the standalone
//  co2-eos-cfd shock-tube suite (families A/B/C) via prob.* params.
//  Extends CO2_B4 with a per-side PHASE selector and a quality field
//  so that saturated / two-phase / cross-critical left/right states
//  can be initialised directly (mirroring the standalone helpers
//  sTP(T,P,phase), sTx(T,x), sL(T)=satL, sV(T)=satV).
//
//  Per-side spec (L = left of diaphragm, R = right):
//    phase_{L,R}  — how to build that side's thermodynamic state:
//        0 = VAPOR      : single-phase vapor at (T, P)      [sTP Vapor]
//        1 = LIQUID     : single-phase liquid at (T, P)     [sTP Liquid]
//        2 = SC / AUTO  : supercritical / auto-root at (T,P)[sTP Vapor hint]
//        3 = TWOPHASE   : two-phase mixture at (T, quality x_qual)
//                          [sTx].  x_qual=0 => sat liquid (sL),
//                          x_qual=1 => sat vapor (sV).
//    T_{L,R}      — temperature [K]
//    p_{L,R}      — pressure [Pa]   (ignored for TWOPHASE; Psat(T) used)
//    x_qual_{L,R} — vapor quality [-] (only used for phase==TWOPHASE)
//    u_{L,R}      — velocity [m/s]
//
//  The 6-equation P-S state is built from these using CAMR's EOS /
//  hem helpers (state_from_T_P_metastable, state_from_T_x, satState*).
//  See prob.H CAMR_initdata for the exact mapping to (alpha_1, m_k,
//  E_k).  All keys are ParmParse-overridable — see prob.cpp.
// =====================================================================

// Phase selector codes (kept as int for easy ParmParse from inputs).
#define PS_PHASE_VAPOR    0
#define PS_PHASE_LIQUID   1
#define PS_PHASE_SC       2
#define PS_PHASE_TWOPHASE 3

struct ProbParmDevice
{
    // ---- Left state -------------------------------------------------
    amrex::Real p_L      = 1.0e7;    // Pa
    amrex::Real T_L      = 270.0;    // K
    amrex::Real u_L      = 0.0;      // m/s
    amrex::Real x_qual_L = 1.0;      // vapor quality (TWOPHASE only)
    int         phase_L  = PS_PHASE_LIQUID;

    // ---- Right state ------------------------------------------------
    amrex::Real p_R      = 5.0e6;    // Pa
    amrex::Real T_R      = 350.0;    // K
    amrex::Real u_R      = 0.0;      // m/s
    amrex::Real x_qual_R = 1.0;      // vapor quality (TWOPHASE only)
    int         phase_R  = PS_PHASE_VAPOR;

    // Trace volume fraction for the (single-phase) minority phase in
    // the 6-eq representation.  1e-6 matches the standalone winning
    // B4 config.
    amrex::Real alpha_trace = 1.0e-6;

    // Diaphragm position (fraction of x-extent).
    amrex::Real x_diaph  = 0.5;

    // Two-face non-reflecting BC ambient targets (see prob.H bcnormal).
    amrex::Real p_amb_lo = 1.0e7;
    amrex::Real p_amb_hi = 5.0e6;
    amrex::Real p_amb    = 5.0e6;    // NSCBC single-ambient fallback
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
