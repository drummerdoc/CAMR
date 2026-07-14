#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// =====================================================================
//  CO2_PipeBreak  —  Phase 0 of the 2-D CCS pipe-break showcase.
//
//  Planar "guillotine" rupture, realised WITHOUT embedded boundary:
//  the x-lo face is the rupture plane.  Over the break gap
//  |y - gap_yc| <= gap_half it acts as a high-pressure dense-CO2
//  RESERVOIR pressure-inlet (stagnation state p_res, T_res); outside
//  the gap it is a reflecting SLIP WALL (the intact pipe wall / orifice
//  plate).  The domain is initialised to the low-pressure ambient
//  (p0, T0) and the other three faces are non-reflecting outflow.
//
//  The rapid vent drives the dense CO2 across the saturation dome ->
//  flashing under-expanded jet (the demonstration money-shot).
//
//  Defaults: dense/near-critical reservoir venting to a low-P ambient.
//    reservoir : p_res = 90 bar, T_res = 290 K  (dense, just off dome)
//    ambient   : p0    =  5 bar, T0    = 260 K
//    gap       : centred at gap_yc with half-width gap_half [m]
// =====================================================================
struct ProbParmDevice
{
    // Ambient (uniform initial domain state + outflow far-field target).
    amrex::Real p0    = 5.0e5;    // Pa   ambient / initial pressure
    amrex::Real T0    = 260.0;    // K    ambient / initial temperature
    amrex::Real u0    = 0.0;      // m/s  initial velocity
    amrex::Real p_amb = 5.0e5;    // Pa   outflow far-field pressure

    // Reservoir (dense CO2 behind the break gap, imposed at x-lo).
    amrex::Real p_res = 9.0e6;    // Pa   reservoir stagnation pressure
    amrex::Real T_res = 290.0;    // K    reservoir stagnation temperature

    // Two-phase (saturated) reservoir option (rung 4 — saturated CO2 blowdown).
    // res_alpha1 = 1 (default): single-phase reservoir at (p_res, T_res) via
    //   PYT2RE (the original behaviour).
    // res_alpha1 in (0,1): inject a SATURATED two-phase mixture at T_res, with
    //   liquid volume fraction res_alpha1 (reservoir P = Psat(T_res), NOT
    //   p_res).  Exercises MT + relaxation on genuine two-phase cells without
    //   flash nucleation.  The per-phase saturated (rho,e) below are filled
    //   once in amrex_probinit (co2_sat_LV) so bcnormal reads them cheaply.
    amrex::Real res_alpha1 = 1.0; //      liquid volume fraction (1 = single-phase)
    amrex::Real res_rhoL = 0.0;   // kg/m3 saturated-liquid density  @ T_res
    amrex::Real res_eL   = 0.0;   // J/kg  saturated-liquid internal energy
    amrex::Real res_rhoV = 0.0;   // kg/m3 saturated-vapor  density  @ T_res
    amrex::Real res_eV   = 0.0;   // J/kg  saturated-vapor  internal energy

    // Dome-consistent lip taper (task #59).  0 (default): legacy linear blend
    //   s_ext = w*R + (1-w)*W of the (dense) two-phase reservoir R against the
    //   vapor wall-reflection W -- this lands taper-edge cells OFF the
    //   saturation dome, where the stiff PR-EOS rings (dense-liquid pockets,
    //   dt collapse) at showcase resolution.
    // 1: build EVERY in-gap ghost cell as a VALID ON-DOME saturated mixture at
    //   P=Psat(T_res) with liquid fraction tapered as res_alpha1*swin (liquid
    //   confined toward the gap core, -> pure saturated vapor at the edge), and
    //   NO linear blend with W.  Keeps the EOS response continuous across the
    //   lip so round-off stays round-off and no off-dome pockets form.
    int res_dome_taper = 0;
    // Reservoir inflow normal velocity [m/s] (tapered by swin).  0 (default):
    //   stagnation ghost (u=0), flow is pressure-driven -- but the dense liquid
    //   then accumulates at the inflow face until P builds, seeding edge
    //   pockets.  A small positive value carries liquid in so it does not
    //   stagnate/concentrate at the orifice.  Only used with res_dome_taper=1.
    amrex::Real res_u = 0.0;

    // Characteristic (non-reflecting) subsonic two-phase inlet (task #41/#59).
    // 0 (default): the dome-consistent ghost fixes velocity (res_u) AND pressure
    //   (Psat) AND composition -- a full-state Dirichlet inlet, which over-
    //   specifies a subsonic inflow and REFLECTS outgoing acoustic waves
    //   (intermittent inlet pressure flares / core-P noise).
    // 1: keep pressure=Psat and the on-dome composition from the reservoir, but
    //   take the inflow VELOCITY from the outgoing (u-c) acoustic invariant of
    //   the interior:  u_b = u_int + (Psat - P_int)/(rho_int c_int).  This lets
    //   the leaving acoustic characteristic pass through instead of reflecting.
    //   Requires res_dome_taper=1.  res_u is then ignored.
    int res_char_inflow = 0;
    amrex::Real res_Psat = 0.0;  // Pa  Psat(T_res); filled in probinit for the
                                 //     characteristic inlet (P target).

    // Break-gap geometry on the x-lo rupture plane.
    amrex::Real gap_yc   = 0.5;   // m    gap centre (y)
    amrex::Real gap_half = 0.05;  // m    gap half-width
    // Lip taper (task #49): smoothly ramp the reservoir inflow strength from
    // full (in the gap) to ambient over gap_taper OUTSIDE the gap edge, so
    // the inlet is a soft orifice profile rather than a hard top-hat.  A hard
    // reservoir|slip-wall jump is an infinite-shear lip that separates and
    // craters the corner pressure to vacuum (dt collapse).  0 = hard top-hat.
    amrex::Real gap_taper = 0.02; // m

    // Half-domain centerline symmetry (rung 5).  When nonzero, the y-lo face
    // is a slip-wall MIRROR (reflect wall-normal momentum), i.e. a symmetry
    // plane at y = prob_lo_y.  Run y in [0.5, 1.0] with the gap at the y-lo
    // plane to solve only the upper half and exclude the antisymmetric
    // (flapping) flash-lip mode.  0 = off (full domain, default).
    int sym_ylo = 0;

    // Ambient phase branch.  0 (default): single-phase ambient on the phase-1
    // (liquid/supercritical) branch.  1: place the ambient on the phase-2
    // (VAPOR) branch -- required for a subcritical GAS ambient so the flash
    // source does not mistake it for superheated liquid (rung 5).
    int amb_vapor = 0;

    // Reservoir soft-start (task #48): ramp the imposed reservoir pressure
    // linearly from ambient (p0) to p_res over ramp_time, so the inlet is
    // not an impulsive p_res|p0 discontinuity at t=0 — that impulsive
    // startup is what drives the growing inlet pressure oscillation.
    // ramp_time <= 0 gives the instantaneous full reservoir (original).
    amrex::Real ramp_time = 5.0e-4;   // s   (0 disables the ramp)
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
