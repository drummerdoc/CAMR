#ifndef PS_HLLC_H
#define PS_HLLC_H
// =====================================================================
//  PS_hllc.H  —  HLLC Riemann fan and wave decomposition for the
//                Pelanti-Shyue six-equation model (the face solver of the
//                wave-propagation interior in PS_umeth.cpp).
//
//  Algorithm: three waves per face (S_L, S_M, S_R).  face_from_state
//  builds the per-side state with branch-locked per-phase (P_k, c_k) under
//  the presence dispatch; wave_speeds gives Davis outer speeds and Toro's
//  HLLC contact speed on the mixture; ps_star_state builds each side's
//  star vector with the relaxed-alpha construction (common pressure rise
//  dP behind each acoustic wave, each phase strained by its own bulk
//  modulus; the COTT "HLLC_v" pairing of mixture P with per-phase rho in
//  the phase energies); fluctuations forms the waves and sums them into
//  A∓ΔQ.  See docs/MODEL_AND_ALGORITHM.md §2.
//
//  Contract: input is the NVAR conserved array of each cell; output the
//  waves, speeds and fluctuations.  Every raw wave satisfies
//  W[URHO] = W[UM1RHO1]+W[UM2RHO2] and W[UEDEN] = W[UE1]+W[UE2].  A face
//  that cannot produce a star state is refused with a named cause
//  (PsFlCause) and the caller falls back to LLF; no state is repaired.
//  Every early return is counted (host builds).  No dials are read here.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "PS_guards.H"
#include "PS_presence.H"
#include "PS_promote.H"
#include "PS_constants.H"
#include "PS_util.H"
#include "PS_counters.H"
#include <AMReX_GpuQualifiers.H>
#include <cmath>

namespace PS_HLLC {

// -------------------------------------------------------------------
//  Per-side face state extracted from a CAMR NVAR array.
// -------------------------------------------------------------------
struct Face {
    amrex::Real alpha_1, alpha_2;
    amrex::Real rho_1, rho_2;    // per-phase densities (branch-locked)
    amrex::Real m_1, m_2;        // α_k ρ_k (partial masses; = NVAR[UM*])
    amrex::Real rho_mix;         // m_1 + m_2
    amrex::Real u_n;             // normal velocity
    amrex::Real ux, uy, uz;      // full velocity components
    amrex::Real e_1, e_2;        // per-phase specific internal energy
    amrex::Real E_1, E_2;        // per-phase specific total energy
    amrex::Real P_1, P_2;        // per-phase branch-locked pressure
    amrex::Real P_mix;           // α₁ P₁ + α₂ P₂ (P-S 2014 eq. 5)
    amrex::Real c_frozen;        // Wallis frozen mixture c
    amrex::Real c_1, c_2;        // per-phase c (floored); read only by the
                                 // relaxed-alpha star closure
    bool        valid;
};

// -------------------------------------------------------------------
//  face_from_state: the Face on one side from that cell's conserved
//  array — the live per-face state of the wp interior.  rho_mix is
//  floored at 1e-6 kg/m³; alpha is clamped to [0,1] with no interior
//  floor (exact 0 and 1 mean Absent).  Per-phase (P_k, c_k) come from
//  the branch-locked EOS under the presence dispatch: the host phase
//  (larger alpha) defines the cell, a non-Independent phase takes the
//  host's P and c, an Absent phase also takes the host's rho, e, E as
//  placeholders multiplied by m_k = 0 downstream.  valid = rho_mix and
//  c_frozen finite and positive.  docs §2.1.
// -------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Face face_from_state(int idir, const amrex::Real U[NVAR],
                     const PsPres& pr) noexcept
{
    using amrex::Real;
    Face f;
    f.valid = false;

    f.m_1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    f.m_2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    f.rho_mix = amrex::max(f.m_1 + f.m_2, Real(ps_const::RHO_MIX_MIN));

    Real a1 = ps_finite_or(U[UALPHA1], Real(1.0));
    // Bounds only; out-of-range values are ps_validate_state's to report.
    if (a1 < Real(0.0)) a1 = Real(0.0);
    if (a1 > Real(1.0)) a1 = Real(1.0);
    f.alpha_1 = a1;
    f.alpha_2 = Real(1.0) - a1;

    const PsRegime rg1 = ps_regime(f.alpha_1, f.m_1, pr);
    const PsRegime rg2 = ps_regime(f.alpha_2, f.m_2, pr);
    // Independent: m/alpha is well-conditioned (alpha >= alpha_cond) and
    // is clamped into the EOS domain (G1, PS_guards.H).  Corridor: the
    // ratio is noise but enters the star-state q_k denominator, where
    // its contribution is bounded by the tiny m_k; formed clamped.
    // Absent: never divided; the host's density is substituted below.
    f.rho_1 = (rg1 == PsRegime::Absent) ? Real(0.0)
            : ps_guard::clamp_phase_density(f.m_1 / amrex::max(f.alpha_1, Real(ps_const::DENOM_TINY)),
                                            EOS::rho_min(), EOS::rho_max());
    f.rho_2 = (rg2 == PsRegime::Absent) ? Real(0.0)
            : ps_guard::clamp_phase_density(f.m_2 / amrex::max(f.alpha_2, Real(ps_const::DENOM_TINY)),
                                            EOS::rho_min(), EOS::rho_max());

    // Velocity components.
    const Real inv_r = Real(1.0) / f.rho_mix;
    f.ux = ps_finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    f.uy = ps_finite_or(U[UMY], Real(0.0)) * inv_r;
#else
    f.uy = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    f.uz = ps_finite_or(U[UMZ], Real(0.0)) * inv_r;
#else
    f.uz = Real(0.0);
#endif
    if      (idir == 0) f.u_n = f.ux;
    else if (idir == 1) f.u_n = f.uy;
    else                f.u_n = f.uz;

    // Per-phase specific total and internal energies.
    const Real ke_spec = ps_kinetic_energy(f.ux, f.uy, f.uz);
    const Real E1_tot = ps_finite_or(U[UE1], Real(0.0));
    const Real E2_tot = ps_finite_or(U[UE2], Real(0.0));
    const Real e_mix  = ps_finite_or(U[UEINT], Real(0.0)) * inv_r;
    f.E_1 = (f.m_1 > ps_const::M_TINY) ? E1_tot / f.m_1 : (e_mix + ke_spec);
    f.E_2 = (f.m_2 > ps_const::M_TINY) ? E2_tot / f.m_2 : (e_mix + ke_spec);
    f.e_1 = f.E_1 - ke_spec;
    f.e_2 = f.E_2 - ke_spec;

    // Species carrier (pure CO₂).
    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);

    // Branch-locked per-phase P and c, one solve per phase, host dispatch
    // as in ps_augment_primitives: the host's query defines the cell; a
    // Corridor/Absent phase takes the host's intensives (the corridor
    // closure is the definition of that phase's thermodynamics, not a
    // repair), so its EOS is never queried at an undefined (rho_k, e_k).
    Real c1, c2;
    const bool host_is_1 = (f.alpha_1 >= f.alpha_2);
    if (host_is_1) {
        EOS::REY2PCs_liquid(f.rho_1, f.e_1, Y, f.P_1, c1);
        if (ps_regime_reach(rg2, f.rho_2, f.e_2, /*liquid=*/false, pr)
                == PsRegime::Independent) {
            EOS::REY2PCs_vapor(f.rho_2, f.e_2, Y, f.P_2, c2);
        } else {
            f.P_2 = f.P_1; c2 = c1;
            if (rg2 == PsRegime::Absent) { f.rho_2 = f.rho_1; f.e_2 = f.e_1; f.E_2 = f.E_1; }
        }
    } else {
        EOS::REY2PCs_vapor(f.rho_2, f.e_2, Y, f.P_2, c2);
        if (ps_regime_reach(rg1, f.rho_1, f.e_1, /*liquid=*/true, pr)
                == PsRegime::Independent) {
            EOS::REY2PCs_liquid(f.rho_1, f.e_1, Y, f.P_1, c1);
        } else {
            f.P_1 = f.P_2; c1 = c2;
            if (rg1 == PsRegime::Absent) { f.rho_1 = f.rho_2; f.e_1 = f.e_2; f.E_1 = f.E_2; }
        }
    }
    // G3 (PS_guards.H): per-phase pressure floor against the host's
    // pressure as the reference, consistent with the ctoprim sites.  A
    // slaved phase's pressure is the host's already, so only an
    // Independent minority phase carries a P worth testing.
    {
        const Real P_ref = (f.alpha_1 >= f.alpha_2) ? f.P_1 : f.P_2;
        const Real P_stock = ps_guard::sanitize_stock_pressure(P_ref);
        if (rg1 == PsRegime::Independent) ps_guard::sanitize_phase_pressure(f.P_1, P_stock);
        if (rg2 == PsRegime::Independent) ps_guard::sanitize_phase_pressure(f.P_2, P_stock);
    }
    f.P_mix = f.alpha_1 * f.P_1 + f.alpha_2 * f.P_2;
    // Both c_k are floored to 1 m/s on a non-finite or non-positive EOS
    // return (docs §8.3, guard inventory).
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = ps_const::C_MIN;
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = ps_const::C_MIN;
    f.c_1 = c1;
    f.c_2 = c2;

    // Wallis frozen mixture sound speed c² = Y₁c₁² + Y₂c₂², Y_k = α_kρ_k/ρ_mix
    // (ps_cmix2, the single construction every consumer uses).  Y₂ is
    // formed as 1 − Y₁ so it is bit-exactly zero at α₁ = 1.
    const Real Y1 = f.alpha_1 * f.rho_1 / f.rho_mix;
    const Real Y2 = Real(1.0) - Y1;
    const Real c2_frozen = ps_cmix2(Y1, Y2, c1, c2);
    f.c_frozen = (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                          : amrex::max(c1, c2);

    f.valid = std::isfinite(f.rho_mix) && f.rho_mix > Real(0.0)
           && std::isfinite(f.c_frozen) && f.c_frozen > Real(0.0);
    return f;
}


// -------------------------------------------------------------------
//  Refusal causes.  fluctuations() refuses from several places and the
//  face falls back to LLF; counting per cause separates a spurious
//  refusal (fix: stop refusing) from a legitimate one (fix: keep the
//  fallback consistent), which are different pieces of work.  docs §2.6.
// -------------------------------------------------------------------
enum PsFlCause {
    PS_FL_OK         = 0,
    PS_FL_FACE       = 1,   // face_from_state invalid (rho_mix <= 0 or c <= 0)
    PS_FL_WS_DENOM   = 2,   // mixture mass-flux denominator degenerate
    PS_FL_WS_SM      = 3,   // contact speed non-finite
    PS_FL_WS_PSTAR   = 4,   // unused slot; never set.
                            // Retire-candidate: see docs/DESIGN_DECISIONS.md §7 (O-8).

    PS_FL_ST_DENOM   = 5,   // star: S_K - S_M degenerate
    PS_FL_ST_MASSNEG = 6,   // star: a partial mass came out negative
    PS_FL_ST_RHO     = 7,   // star: mixture density <= 0
    PS_FL_ST_Q       = 8,   // star: a per-phase mass-flux denominator vanished
    PS_FL_ST_ESTAR   = 9,   // star: a phase energy came out non-finite
    PS_FL_NCAUSE     = 10
};
static_assert(PS_FL_NCAUSE == PS_COUNTERS_FL_NCAUSE,
              "PsCounters::fl_cause extent must match PsFlCause");

// -------------------------------------------------------------------
//  wave_speeds: Davis outer speeds and the HLLC contact speed for a face.
//  Each side uses its own frozen c, so the fan is at least as wide as
//  the local sound speed whatever the composition asymmetry; the fan
//  only has to bound the waves (docs §2.2).  S_M from the mixture mass
//  and momentum balance (Toro §10.4) using P_mix.  cause names the
//  refusing test (PS_FL_OK on success) on every return.
// -------------------------------------------------------------------
struct WaveSpeeds {
    amrex::Real S_L, S_M, S_R;
    bool        valid;
    int         cause;
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
WaveSpeeds wave_speeds(const Face& fL, const Face& fR) noexcept
{
    using amrex::Real;
    WaveSpeeds w;
    w.valid = false;
    w.cause = PS_FL_FACE;
    if (!fL.valid || !fR.valid) return w;

    const Real S_L = amrex::min(fL.u_n - fL.c_frozen, fR.u_n - fR.c_frozen);
    const Real S_R = amrex::max(fL.u_n + fL.c_frozen, fR.u_n + fR.c_frozen);

    const Real denom = fL.rho_mix * (S_L - fL.u_n)
                     - fR.rho_mix * (S_R - fR.u_n);
    w.cause = PS_FL_WS_DENOM;
    if (!std::isfinite(denom) || std::abs(denom) < Real(1.0e-30)) return w;

    const Real S_M = (fR.P_mix - fL.P_mix
                    + fL.rho_mix * fL.u_n * (S_L - fL.u_n)
                    - fR.rho_mix * fR.u_n * (S_R - fR.u_n)) / denom;
    w.cause = PS_FL_WS_SM;
    if (!std::isfinite(S_M)) return w;

    //  No contact-pressure positivity test here by decision: the PVRS
    //  estimate fails on strong double rarefactions and its value is read
    //  nowhere; ps_star_state's tests inspect the quantities actually used.
    //  docs §2.2.
    w.S_L    = S_L;
    w.S_M    = S_M;
    w.S_R    = S_R;
    w.valid  = true;
    w.cause  = PS_FL_OK;
    return w;
}


// -------------------------------------------------------------------
//  ps_star_masses: m_k* = r_K m_k for both phases — the per-phase mass
//  Rankine-Hugoniot jump across one discontinuity at S_K in a
//  one-velocity model, m_k*(S_K − S_M) = m_k(S_K − u).  Forced whatever
//  the volume or energy partition does; single-sourced so no star
//  construction can disagree (docs §2.3, C5).
// -------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void ps_star_masses(const Face& fK, amrex::Real r_K,
                    amrex::Real& m1s, amrex::Real& m2s) noexcept
{
    using amrex::Real;
    m1s = fK.m_1 * r_K;
    m2s = fK.m_2 * r_K;
}

//  PsCounters::relax_fb counts the faces where the relaxed-alpha
//  construction fell back to the equal-strain form; reported as
//  [PS-RELAXFB] and reset with the face audit in CAMR_advance.cpp.

// ===================================================================
//  ps_star_state: the HLLC star vector U_star (NVAR layout) on side K
//  behind the wave at S_K, from that side's Face and conserved array.
//  Returns false with *cause set on every refusal (degenerate S_K − S_M,
//  negative star partial mass, rho* <= 0, vanishing q_k, non-finite
//  E_k*); nothing is clamped.  Conserved-slot components follow the
//  standard mixture HLLC algebra; alpha and the phase energies come from
//  the relaxed-alpha construction below.  docs §2.3, §2.6.
// ===================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool ps_star_state(int idir, const Face& fK, const amrex::Real UK[NVAR],
                   amrex::Real S_K, amrex::Real S_M,
                   amrex::Real U_star[NVAR],
                   int* cause = nullptr) noexcept
{
    using amrex::Real;
#define PS_ST_REFUSE(c) do { if (cause) *cause = (c); return false; } while (0)
    const Real denom = S_K - S_M;
    if (!std::isfinite(denom) || std::abs(denom) < Real(1.0e-30))
        PS_ST_REFUSE(PS_FL_ST_DENOM);
    const Real r_K = (S_K - fK.u_n) / denom;
    Real m1s, m2s;
    ps_star_masses(fK, r_K, m1s, m2s);
    if (m1s < Real(0.0) || m2s < Real(0.0)) PS_ST_REFUSE(PS_FL_ST_MASSNEG);
    const Real rhos = m1s + m2s;
    if (rhos <= Real(0.0)) PS_ST_REFUSE(PS_FL_ST_RHO);
    const Real q1 = fK.rho_1 * (S_K - fK.u_n);
    const Real q2 = fK.rho_2 * (S_K - fK.u_n);
    if (std::abs(q1) < Real(1.0e-30) || std::abs(q2) < Real(1.0e-30))
        PS_ST_REFUSE(PS_FL_ST_Q);
    const Real Pe1 = fK.P_mix;   // mixture P in the phase energies (H-3)
    const Real Pe2 = fK.P_mix;

    //  Relaxed-alpha star state (docs §2.3).  Constraints per side, all
    //  inherited from the mixture HLLC algebra: (C1) m1*+m2* = r_K rho,
    //  (C2) u* = S_M, (C3) sum m_k* E_k* = rho* E*_HLLC, (C4) a1*+a2* = 1,
    //  (C5) m_k* = r_K m_k.  Closure: behind the wave both phases sit at a
    //  common pressure rise dP, each strained by its own bulk modulus
    //  B_k = rho_k c_k²:
    //      s_k = dP/B_k,  rho_k* = rho_k(1+s_k),  alpha_k* = alpha_k r_K/(1+s_k),
    //  and (C4) is one monotone scalar equation in dP,
    //      f(dP) = alpha_1 r_K/(1+dP/B1) + alpha_2 r_K/(1+dP/B2) − 1 = 0,
    //  seeded by the Wood-weighted linearisation and solved with three
    //  Newton steps at a fixed count (no tolerance to tune).  Energy: each
    //  phase keeps the kinetic term and is paid P dV for its own volume
    //  change, E_k* = E_k + (S_M−u) S_M + P s_k/(rho_k(1+s_k)); with
    //  f(dP) = 0 the sum m_k* E_k* equals the equal-strain mixture sum
    //  exactly, so (C3) holds and every conserved-slot star component is
    //  unchanged — only the UALPHA1/UE1/UE2 wave components move.  Equal
    //  bulk moduli (every single-phase cell, c-slaved corridor faces)
    //  reduce to the equal-strain form; the equal-strain form is what the
    //  literal Pelanti B.14 gives and it over-compresses the stiff phase
    //  by the bulk-modulus ratio.  Degeneracy (non-finite B_k, 1+s_k <= 0,
    //  vanishing derivative, alpha* outside (0,1)) falls back to the
    //  equal-strain form for that face and is counted, never clamped.
    //  Exact alpha in {0,1} skips the solve: the forms coincide there and
    //  presence needs those values bit-exact.
    //  Equal-strain values first; the relaxed solve overwrites on success.
    Real a1s = fK.alpha_1;
    Real E1s = fK.E_1 + (S_M - fK.u_n) * (S_M + Pe1 / q1);
    Real E2s = fK.E_2 + (S_M - fK.u_n) * (S_M + Pe2 / q2);
    bool relaxed_done = false;
    if (fK.alpha_1 > Real(0.0) && fK.alpha_1 < Real(1.0)) {
        const Real B1 = fK.rho_1 * fK.c_1 * fK.c_1;
        const Real B2 = fK.rho_2 * fK.c_2 * fK.c_2;
        if (std::isfinite(B1) && B1 > Real(1.0e-30) &&
            std::isfinite(B2) && B2 > Real(1.0e-30)) {
            const Real sumc = fK.alpha_1 / B1 + fK.alpha_2 / B2;   // Wood
            Real dP = (sumc > Real(1.0e-30))
                    ? (Real(1.0) - Real(1.0) / r_K) / sumc : Real(0.0);
            bool ok_r = std::isfinite(dP);
            for (int it = 0; ok_r && it < 3; ++it) {
                const Real t1 = Real(1.0) + dP / B1;
                const Real t2 = Real(1.0) + dP / B2;
                if (!(t1 > Real(1.0e-30)) || !(t2 > Real(1.0e-30))) {
                    ok_r = false; break;
                }
                const Real fv = fK.alpha_1 * r_K / t1
                              + fK.alpha_2 * r_K / t2 - Real(1.0);
                const Real fd = -(fK.alpha_1 * r_K / (B1 * t1 * t1)
                                + fK.alpha_2 * r_K / (B2 * t2 * t2));
                if (!std::isfinite(fv) || !(std::abs(fd) > Real(1.0e-30))) {
                    ok_r = false; break;
                }
                dP -= fv / fd;
                if (!std::isfinite(dP)) { ok_r = false; break; }
            }
            if (ok_r) {
                const Real t1 = Real(1.0) + dP / B1;
                const Real t2 = Real(1.0) + dP / B2;
                const Real a1r = fK.alpha_1 * r_K / t1;
                if (t1 > Real(1.0e-30) && t2 > Real(1.0e-30) &&
                    a1r > Real(0.0) && a1r < Real(1.0)) {
                    const Real E1r = fK.E_1 + (S_M - fK.u_n) * S_M
                                   + Pe1 * (dP / B1) / (fK.rho_1 * t1);
                    const Real E2r = fK.E_2 + (S_M - fK.u_n) * S_M
                                   + Pe2 * (dP / B2) / (fK.rho_2 * t2);
                    if (std::isfinite(E1r) && std::isfinite(E2r) &&
                        std::isfinite(a1r)) {
                        a1s = a1r; E1s = E1r; E2s = E2r;
                        relaxed_done = true;
                    }
                }
            }
        }
    }
#if !defined(AMREX_USE_GPU)
    if (!relaxed_done &&
        fK.alpha_1 > Real(0.0) && fK.alpha_1 < Real(1.0)) {
        ++ps_counters().relax_fb;    // equal-strain fallback, counted
    }
#endif
    amrex::ignore_unused(relaxed_done);
    if (!std::isfinite(E1s) || !std::isfinite(E2s)) PS_ST_REFUSE(PS_FL_ST_ESTAR);
#undef PS_ST_REFUSE

    Real ux = fK.ux, uy = fK.uy, uz = fK.uz;   // transverse advected; normal → S_M
    if      (idir == 0) ux = S_M;
    else if (idir == 1) uy = S_M;
    else                uz = S_M;
    const Real ke = ps_kinetic_energy(ux, uy, uz);
    const Real e_mix_s = (m1s * (E1s - ke) + m2s * (E2s - ke)) / rhos;

    for (int n = 0; n < NVAR; ++n) U_star[n] = UK[n];
    U_star[URHO]    = rhos;
    U_star[UMX]     = rhos * ux;
#if (AMREX_SPACEDIM >= 2)
    U_star[UMY]     = rhos * uy;
#endif
#if (AMREX_SPACEDIM == 3)
    U_star[UMZ]     = rhos * uz;
#endif
    U_star[UEINT]   = rhos * e_mix_s;
    U_star[UEDEN]   = rhos * (e_mix_s + ke);   // = m1s·E1s + m2s·E2s
    U_star[UALPHA1] = a1s;   // alpha jumps across the acoustic waves too
    U_star[UM1RHO1] = m1s;
    U_star[UM2RHO2] = m2s;
    U_star[UE1]     = m1s * E1s;
    U_star[UE2]     = m2s * E2s;
    const Real rr = rhos / fK.rho_mix;
    for (int n = 0; n < NUM_SPECIES; ++n)
        U_star[UFS + n] = ps_finite_or(UK[UFS + n], Real(0.0)) * rr;
    return true;
}

// -------------------------------------------------------------------
//  Face audit (host-only, compiled out on GPU, no effect on the
//  numerics; printed as [PS-FACE] / [PS-FLCAUSE] / [PS-RELAXFB] /
//  [PS-W21] under CAMR.ps_face_diag, see CAMR_advance.cpp).  Face
//  classes: 0 = every (side, phase) Independent; 1 = at least one
//  Corridor, no Absent; 2 = at least one Absent (a presence front).
//  Storage is the face section of PsCounters (PS_counters.H):
//    fl_seen / fl_fail / max_incmis, per class: faces seen by
//      fluctuations(), faces refused, and the largest per-fluctuation
//      identity mismatch |Δ(UE1+UE2) − ΔUEDEN|;
//    mis_i / mis_j / mis_k / mis_dir: location of max_incmis, the face's
//      high-side cell (i,j,k) and idir as ps_wp_face sees it; -9999 =
//      nothing recorded;
//    fl_cause: refusals per cause (PsFlCause), not keyed by class;
//    fl_side: per side whose star state refused (0 = L, 1 = R);
//    max_w21_mass / max_w21_energy: worst relative residual of the
//      correction-flux sum identities before the derived assignments in
//      PS_umeth.cpp; must read round-off.
//  Reset by PsCounters::reset_face().
// -------------------------------------------------------------------

//  wp_face_class: audit class of a face from the regimes of both phases
//  on both sides (0 all Independent, 1 Corridor present, 2 Absent present).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int wp_face_class(const amrex::Real UL[NVAR], const amrex::Real UR[NVAR],
                  const PsPres& pr) noexcept
{
    int cls = 0;
    for (int side = 0; side < 2; ++side) {
        const amrex::Real* U = (side == 0) ? UL : UR;
        const amrex::Real a1 = U[UALPHA1];
        const PsRegime r1 = ps_regime(a1, U[UM1RHO1], pr);
        const PsRegime r2 = ps_regime(amrex::Real(1.0) - a1, U[UM2RHO2], pr);
        if (r1 == PsRegime::Absent   || r2 == PsRegime::Absent)   return 2;
        if (r1 == PsRegime::Corridor || r2 == PsRegime::Corridor) cls = 1;
    }
    return cls;
}



struct Fluctuations {
    amrex::Real Am[NVAR];   // A^-ΔQ  (left-going;  updates the L cell)
    amrex::Real Ap[NVAR];   // A^+ΔQ  (right-going; updates the R cell)
    amrex::Real W[3][NVAR]; // waves  (S_L, S_M, S_R)
    amrex::Real s[3];       // wave speeds
    bool        valid;
};

// -------------------------------------------------------------------
//  fluctuations: the wave decomposition of one face from the two cell
//  arrays.  W1 = U_L* − U_L (S_L), W2 = U_R* − U_L* (S_M), W3 = U_R − U_R*
//  (S_R); s·W sums into Am for s < 0 and Ap for s > 0 (a zero-speed wave
//  contributes to neither).  Returns false, with the cause and side
//  counted on host builds, when the face state, the wave speeds or
//  either star state is refused.  For a non-conserved slot Ap ≠ −Am and
//  the non-conservative product is embedded at the face; for a conserved
//  slot Am + Ap = ΔF.  fi/fj/fk locate the face for the audit only.
// -------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool fluctuations(int idir, const amrex::Real UL[NVAR], const amrex::Real UR[NVAR],
                  Fluctuations& f, const PsPres& pr,
                  int fi = -9999, int fj = -9999, int fk = -9999) noexcept
{
    amrex::ignore_unused(fi, fj, fk);
    using amrex::Real;
    f.valid = false;
    for (int n = 0; n < NVAR; ++n) { f.Am[n] = Real(0.0); f.Ap[n] = Real(0.0); }

#if !defined(AMREX_USE_GPU)
    const int fcls_fl = wp_face_class(UL, UR, pr);
    ++ps_counters().fl_seen[fcls_fl];
#endif
#if !defined(AMREX_USE_GPU)
#define PS_W0_FL_FAIL(c) do { ++ps_counters().fl_fail[fcls_fl]; \
                              ++ps_counters().fl_cause[c]; } while (0)
#define PS_W0_FL_SIDE(s_) (++ps_counters().fl_side[s_])
#else
#define PS_W0_FL_FAIL(c)  ((void)0)
#define PS_W0_FL_SIDE(s_) ((void)0)
#endif

    const Face fL = face_from_state(idir, UL, pr);
    const Face fR = face_from_state(idir, UR, pr);
    if (!fL.valid || !fR.valid) { PS_W0_FL_FAIL(PS_FL_FACE); return false; }
    const WaveSpeeds w = wave_speeds(fL, fR);
    if (!w.valid) { PS_W0_FL_FAIL(w.cause); return false; }

    Real UsL[NVAR], UsR[NVAR];
    int stc = PS_FL_OK;
    if (!ps_star_state(idir, fL, UL, w.S_L, w.S_M, UsL, &stc)) {
        PS_W0_FL_FAIL(stc); PS_W0_FL_SIDE(0); return false; }
    if (!ps_star_state(idir, fR, UR, w.S_R, w.S_M, UsR, &stc)) {
        PS_W0_FL_FAIL(stc); PS_W0_FL_SIDE(1); return false; }
#undef PS_W0_FL_FAIL
#undef PS_W0_FL_SIDE

    f.s[0] = w.S_L; f.s[1] = w.S_M; f.s[2] = w.S_R;
    for (int n = 0; n < NVAR; ++n) {
        f.W[0][n] = UsL[n] - UL[n];    // S_L wave
        f.W[1][n] = UsR[n] - UsL[n];   // contact (S_M) wave
        f.W[2][n] = UR[n]  - UsR[n];   // S_R wave
    }
    for (int l = 0; l < 3; ++l) {
        const Real sl = f.s[l];
        if      (sl < Real(0.0)) for (int n = 0; n < NVAR; ++n) f.Am[n] += sl * f.W[l][n];
        else if (sl > Real(0.0)) for (int n = 0; n < NVAR; ++n) f.Ap[n] += sl * f.W[l][n];
    }
    for (int n = 0; n < NVAR; ++n) {
        f.Am[n] = ps_finite_or(f.Am[n], Real(0.0));
        f.Ap[n] = ps_finite_or(f.Ap[n], Real(0.0));
    }
#if !defined(AMREX_USE_GPU)
    {   // Per-fluctuation identity mismatch |Δ(UE1+UE2) − ΔUEDEN|; the star
        // states satisfy the identity by (C3), so this must read round-off
        // on all-Independent faces (docs §2.5).
        const double misM = std::abs(double(f.Am[UE1]) + double(f.Am[UE2])
                                     - double(f.Am[UEDEN]));
        const double misP = std::abs(double(f.Ap[UE1]) + double(f.Ap[UE2])
                                     - double(f.Ap[UEDEN]));
        const double mis = std::max(misM, misP);
        PsCounters& pc = ps_counters();
        if (mis > pc.max_incmis[fcls_fl]) {
            pc.max_incmis[fcls_fl] = mis;
            pc.mis_i[fcls_fl]      = fi;
            pc.mis_j[fcls_fl]      = fj;
            pc.mis_k[fcls_fl]      = fk;
            pc.mis_dir[fcls_fl]    = idir;
        }
    }
#endif
    f.valid = true;
    return true;
}

}  // namespace PS_HLLC

#endif // USE_PS_HYDRO
#endif // PS_HLLC_H
