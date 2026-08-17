#ifndef PS_PRESENCE_H
#define PS_PRESENCE_H
// =====================================================================
//  PS_presence.H — the presence-discrete phase-state parameters and
//  regime classifier (S1 of DESIGN_ps_presence_discrete.md).
//
//  THE three constants of the design, single-sourced (§2 of the note):
//
//    alpha_vanish (1e-8)  death threshold — below this a phase is folded
//                         to exactly zero (standalone-validated value;
//                         decade sweep {1e-9, 1e-7} measured bit-identical
//                         on the 20-case battery, WORKLOG 2026-08-15 M-B).
//    alpha_cond   (2e-2)  conditioning bound — at or above this
//                         rho_k = m_k/alpha_k is trustworthy.  RE-DERIVED
//                         2026-08-15 (Stage 5, after the design's own
//                         insensitivity sweep finally ran and PASSED):
//                         eta_max/eps with the 2-D PRODUCTION eta_max =
//                         9.6e-4 (design §14) and eps = 5 % gives 1.92e-2
//                         -> 2e-2.  The 1-D-only value was 1e-2 (eta_max =
//                         5.2e-4); the sweep over [4.2e-3, 2e-2] moved no
//                         case above scheme error (largest 0.014 abs on
//                         B7's 0.66 P error), so the deployment-wide bound
//                         is adopted.
//    alpha_birth  (4e-2)  flash-nucleation seed level (= 2*alpha_cond,
//                         factor unchanged; the factor sweep {1.5, 2, 4}
//                         was measured BIT-IDENTICAL on the battery).
//                         Note this departs from the flash kernel's
//                         historical alpha_seed_target = 0.02; the seed is
//                         single-sourced from here (PS_sources.H S3), and
//                         the 4e-2 point was measured directly.
//
//  Regimes per phase (design §1):
//    ABSENT      alpha_k <= alpha_vanish, or m_k <= 0: the phase does
//                not exist; no intensive quantity is ever formed.
//    CORRIDOR    alpha_vanish < alpha_k < alpha_cond: transport-only;
//                conserved quantities advect, intensives are the HOST
//                phase's (definition, not repair).
//    INDEPENDENT alpha_k >= alpha_cond: full six-equation phase.
//
//  GPU rules (design §12.2): this is a by-value POD; the getter does its
//  ParmParse read ONCE on the host; kernels receive the struct by value
//  (the PS_umeth capture pattern).  No getenv, no statics in kernels.
//  The classifier reads alpha_k and phase totals ONLY — never species
//  (design §13, multicomponent discipline).
//
//  CAMR.ps_presence = 0 (default): every consumer follows its legacy
//  path bit-for-bit.  = 1: presence-aware paths active.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuQualifiers.H>

struct PsPres {
    amrex::Real alpha_vanish = amrex::Real(1.0e-8);
    amrex::Real alpha_cond   = amrex::Real(2.0e-2);   // re-derived 2026-08-15
    amrex::Real alpha_birth  = amrex::Real(4.0e-2);   // = 2*alpha_cond
    //  The legacy (enabled == 0) path is DELETED.  Presence semantics are
    //  the only semantics: alpha free in [0,1] so ABSENT is reachable, phase
    //  states formed by checked construction, branch-locked queries only for
    //  phases that exist.  Every manufacture found in the 2026-08-10 audit
    //  lived in an `enabled == 0` fork whose presence-side twin was already
    //  correct.  This field stays 1 and is retained only until the remaining
    //  forks are removed; `CAMR.ps_presence` is no longer read.
    int         enabled      = 1;
    //  CAMR.ps_cmix_model -- WHICH MIXTURE SOUND SPEED a two-phase cell gets.
    //  A NUMERICS dial riding the presence params struct because that struct is
    //  already the by-value carrier reaching all three live c_mix sites
    //  (PS_hllc.H face_from_state, PS_wavespeed.H, PS_umeth.cpp) with no
    //  signature churn and no ParmParse inside a device-callable function
    //  (GPU rules, DESIGN_ps_presence_discrete.md 12.2 rule 1).  The naming is
    //  a compromise and is flagged here rather than hidden.
    //
    //  WHY IT EXISTS.  Blending two sound speeds into one asserts the phases
    //  are finely intermixed, so that a wave crossing the cell samples many
    //  alternating layers and feels an average.  That is a MORPHOLOGY
    //  assumption, the same one measured wrong for theta by B11 -- but applied
    //  to EVERY two-phase cell of EVERY case, and feeding S_L/S_R, hence every
    //  flux and the timestep.  See LITERATURE_relaxation_rates.md 9.
    //
    //    0  FROZEN mass-weighted (default, bit-identical):
    //         c^2 = Y_1 c_1^2 + Y_2 c_2^2
    //    1  MAX: c = max(c_1, c_2).  The contact-appropriate choice -- at a
    //         resolved interface a wave travels in ONE phase, so the fastest
    //         real signal is the faster phase speed.  Also the safest bound
    //         for S_L/S_R, hence more diffusive, never narrower than truth.
    //    2  WOOD / equilibrium:  1/(rho c^2) = alpha_1/(rho_1 c_1^2)
    //                                        + alpha_2/(rho_2 c_2^2)
    //         The strongly-dispersed choice, and famously LOWER than either
    //         phase at intermediate alpha.  Physically right for a fine
    //         mixture in equilibrium; as an S_L/S_R estimate it is narrower
    //         than the frozen speed, so if the fan stops bounding the true
    //         waves this is where it shows.
    int         cmix_model   = 0;
    //  CAMR.ps_rk_model -- HOW THE STAR-STATE PARTIAL MASSES CONTRACT through
    //  the acoustic fan (D7 instrument; PLAN_measurements_and_fixes.md S1 C5,
    //  REVIEW_RESPONSE.md D7).  Rides this struct for the same reason
    //  cmix_model does: it is already the by-value carrier reaching the star
    //  builders (PS_hllc.H hllc_flux / ps_star_state) with no ParmParse in a
    //  device-callable function.
    //
    //    0  SHARED (default, bit-identical): both phase masses scale by the
    //         one ratio r_K = (S_K - u)/(S_K - S_M) (Pelanti 2022 B.14).
    //         Encodes "acoustic compression does not change composition" --
    //         true for a dispersed mixture, questionable at a contact where
    //         the acoustic wave lives in one phase only.
    //    1  PER-PHASE (instrument, NOT a design): each phase contracts by its
    //         own acoustic ratio r_k built from S_{K,k} = u -/+ c_k, then both
    //         are renormalised by a common factor so that
    //         m_1 r_1' + m_2 r_2' = (m_1+m_2) r_K exactly -- the MIXTURE mass
    //         keeps the HLLC Rankine-Hugoniot contraction, conservation is
    //         untouched, and only the PARTITION of the compression between
    //         the phases moves.  The per-phase energy fluxes (q_k path) are
    //         deliberately unchanged so the dial isolates the mass-partition
    //         question.  Any degeneracy falls back to SHARED.
    int         rk_model     = 0;
};

//  Single-sourced mixture sound speed SQUARED.  One definition for all three
//  live call sites -- the four-way duplication of exactly this expression is
//  what left three copies carrying the wrong task-#199 alpha factor
//  (STATUS_multiphase.md 6.1), so the replacement is not duplicated again.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ps_cmix2(int model, amrex::Real Y1, amrex::Real Y2,
                     amrex::Real c1, amrex::Real c2,
                     amrex::Real alpha_1, amrex::Real alpha_2,
                     amrex::Real rho_1, amrex::Real rho_2, amrex::Real rho_mix)
{
    using amrex::Real;
    if (model == 1) {                    // MAX -- contact-appropriate
        const Real cm = amrex::max(c1, c2);
        return cm * cm;
    }
    if (model == 2) {                    // WOOD / equilibrium
        //  Guarded: a vanishing alpha_k contributes nothing (its phase owns no
        //  volume), which is the correct limit and avoids 0/0.  If neither
        //  phase contributes, fall through to the frozen form rather than
        //  invent a number.
        Real inv = Real(0.0);
        if (alpha_1 > Real(0.0) && rho_1 > Real(0.0) && c1 > Real(0.0)) {
            inv += alpha_1 / (rho_1 * c1 * c1);
        }
        if (alpha_2 > Real(0.0) && rho_2 > Real(0.0) && c2 > Real(0.0)) {
            inv += alpha_2 / (rho_2 * c2 * c2);
        }
        if (inv > Real(0.0) && rho_mix > Real(0.0)) {
            return Real(1.0) / (rho_mix * inv);
        }
    }
    return Y1 * c1 * c1 + Y2 * c2 * c2;  // 0 = FROZEN mass-weighted (default)
}

// Per-phase regime.  alpha is THIS phase's volume fraction, m its
// partial mass.  Degenerate inputs (NaN alpha) classify as ABSENT so a
// bad cell contributes nothing rather than garbage; ps_validate_state
// is the instrument that reports it.
enum class PsRegime : int { Absent = 0, Corridor = 1, Independent = 2 };

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PsRegime ps_regime(amrex::Real alpha, amrex::Real m, const PsPres& pr) noexcept
{
    if (!(alpha > pr.alpha_vanish) || !(m > amrex::Real(0.0)))
        return PsRegime::Absent;
    if (alpha < pr.alpha_cond) return PsRegime::Corridor;
    return PsRegime::Independent;
}

// =====================================================================
//  Contract 3: a phase state is a CHECKED construction.
//
//  rho_k = m_k / alpha_k  and  e_k = E_k / m_k  have DIFFERENT denominators
//  -- alpha_k partitions volume, m_k partitions mass -- so neither
//  constrains the other and a test on one does NOT license the other.  A
//  phase therefore either HAS a state (both partitions meaningfully
//  non-zero, both quotients well posed) or it has NONE, and consumers of
//  intensive properties must take the no-state branch rather than receive a
//  clamped number.
//
//  This is the predicate the code already owns -- ps_regime's Absent edge --
//  applied at the point where the quotients are FORMED rather than repaired
//  afterwards.  See FINDINGS Addenda 10k-10l.
// =====================================================================
struct PsPhaseQuot {
    amrex::Real rho    = amrex::Real(0.0);
    amrex::Real e      = amrex::Real(0.0);
    bool        exists = false;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PsPhaseQuot ps_phase_quot(amrex::Real alpha, amrex::Real m, amrex::Real E,
                          amrex::Real ke, const PsPres& pr) noexcept
{
    PsPhaseQuot q{};
    if (ps_regime(alpha, m, pr) == PsRegime::Absent) return q;   // no state exists
    q.rho = m / alpha;
    q.e   = E / m - ke;
    if (!(q.rho > amrex::Real(0.0)) || !(q.rho == q.rho) || !(q.e == q.e)) {
        q.rho = amrex::Real(0.0); q.e = amrex::Real(0.0);
        return q;                                                // still exists=false
    }
    q.exists = true;
    return q;
}

// Host-side one-shot read (ParmParse), captured by value into kernels.
inline const PsPres& ps_presence_params()
{
    static const PsPres p = []() {
        PsPres v;
        amrex::ParmParse pp("CAMR");
        //  CAMR.ps_presence deliberately NOT read: the legacy path is gone.
        pp.query("ps_alpha_cond",       v.alpha_cond);
        pp.query("ps_alpha_birth",      v.alpha_birth);
        pp.query("ps_cmix_model",       v.cmix_model);
        pp.query("ps_rk_model",         v.rk_model);
        pp.query("ps_presence_vanish",  v.alpha_vanish);   // unified with the
        // fold's CAMR.ps_alpha_vanish at S3 (that key currently means
        // "fold threshold, 0 = fold off"; this one is the presence
        // classifier's death edge and must stay > 0).
        return v;
    }();
    return p;
}

#endif  // USE_PS_HYDRO
#endif  // PS_PRESENCE_H
