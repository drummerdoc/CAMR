#ifndef PS_RETIRED_KEYS_H
#define PS_RETIRED_KEYS_H
// =====================================================================
//  PS_retired_keys.H  --  the one table of retired CAMR.ps_* input keys.
//
//  A deck that sets a key whose behaviour no longer exists must fail
//  loudly rather than run a different physics than it names.  Every
//  retired key is listed here with the reason (or the replacement), and
//  ps_abort_on_retired_keys() is called once at startup from
//  CAMR::read_params.  Retired *values* of live keys (CAMR.ps_flux other
//  than wp, CAMR.ps_relax_mode=3) are refused at their accessors.
//  Docs/DESIGN_DECISIONS.md §6.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_ParmParse.H>
#include <AMReX.H>
#include <string>

inline void ps_abort_on_retired_keys()
{
    struct RetiredKey { const char* key; const char* reason; };
    static const RetiredKey table[] = {
        { "ps_recon",          "reconstruction is not selectable; wp is the only interior scheme" },
        { "ps_alpha_limiter",  "reconstruction is not selectable; wp is the only interior scheme" },
        { "ps_pk_energy_flux", "the mixture-pressure phase-energy flux is the only form" },
        { "ps_ctu",            "no CTU path" },
        { "ps_star_relaxed",   "the relaxed-alpha star state is the only star state" },
        { "ps_rk_model",       "the relaxed-alpha star state is the only star state" },
        { "ps_alpha_vanish",   "use ps_presence_vanish" },
        { "ps_presence",       "presence is always on" },
        { "ps_cmix_model",     "corridor sound speed is host-slaved (no selector)" },
        { "ps_theta_model",    "use ps_theta_tau" },
        { "ps_hrm_theta",      "the HRM correlation is not available" },
        { "ps_bc_nscbc_v2",    "the characteristic outflow has one formulation; use ps_bc_use_nscbc / ps_bc_nscbc_*" },
    };
    amrex::ParmParse pp("CAMR");
    for (const RetiredKey& rk : table) {
        if (pp.contains(rk.key)) {
            amrex::Abort(std::string("CAMR.") + rk.key + " is retired: " + rk.reason);
        }
    }
}

#endif  // USE_PS_HYDRO
#endif  // PS_RETIRED_KEYS_H
