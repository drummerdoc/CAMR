#ifndef PS_RETIRED_KEYS_H
#define PS_RETIRED_KEYS_H
// =====================================================================
//  PS_retired_keys.H  --  the one table of retired CAMR.ps_* input keys.
//
//  A deck that sets a key whose behaviour no longer exists must fail
//  loudly rather than run a different physics than it names.  Every
//  retired key is listed here with the reason (or the replacement), and
//  ps_abort_on_retired_keys() is called once at startup from
//  CAMR::read_params.  The EOS backends read their dials from the same
//  CAMR prefix, so their retired keys (eos_*) live here too.  Retired *values* of live keys (CAMR.ps_flux other
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
        { "ps_shear_diss",     "the Jameson shear dissipation is deleted; use ps_mu (physical viscosity)" },
        { "ps_harvest",        "the EOS state harvester is deleted" },
        { "ps_harvest_trace",  "the EOS state harvester is deleted" },
        { "ps_harvest_afloor", "the EOS state harvester is deleted" },
        { "ps_harvest_dlr",    "the EOS state harvester is deleted" },
        { "ps_harvest_de",     "the EOS state harvester is deleted" },
        { "ps_harvest_file",   "the EOS state harvester is deleted" },
        { "ps_harvest_flush",  "the EOS state harvester is deleted" },
        { "eos_warmstart",       "the Newton warm-start is deleted; the bracketed solvers seed themselves" },
        { "eos_warmstart_fixed", "the fixed-iteration seeded solve is deleted" },
        { "eos_mlp",             "use eos_table" },
        { "eos_mlp_auto",        "use eos_table_auto" },
        { "ps_prehydro_diag",  "the [ps_Tdiag] temperature walker is deleted" },
        { "ps_psat_diag",      "the [ps_psat] walker is deleted" },
        { "ps_t2_diag",        "the [ps_T2] stage walker is deleted" },
        { "ps_diag_morph",     "the [ps_morph] census is deleted" },
        { "ps_eovs_diag",      "the energy-overshoot walker is deleted" },
        { "ps_coexit_diag",    "the [PS-COEXIT-TH] report is deleted; ps_coexist_action=1 still aborts on a band exit" },
        { "ps_prdiag",         "the [ps_prdiag] report is deleted" },
        { "ps_relax_diag",     "the [RELAX-PROBE] cell trace is deleted; use ps_cell_diag" },
        { "ps_floor_budget",   "the [PS-FLOORBUDGET] report is deleted; [PS-PROMOTE] still counts floor_skip" },
        { "ps_dm_census",      "the [PS-DM] census is deleted" },
        { "ps_relax_verbose",  "the [ps_relax] anomaly line is deleted; use ps_pres_diag / ps_diag_alpha" },
        { "ps_flash_ev_diag",  "the [PS-FLASH-EV] hook is deleted; [PS-FLASH-REF] prints under ps_pres_diag" },
        { "ps_dilute_closure", "the dilute energy closure is deleted" },
        { "ps_dilute_alpha0",  "the dilute energy closure is deleted" },
        { "ps_dilute_probe",   "the 0-D probe is deleted; the kept self-tests are ps_ptg_selftest, ps_relax_sweep, ps_x3_test" },
        { "ps_m2_test",        "the 0-D probe is deleted; the kept self-tests are ps_ptg_selftest, ps_relax_sweep, ps_x3_test" },
        { "ps_asy1_probe",     "the 0-D probe is deleted; the kept self-tests are ps_ptg_selftest, ps_relax_sweep, ps_x3_test" },
        { "ps_llf_identity",     "the identity-consistent LLF fallback is the only form" },
        { "ps_wp_proj_scale",    "the nondimensional limiter projection is the only form" },
        { "ps_src_p_reproject",  "sources always re-establish P1 = P2 where the state changed" },
        { "ps_mt_target",        "the mass-transfer target is always consistent with the step's carrier and (E.1) path" },
        { "ps_flash_from_absent","the flash always nucleates from a one-absent cell" },
        { "ps_flash_project_sat","the flash always seeds the newborn phase at its saturation density" },
        { "ps_relax_hyst",       "the relax gate's hysteresis band is always on" },
        { "ps_resync_mass",      "the mixture-mass resync is always on" },
        { "ps_bc_nscbc_flash",   "the NSCBC outflow fan is always the flash-aware choked fan" },
        { "ps_p_mode",           "ps_state_from_cons forms the volume-weighted mixture pressure only" },
        { "ps_p_clip",           "ps_state_from_cons forms the volume-weighted mixture pressure only" },
        { "ps_c_mode",           "ps_state_from_cons forms the Wood mixture sound speed only" },
        { "ps_y_floor",          "ps_state_from_cons forms the Wood mixture sound speed only" },
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
