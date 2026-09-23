#ifndef PS_DIALS_H
#define PS_DIALS_H
// =====================================================================
//  PS_dials.H — the one way a Pelanti-Shyue dial is read.  Each helper
//  reads CAMR.<key> with its default and adds the resolved value back to
//  the ParmParse table, so job_info records the configuration that ran
//  even when the deck was silent.  Host-only (ParmParse, std::string).
//
//  Contract: a dial has exactly one read site, and that site calls the
//  helper from a function-local static initialiser (or an equivalent
//  once-per-process cache), so the table is queried and extended once
//  per key.  Accessors that validate, resolve an auto value or couple
//  several keys keep their logic and use the helper for the raw read.
//  No environment fallbacks (rule 20).
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <string>

//  CAMR.<key> as an int, default `def`; the resolved value is added back.
inline int ps_dial_int(const char* key, int def)
{
    int v = def;
    amrex::ParmParse pp("CAMR");
    pp.query(key, v);
    pp.add(key, v);
    return v;
}

//  CAMR.<key> as a Real, default `def`; the resolved value is added back.
inline amrex::Real ps_dial_real(const char* key, amrex::Real def)
{
    amrex::Real v = def;
    amrex::ParmParse pp("CAMR");
    pp.query(key, v);
    pp.add(key, v);
    return v;
}

//  CAMR.<key> as a string, default `def`; the resolved value is added back.
inline std::string ps_dial_string(const char* key, const char* def)
{
    std::string v(def);
    amrex::ParmParse pp("CAMR");
    pp.query(key, v);
    pp.add(key, v);
    return v;
}

//  CAMR.ps_sigma (0): the transported interfacial-area density Sigma
//  (Sigma_Type, docs/PLAN_sigma.md).  0 = allocated, zero, never
//  transported, not plotted; 1 = transported on the HLLC contact speed,
//  plotted and checkpointed, read by nothing; 2 = as 1 and read by the
//  sources (commit 4).  Any other value aborts.  One read per process.
inline int ps_sigma_mode()
{
    static const int cached = []() -> int {
        const int v = ps_dial_int("ps_sigma", 0);
        if (v < 0 || v > 2) {
            amrex::Abort("CAMR.ps_sigma must be 0, 1 or 2");
        }
        return v;
    }();
    return cached;
}

//  CAMR.ps_sigma_restart_init (abort): how Sigma is filled when a run at
//  ps_sigma >= 1 restarts from a checkpoint that does not hold it.
//  abort = stop with a message naming this key; zero = Sigma = 0; alg =
//  the algebraic closure (commit 3).  Read only on that path.
inline const std::string& ps_sigma_restart_init()
{
    static const std::string cached = []() -> std::string {
        std::string v = ps_dial_string("ps_sigma_restart_init", "abort");
        if (v != "abort" && v != "zero" && v != "alg") {
            amrex::Abort("CAMR.ps_sigma_restart_init must be abort, zero or alg");
        }
        return v;
    }();
    return cached;
}

//  CAMR.ps_sigma_rnuc (1e-5 m): the nucleation radius of Sigma's
//  production term, dSigma = (3/r_nuc) d_alpha_seed per committed
//  nucleation event (docs/PLAN_sigma.md S-c).  Read only at
//  CAMR.ps_sigma >= 1.  Must be > 0.
inline amrex::Real ps_sigma_rnuc()
{
    static const amrex::Real cached = []() -> amrex::Real {
        const amrex::Real v = ps_dial_real("ps_sigma_rnuc", amrex::Real(1.0e-5));
        if (!(v > amrex::Real(0.0))) {
            amrex::Abort("CAMR.ps_sigma_rnuc must be > 0");
        }
        return v;
    }();
    return cached;
}

//  CAMR.ps_sigma_corridor (1) and CAMR.ps_sigma_bc_mirror (1): TEST-ONLY
//  (the B7/FF diagnostic, docs/PLAN_sigma.md §II.5-D; retired in commit 6).
//  Read only at CAMR.ps_sigma = 2.  ps_sigma_corridor = 0: the corridor
//  transfer counts but does not commit, as at level 1.  ps_sigma_bc_mirror
//  = 0: the relaxed-boundary trial runs without Sigma.
inline int ps_sigma_test_switch(const char* key)
{
    const int v = ps_dial_int(key, 1);
    if (v != 0 && v != 1) { amrex::Abort(std::string("CAMR.") + key + " must be 0 or 1"); }
    return v;
}
inline int ps_sigma_corridor()
{
    static const int cached = ps_sigma_test_switch("ps_sigma_corridor");
    return cached;
}
inline int ps_sigma_bc_mirror()
{
    static const int cached = ps_sigma_test_switch("ps_sigma_bc_mirror");
    return cached;
}

//  CAMR.ps_sigma_growth (1): TEST-ONLY (retired in commit 6).  0 switches
//  off Sigma's growth source (P1) and nothing else, so the production
//  identity of witness G3b can be asserted on its own.  Read only at
//  CAMR.ps_sigma >= 1.
inline int ps_sigma_growth()
{
    static const int cached = []() -> int {
        const int v = ps_dial_int("ps_sigma_growth", 1);
        if (v != 0 && v != 1) {
            amrex::Abort("CAMR.ps_sigma_growth must be 0 or 1");
        }
        return v;
    }();
    return cached;
}

#endif  // USE_PS_HYDRO
#endif  // PS_DIALS_H
