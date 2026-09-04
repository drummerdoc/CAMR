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

#endif  // USE_PS_HYDRO
#endif  // PS_DIALS_H
