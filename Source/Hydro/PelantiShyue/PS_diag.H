#ifndef PS_DIAG_H
#define PS_DIAG_H
// =====================================================================
//  PS_diag.H — host-only diagnostic dials of the six-equation
//  Pelanti–Shyue model on the CAMR grid that are read from more than one
//  report site.  Every dial defaults off and is documented once at its
//  accessor.  Included through the PS_relaxation.H umbrella.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "PS_dials.H"      // ps_dial_int (CAMR.* reads)

// CAMR.ps_diag_alpha (0): report the alpha_1 range at labelled points
// (CAMR_advance.cpp, CAMR.cpp) so a drive toward a pure phase can be
// attributed to mass transfer, advection/C-F, interpolation or reflux.
inline int ps_diag_alpha()
{
    static const int v = ps_dial_int("ps_diag_alpha", 0);
    return v;
}

#endif // USE_PS_HYDRO
#endif // PS_DIAG_H
