// =====================================================================
//  gerg_coeffs.cpp — the one definition of the GERG-2008 coefficient and
//  saturation-spline arrays on a GPU build (device constant memory);
//  every other translation unit sees extern declarations (gerg_co2.H).
//  On a host build the arrays are header constexpr and this file is inert.
// =====================================================================
#define GERG_DEFINE_COEFFS
#include "gerg_co2_sat.H"
