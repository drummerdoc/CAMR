#ifndef GODUNOV_UTILS_1D_H
#define GODUNOV_UTILS_1D_H

// =====================================================================
//  Godunov_utils_1D.H  —  1D STUB.
//
//  The Godunov (directionally-split PLM/PPM) hydro path is NOT
//  implemented for AMREX_SPACEDIM==1.  A DIM=1 build is only supported
//  under USE_PS_HYDRO, where the Pelanti-Shyue wave-propagation solver
//  (PS_umeth) replaces the Godunov path entirely, so none of the 2D/3D
//  transverse-correction helpers that live in Godunov_utils_2D/3D.H are
//  needed.  This header exists so the Make.package
//  `CEXE_headers += Godunov_utils_$(DIM)D.H` dependency resolves in 1D;
//  it intentionally declares nothing.
// =====================================================================

#endif
