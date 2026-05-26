#ifndef Time_Centered_CTU_Fluxes_3D_H
#define Time_Centered_CTU_Fluxes_3D_H

#include "AMReX_FArrayBox.H"

using namespace amrex;

void Time_Centered_CTU_Fluxes_3D(amrex::FArrayBox& prim, const amrex::Box& bx,
                                 amrex::GpuArray<amrex::FArrayBox, AMREX_SPACEDIM>& flux,
                                 amrex::Real dt, const GpuArray<double, AMREX_SPACEDIM>& dx);

#endif
