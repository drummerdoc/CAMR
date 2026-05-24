#include "CAMR.H"
#include "Hydro.H"
#include "Hydro_ctoprim.H"
#include "CAMR_Constants.H"

using namespace amrex;

// Device helper to apply the van Leer limiter
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real van_leer_limiter(amrex::Real sl_L, amrex::Real sl_R) {
    if (sl_L * sl_R <= 0.0) {
        return 0.0;
    }
    return (2.0 * sl_L * sl_R) / (sl_L + sl_R);
}

struct SplitFluxContribution {
    amrex::Real wave_flux_L[3]; // Flux mapped to Left grid cell
    amrex::Real wave_flux_R[3]; // Flux mapped to Right grid cell
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SplitFluxContribution evaluate_bct_sonic_split(
    const amrex::GpuArray<amrex::Real, 3>& q_L,  // Primitive vector [rho, u, p] Left
    const amrex::GpuArray<amrex::Real, 3>& q_R,  // Primitive vector [rho, u, p] Right
    amrex::Real c_L, amrex::Real c_R,
    int wave_index) // wave_index: 1 = (u-c), 2 = (u+c)
{
    SplitFluxContribution sfc = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Calculate wave speeds at both left and right face states
    amrex::Real lambda_L = (wave_index == 1) ? (q_L[1] - c_L) : (q_L[1] + c_L);
    amrex::Real lambda_R = (wave_index == 2) ? (q_R[1] - c_R) : (q_R[1] + c_R);

    // Compute characteristic wave strength vector via Left Eigenvectors
    // (Assuming compute_bct_eigenvectors is defined as shown previously)
    auto mats_L = EOS::compute_bct_eigenvectors(q_L[0], q_L[1], c_L);
    
    amrex::Real dq[3] = { q_R[0] - q_L[0], q_R[1] - q_L[1], q_R[2] - q_L[2] };
    amrex::Real alpha_k = 0.0;
    for(int m=0; m<3; ++m) {
        alpha_k += mats_L.L[wave_index][m] * dq[m];
    }

    // CHECK FOR SONIC TRANSIT (Sign changes across the face)
    if (lambda_L * lambda_R < 0.0) {
        // Linearized sonic state approximation to locate where lambda == 0
        amrex::Real theta = lambda_L / (lambda_L - lambda_R); // Interpolation parameter
        
        amrex::Real q_sonic[3];
        for(int m=0; m<3; ++m) {
            q_sonic[m] = q_L[m] + theta * (q_R[m] - q_L[m]);
        }

        // Decompose the wave step across the zero boundary
        // Signal 1: From Left to Sonic State
        amrex::Real alpha_k_Left = theta * alpha_k;
        // Signal 2: From Sonic State to Right
        amrex::Real alpha_k_Right = (1.0 - theta) * alpha_k;

        // Map fluxes conditionally based on the orientation of the entry/exit speeds
        if (lambda_L > 0.0) {
            for(int m=0; m<3; ++m) {
                sfc.wave_flux_L[m] += alpha_k_Left * mats_L.R[m][wave_index];
                sfc.wave_flux_R[m] += alpha_k_Right * mats_L.R[m][wave_index]; 
            }
        } else {
            for(int m=0; m<3; ++m) {
                sfc.wave_flux_R[m] += alpha_k_Left * mats_L.R[m][wave_index];
                sfc.wave_flux_L[m] += alpha_k_Right * mats_L.R[m][wave_index];
            }
        }
    } 
    else {
        // NO SONIC POINT: Pure upwind assignment based on sign
        amrex::Real wave_speed_avg = 0.5 * (lambda_L + lambda_R);
        if (wave_speed_avg >= 0.0) {
            for(int m=0; m<3; ++m) sfc.wave_flux_L[m] += alpha_k * mats_L.R[m][wave_index];
        } else {
            for(int m=0; m<3; ++m) sfc.wave_flux_R[m] += alpha_k * mats_L.R[m][wave_index];
        }
    }
    return sfc;
}

void
CAMR::construct_hydro_source (const MultiFab& S,
                              MultiFab& src_to_fill,
                              Real /*time*/,
                              Real dt)
{
    src_to_fill.setVal(0);

    if (verbose) {
        if (do_mol) {
            amrex::Print() << "... Computing MOL-based hydro advance" << std::endl;
        } else {
            amrex::Print() << "... Computing Godunov-based hydro advance" << std::endl;
        }
    }

    Real fac_for_reflux = (do_mol) ? Real(0.5) : Real(1.0);

    AMREX_ASSERT(S.nGrow() == numGrow());

    // Fill the source terms to go into the hydro with only the old-time sources
    int ng = 0;

    for (int n = 0; n < src_list.size(); ++n) {
        MultiFab::Saxpy(sources_for_hydro, 1.0, *old_sources[src_list[n]], 0, 0, NVAR, ng);
    }
    sources_for_hydro.FillBoundary(geom.periodicity());

    int finest_level = parent->finestLevel();

    const auto& dx    = geom.CellSizeArray();

    const PassMap* lpmap = d_pass_map;

    Real dx1 = dx[0];
    for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
      dx1 *= dx[dir];
    }

    std::array<Real, AMREX_SPACEDIM> dxD = {
      {AMREX_D_DECL(dx1, dx1, dx1)}};
    const Real* dxDp = &(dxD[0]);

    MultiFab& S_new = get_new_data(State_Type);

    BL_PROFILE_VAR("CAMR::advance_hydro_umdrv()", PC_UMDRV);

#ifdef USE_PR_EOS


#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
      amrex::MFItInfo tiling = amrex::TilingIfNotGPU() ? amrex::MFItInfo().EnableTiling(hydro_tile_size) : amrex::MFItInfo();
      for (MFIter mfi(S_new, tiling); mfi.isValid(); ++mfi)
		{
		  const Box& bx = mfi.tilebox();
		  const Box& qbx = amrex::grow(bx, numGrow());

		  amrex::GpuArray<amrex::FArrayBox, AMREX_SPACEDIM> flux;
		  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
			const Box& efbx = amrex::surroundingNodes(bx, dir);
			flux[dir].resize(efbx, NVAR, amrex::The_Async_Arena());
			flux[dir].setVal<RunOn::Device>(0.);
		  }

		  auto const& sarr    = S.const_array(mfi);
		  auto const& hyd_src = src_to_fill.array(mfi);

		  // Resize Temporary Fabs
		  FArrayBox q(qbx, QVAR, amrex::The_Async_Arena());
		  FArrayBox qaux(qbx, NQAUX, amrex::The_Async_Arena());
		  FArrayBox src_q(qbx, QVAR, amrex::The_Async_Arena());

		  // Get Arrays to pass to the gpu.
		  auto const& qarr    = q.array();
		  auto const& qauxar  = qaux.array();
		  auto const& srcqarr = src_q.array();

		  BL_PROFILE_VAR("ctoprim()", ctop);
		  const Real small_num        = CAMRConstants::small_num;
		  const Real dual_energy_eta  = CAMR::dual_energy_eta1;
		  int l_allow_negative_energy = CAMR::allow_negative_energy;
		  ParallelFor(
            qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  hydro_ctoprim(i, j, k, sarr, qarr, qauxar, *lpmap,
								small_num, dual_energy_eta, l_allow_negative_energy);
            });
		  BL_PROFILE_VAR_STOP(ctop);

		  const GpuArray<const Array4<      Real>, AMREX_SPACEDIM>
			flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};
		  const amrex::GpuArray<const Array4<const Real>, AMREX_SPACEDIM>
			a{{AMREX_D_DECL(area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};

		  // Create source terms for primitive variables
		  if (!do_mol) {
            const auto& src_in = sources_for_hydro.array(mfi);
            ParallelFor(
			  qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
				hydro_srctoprim(i, j, k, qarr, qauxar, src_in, srcqarr, *lpmap);
              });
		  }


		  // Parallel thread launch over the grid face interfaces (i+1/2)
		  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		  {
            // ---------------------------------------------------------------
            // 1. EXTRACT STENCIL VARIABLES [-2, -1, 0, 1, 2] RELATIVE TO INTERFACE
            //    Mapping: idx 0=(i-2), 1=(i-1), 2=(i), 3=(i+1), 4=(i+2)
            // ---------------------------------------------------------------
 		    amrex::GpuArray<amrex::Real, 5> r, u, v, w, e, T, p, c;
            
            int idx = 0;
            for (int s = -2; s <= 2; ++s) {
			  r[idx] = qarr(i, j, k, QRHO);
			  AMREX_D_TERM(u[idx] = qarr(i, j, k, QU);,
						   v[idx] = qarr(i, j, k, QV);,
						   w[idx] = qarr(i, j, k, QW););
			  e[idx] = qarr(i,j,k,QREINT) / r[idx];
			  T[idx] = qarr(i, j, k, QTEMP);
			  p[idx] = qarr(i, j, k, QPRES);
			  c[idx] = qauxar(i, j, k, QC);
			  idx++;
            }

            // ---------------------------------------------------------------
            // 2. CHARACTERISTIC SPACE LIMITING FOR CELL (i) [idx 2]
            // ---------------------------------------------------------------
            auto mats_i = EOS::compute_bct_eigenvectors(r[2], u[2], c[2]);
            
            amrex::GpuArray<amrex::Real, 3> dq_minus_i = { r[2] - r[1], u[2] - u[1], p[2] - p[1] };
            amrex::GpuArray<amrex::Real, 3> dq_plus_i  = { r[3] - r[2], u[3] - u[2], p[3] - p[2] };
            amrex::GpuArray<amrex::Real, 3> dq_limited_i = { 0.0, 0.0, 0.0 };

            for (int w = 0; w < 3; ++w) {
                amrex::Real dw_minus
				  = mats_i.L[w][0]*dq_minus_i[0]
				  + mats_i.L[w][1]*dq_minus_i[1]
				  + mats_i.L[w][2]*dq_minus_i[2];
                amrex::Real dw_plus
				  = mats_i.L[w][0]*dq_plus_i[0]
				  + mats_i.L[w][1]*dq_plus_i[1]
				  + mats_i.L[w][2]*dq_plus_i[2];

                amrex::Real dw_lim   = van_leer_limiter(dw_minus, dw_plus);
				
                for (int m = 0; m < 3; ++m) {
                    dq_limited_i[m] += dw_lim * mats_i.R[m][w];
                }
            }
			
            // ---------------------------------------------------------------
            // 3. CHARACTERISTIC SPACE LIMITING FOR CELL (i+1) [idx 3]
            // ---------------------------------------------------------------
            auto mats_ip1 = EOS::compute_bct_eigenvectors(r[3], u[3], c[3]);
            
            amrex::GpuArray<amrex::Real, 3> dq_minus_ip1 = { r[3] - r[2], u[3] - u[2], p[3] - p[2] };
            amrex::GpuArray<amrex::Real, 3> dq_plus_ip1  = { r[4] - r[3], u[4] - u[3], p[4] - p[3] };
            amrex::GpuArray<amrex::Real, 3> dq_limited_ip1 = { 0.0, 0.0, 0.0 };

            for (int w = 0; w < 3; ++w) {
                amrex::Real dw_minus
				  = mats_ip1.L[w][0]*dq_minus_ip1[0]
				  + mats_ip1.L[w][1]*dq_minus_ip1[1]
				  + mats_ip1.L[w][2]*dq_minus_ip1[2];
                amrex::Real dw_plus
				  = mats_ip1.L[w][0]*dq_plus_ip1[0]
				  + mats_ip1.L[w][1]*dq_plus_ip1[1]
				  + mats_ip1.L[w][2]*dq_plus_ip1[2];

                amrex::Real dw_lim   = van_leer_limiter(dw_minus, dw_plus);

                for (int m = 0; m < 3; ++m) {
                    dq_limited_ip1[m] += dw_lim * mats_ip1.R[m][w];
                }
            }

            // ---------------------------------------------------------------
            // 4. RESOLVE FINAL RECONSTRUCTED INTERFACE EDGE STATES (i+1/2)
            // ---------------------------------------------------------------
            amrex::Real rho_L   = r[2] + 0.5 * dq_limited_i[0];
            amrex::Real u_L     = u[2] + 0.5 * dq_limited_i[1];
            amrex::Real p_L     = p[2] + 0.5 * dq_limited_i[2];
            amrex::Real v_L     = v[2]; // Transverse components tracked via upwind
            amrex::Real w_L     = w[2];

            amrex::Real rho_R   = r[3] - 0.5 * dq_limited_ip1[0];
            amrex::Real u_R     = u[3] - 0.5 * dq_limited_ip1[1];
            amrex::Real p_R     = p[3] - 0.5 * dq_limited_ip1[2];
            amrex::Real v_R     = v[3];
            amrex::Real w_R     = w[3];

            // Reconstruct internal energy from limited edge states
            amrex::Real T_L_est = T[2]; amrex::Real e_int_L = e[2];
            amrex::Real T_R_est = T[3]; amrex::Real e_int_R = e[3];
            EOS::apply_spinodal_guard(fluid,rho_L, e_int_L, T_L_est);
            EOS::apply_spinodal_guard(fluid,rho_R, e_int_R, T_R_est);

			amrex::Real c_L, chi_L, kappa_L;
			amrex::Real c_R, chi_R, kappa_R;
			RET2C(fluid, e_int_L, rho_L, T_L_est, c_L, p_L, chi_L, kappa_L);
			RET2C(fluid, e_int_R, rho_R, T_R_est, c_R, p_R, chi_R, kappa_R);
			
            // ---------------------------------------------------------------
            // 5. ASSEMBLE RIEMANN INTERFACE REFERENCE MATRICES
            // ---------------------------------------------------------------
            amrex::Real rho_avg   = 0.5 * (rho_L + rho_L);
            amrex::Real u_avg     = 0.5 * (u_L + u_R);
            amrex::Real c_avg     = 0.5 * (c_L + c_R);
            amrex::Real e_avg     = 0.5 * (e_int_L + e_int_R);
            amrex::Real p_avg     = 0.5 * (p_L + p_R);
            amrex::Real chi_avg   = 0.5 * (chi_L + chi_R);
            amrex::Real kappa_avg = 0.5 * (kappa_L + kappa_R);

            auto mats = EOS::compute_bct_eigenvectors(rho_avg, u_avg, c_avg);

            amrex::GpuArray<amrex::Real, 3> dq = { rho_R - rho_L, u_R - u_L, p_R - p_L };

            // ---------------------------------------------------------------
            // 6. SOLVE BCT CHARACTERISTIC WAVE PROPAGATION
            // ---------------------------------------------------------------
			amrex::GpuArray<amrex::Real, 3> lambda = { u_avg, u_avg + c_avg, u_avg - c_avg };
			amrex::GpuArray<amrex::Real, 3> F_prim_L = { 0.0, 0.0, 0.0 };
			amrex::GpuArray<amrex::Real, 3> F_prim_R = { 0.0, 0.0, 0.0 };
			for (int w = 0; w < 3; ++w) {
			  amrex::Real alpha_w = mats.L[w][0]*dq[0] + mats.L[w][1]*dq[1] + mats.L[w][2]*dq[2];
			  amrex::Real lambda_L = (w == 0) ? u_L : ((w == 1) ? (u_L + c_L) : (u_L - c_L));
			  amrex::Real lambda_R = (w == 0) ? u_R : ((w == 1) ? (u_R + c_R) : (u_R - c_R));
			  if (lambda_L * lambda_R < 0.0 && w > 0) {
				amrex::GpuArray<amrex::Real, 3> q_L = { rho_L, u_L, p_L };
				amrex::GpuArray<amrex::Real, 3> q_R = { rho_R, u_R, p_R };
				auto split = evaluate_bct_sonic_split(q_L, q_R, c_L, c_R, w);
				for (int m = 0; m < 3; ++m) {
				  F_prim_L[m] += split.wave_flux_L[m];
				  F_prim_R[m] += split.wave_flux_R[m];
				}
			  } else {
				if (lambda[w] >= 0.0) {
				  for (int m = 0; m < 3; ++m) {
					F_prim_L[m] += alpha_w * mats.R[m][w];
				  }
				} else {
				  for (int m = 0; m < 3; ++m) {
					F_prim_R[m] += alpha_w * mats.R[m][w];
				  }
				}
			  }
			}

			// ---------------------------------------------------------------
			// 7. TRANSFORM PRIMITIVE JUMPS INTO CONSERVATIVE FLUXES
			// ---------------------------------------------------------------
			amrex::Real E_avg = rho_avg * e_avg + 0.5 * rho_avg * (u_avg*u_avg + v_L*v_L + w_L*w_L);
			// Compute exact total energy flux components via chain rules
			amrex::Real E_flux_L = F_prim_L[0] * (E_avg + p_avg + rho_avg * u_avg * u_avg)
			  + F_prim_L[1] * u_avg * (e_avg + 0.5 * u_avg * u_avg - (chi_avg / kappa_avg))
			  + F_prim_L[2] * u_avg * (1.0 + (rho_avg / kappa_avg));
			amrex::Real E_flux_R = F_prim_R[0] * (E_avg + p_avg + rho_avg * u_avg * u_avg)
			  + F_prim_R[1] * u_avg * (e_avg + 0.5 * u_avg * u_avg - (chi_avg / kappa_avg))
			  + F_prim_R[2] * u_avg * (1.0 + (rho_avg / kappa_avg));
			amrex::GpuArray<amrex::Real, 5> base_F_L;
			base_F_L[URHO]  = F_prim_L[0] * u_avg + rho_avg * F_prim_L[1];
			// Mass
			base_F_L[UMX]   = base_F_L[URHO] * u_avg + rho_avg * u_avg * F_prim_L[1];
			// Normal Mom
			base_F_L[UMY]   = 0.0;
			// Placeholders
			//base_F_L[UMZ]   = 0.0;
			base_F_L[UEDEN] = E_flux_L;
			// Energy
			amrex::GpuArray<amrex::Real, 5> base_F_R;
			base_F_R[URHO]  = F_prim_R[0] * u_avg + rho_avg * F_prim_R[1];
			base_F_R[UMX]   = base_F_R[URHO] * u_avg + rho_avg * u_avg * F_prim_R[1];
			base_F_R[UMY]   = 0.0;
			//base_F_R[UMZ]   = 0.0;
			base_F_R[UEDEN] = E_flux_R;

			// Upwind passive transverse momentum variables based on advection velocity
			if (lambda[0] >= 0.0) {
			  base_F_L[UMY] = base_F_L[URHO] * v_L;
			  //base_F_L[UMZ] = base_F_L[URHO] * w_L;
			} else {
			  base_F_R[UMY] = base_F_R[URHO] * v_R;
			  //base_F_R[UMZ] = base_F_R[URHO] * w_R;
			}

			// ----------------------------
			// 8. WRITE CONSERVATIVE FLUXES
			// ----------------------------
			/*
			Fxarr(i, j, k, URHO)  = base_F_L[URHO]  + base_F_R[URHO];
			Fxarr(i, j, k, UMX)   = base_F_L[UMX]   + base_F_R[UMX] + 0.5 * (p_L + p_R);
			Fxarr(i, j, k, UMY)   = base_F_L[UMY]   + base_F_R[UMY];
			//Fxarr(i, j, k, UMZ)   = base_F_L[UMZ]   + base_F_R[UMZ];
			Fxarr(i, j, k, UEDEN) = base_F_L[UEDEN] + base_F_R[UEDEN];
			*/
		});
		
		//
		// Here fac_for_reflux = 1.0 if doing Godunov, 0.5 if doing MOL
		//
		if (do_reflux) {
		  if (level < finest_level) {
			getFluxReg(level + 1).CrseAdd(mfi,
              {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
			  dxDp, fac_for_reflux*dt, amrex::RunOn::Device);
		  }
		  if (level > 0) {
			getFluxReg(level).FineAdd(mfi,
              {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
			  dxDp, fac_for_reflux*dt, amrex::RunOn::Device);
		  }
		} // do_reflux
		} // mfi
    } // openmp

#else
	
#ifdef AMREX_USE_EB
    const auto& ebfact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
#ifdef AMREX_USE_EB
      int ncomp = src_to_fill.nComp();
      FArrayBox dm_as_fine(Box::TheUnitBox(),ncomp);
      FArrayBox fab_drho_as_crse(Box::TheUnitBox(),ncomp);
      IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
#endif

      amrex::MFItInfo tiling = amrex::TilingIfNotGPU() ? amrex::MFItInfo().EnableTiling(hydro_tile_size) : amrex::MFItInfo();
      for (MFIter mfi(S_new, tiling); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
        EBCellFlagFab const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        auto const& flag_arr = flagfab.const_array();

        if (flagfab.getType(bx) != FabType::covered) {
            auto const& vfrac_arr = volfrac->const_array(mfi);
#endif
            const Box& qbx = amrex::grow(bx, numGrow());

            amrex::GpuArray<amrex::FArrayBox, AMREX_SPACEDIM> flux;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              const Box& efbx = amrex::surroundingNodes(bx, dir);
              flux[dir].resize(efbx, NVAR, amrex::The_Async_Arena());
              flux[dir].setVal<RunOn::Device>(0.);
            }

            auto const& sarr    = S.const_array(mfi);
            auto const& hyd_src = src_to_fill.array(mfi);

            // Resize Temporary Fabs
            FArrayBox q(qbx, QVAR, amrex::The_Async_Arena());
            FArrayBox qaux(qbx, NQAUX, amrex::The_Async_Arena());
            FArrayBox src_q(qbx, QVAR, amrex::The_Async_Arena());

            // Get Arrays to pass to the gpu.
            auto const& qarr    = q.array();
            auto const& qauxar  = qaux.array();
            auto const& srcqarr = src_q.array();

            BL_PROFILE_VAR("ctoprim()", ctop);
            const Real small_num        = CAMRConstants::small_num;
            const Real dual_energy_eta  = CAMR::dual_energy_eta1;
            int l_allow_negative_energy = CAMR::allow_negative_energy;
            ParallelFor(
              qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef AMREX_USE_EB
                if (!flag_arr(i,j,k).isCovered()) {
#endif
                    hydro_ctoprim(i, j, k, sarr, qarr, qauxar, *lpmap,
                                  small_num, dual_energy_eta, l_allow_negative_energy);
#ifdef AMREX_USE_EB
                } else {
                   for (int n=0; n<QVAR; n++) qarr(i,j,k,n) = 0.;
                }
#endif
            });
            BL_PROFILE_VAR_STOP(ctop);

            const GpuArray<const Array4<      Real>, AMREX_SPACEDIM>
              flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};
            const amrex::GpuArray<const Array4<const Real>, AMREX_SPACEDIM>
              a{{AMREX_D_DECL(area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};

        // Create source terms for primitive variables
        if (!do_mol) {
            const auto& src_in = sources_for_hydro.array(mfi);
            ParallelFor(
              qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef AMREX_USE_EB
                if (!flag_arr(i,j,k).isCovered()) {
#endif
                    hydro_srctoprim(i, j, k, qarr, qauxar, src_in, srcqarr, *lpmap);
#ifdef AMREX_USE_EB
                } else {
                   for (int n=0; n<QVAR; n++) srcqarr(i,j,k,n) = 0.;
                }
#endif
              });
        }

#ifdef AMREX_USE_EB
        int ngrow_bx;
        if (redistribution_type == "StateRedist") {
           ngrow_bx = 3;
        } else {
           ngrow_bx = 2;
        }
        const Box& bxg_i  = grow(bx,ngrow_bx);
        if (flagfab.getType(bxg_i) != FabType::regular) {

            EBFluxRegister* fr_as_crse = nullptr;
            if (do_reflux && level < parent->finestLevel()) {
                CAMR& fine_level = getLevel(level+1);
                fr_as_crse = &fine_level.flux_reg;
            }

            EBFluxRegister* fr_as_fine = nullptr;
            if (do_reflux && level > 0) {
                fr_as_fine = &flux_reg;
            }

            int as_crse = (fr_as_crse != nullptr);
            int as_fine = (fr_as_fine != nullptr);

            FArrayBox* p_drho_as_crse = (fr_as_crse) ?
                    fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
            const IArrayBox* p_rrflag_as_crse = (fr_as_crse) ?
                   fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;

            if (fr_as_fine) {
                const Box dbox1 = geom.growPeriodicDomain(1);
                Box bx_for_dm(amrex::grow(bx,1) & dbox1);
                dm_as_fine.resize(bx_for_dm,ncomp);
                dm_as_fine.setVal<RunOn::Device>(0.0);
            }

            const amrex::StateDescriptor* desc = state[State_Type].descriptor();
            const auto& bcs = desc->getBCs();
            amrex::Gpu::DeviceVector<amrex::BCRec> bcs_d(desc->nComp());
            amrex::Gpu::copy(
              amrex::Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_d.begin());

            const auto& dxInv = geom.InvCellSizeArray();

            // Return hyd_src - centered at half-time if using Godunov method
            //                - centered at  old-time if using MOL method
            //
            // The dt we pass in here is used if (do_mol == 0), i.e.
            //      in the Godunov prediction, but also if we do StateRedistribution
            //
            hydro_umdrv_eb(do_mol, bx, bxg_i, mfi, geom, &ebfact,
                          phys_bc.lo(), phys_bc.hi(),
                          sarr, hyd_src, qarr, qauxar, srcqarr,
                          vfrac_arr, flag_arr, dx, dxInv, flx_arr, a,
                          as_crse, p_drho_as_crse->array(), p_rrflag_as_crse->array(),
                          as_fine, dm_as_fine.array(), level_mask.const_array(mfi),
                          dt, ppm_type, plm_iorder, use_pslope,
                          use_flattening, transverse_reset_density,
                          small, small_dens, small_pres, CAMRConstants::smallu, difmag,
                          bcs_d.data(), redistribution_type, lpmap, eb_weights_type);

            //
            // Here fac_for_reflux = 1.0 if doing Godunov, 0.5 if doing MOL
            //
            if (do_reflux) {
                if (level < finest_level) {
                     getFluxReg(level + 1).CrseAdd(mfi,
                        {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                        dxDp, fac_for_reflux*dt, (*volfrac)[mfi],
                        {AMREX_D_DECL(&(*areafrac[0])[mfi], &(*areafrac[1])[mfi], &(*areafrac[2])[mfi])},
                        amrex::RunOn::Device);
                }
                if (level > 0) {
                    getFluxReg(level).FineAdd(mfi,
                       {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                       dxDp, fac_for_reflux*dt, (*volfrac)[mfi],
                       {AMREX_D_DECL(&(*areafrac[0])[mfi], &(*areafrac[1])[mfi], &(*areafrac[2])[mfi])},
                       dm_as_fine, amrex::RunOn::Device);
                } // level > 0
            } // do_reflux
        } else {
#endif
            // Return hyd_src - centered at half-time if using Godunov method
            //                - centered at  old-time if using MOL method
            //
            // Note that the dt here is only used if (do_mol == 0), i.e.
            //      in the Godunov prediction
            //
            hydro_umdrv(do_mol, bx, geom, phys_bc.lo(), phys_bc.hi(),
                       sarr, hyd_src, qarr, qauxar, srcqarr, dx,
                       dt, ppm_type, plm_iorder, use_pslope,
                       use_flattening, transverse_reset_density,
                       small, small_dens, small_pres, CAMRConstants::smallu, difmag,
                       flx_arr, a, volume.array(mfi), lpmap);

            //
            // Here fac_for_reflux = 1.0 if doing Godunov, 0.5 if doing MOL
            //
            if (do_reflux) {
                if (level < finest_level) {
                    getFluxReg(level + 1).CrseAdd(mfi,
                        {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                        dxDp, fac_for_reflux*dt, amrex::RunOn::Device);
                }
                if (level > 0) {
                    getFluxReg(level).FineAdd(mfi,
                       {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                       dxDp, fac_for_reflux*dt, amrex::RunOn::Device);
                }
            } // do_reflux

#ifdef AMREX_USE_EB
        } // regular
#endif

#ifdef AMREX_USE_EB
        } // not covered
#endif
      } // mfi
    } // openmp

#endif //PR_EOS

    BL_PROFILE_VAR_STOP(PC_UMDRV);
}
