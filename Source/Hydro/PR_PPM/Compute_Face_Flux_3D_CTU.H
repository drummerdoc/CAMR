#ifndef Compute_Face_Flux_3D_CTU_H
#define Compute_Face_Flux_3D_CTU_H

#include "IndexDefines.H"
#include "Perform_1D_PPM_Time_Centered_Reconstruction.H"
#include "solve_exact_pr_hem_riemann.H"
#include "compute_pr_speed_of_sound.H"
#include "compute_phase_aware_viscosity_x.H"
#include "CAMR.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Compute_Face_Flux_3D_CTU(int i, int j, int k,
                              int di, int dj, int dk,          // Normal vector strides (e.g., 1,0,0)
                              int t1_di, int t1_dj, int t1_dk, // Transverse Direction 1 (e.g., 0,1,0) - if all zeros, effectively skips correction
                              int t2_di, int t2_dj, int t2_dk, // Transverse Direction 2 (e.g., 0,0,1) - if all zeros, effectively skips correction
                              int u_idx, int v_idx, int w_idx, // Velocity component index mappings
                              amrex::Array4<const amrex::Real> const& prim,
                              amrex::Array4<const amrex::Real> const& trans_flux_1,
                              amrex::Array4<const amrex::Real> const& trans_prim_1,
                              amrex::Array4<const amrex::Real> const& trans_flux_2,
                              amrex::Array4<const amrex::Real> const& trans_prim_2,
                              amrex::Array4<amrex::Real> const& gd_flux,
                              amrex::Array4<amrex::Real> const& gd_state,
                              amrex::Real dt, const amrex::GpuArray<amrex::Real, 3>& dx)
{
    using namespace amrex;

    // =========================================================================
    // 1. Time-Centered 1D PPM Normal Reconstructions onto Interface Edges
    // =========================================================================
    Real rho_L,  rho_R;
    Real p_L,    p_R;
    Real u_L,    u_R; // Normal velocity
    Real v_L,    v_R; // Tangential velocity 1
    Real w_L,    w_R; // Tangential velocity 2
    Real Y_L,    Y_R;

    Real dx_normal = dx[u_idx - UMX];
    Real dx_trans1 = dx[v_idx - UMX];
    Real dx_trans2 = dx[w_idx - UMX];

    Perform_1D_PPM_Time_Centered_Reconstruction(i, j, k, di, dj, dk, u_idx, prim, rho_L,  rho_R,  QRHO,   dt, dx_normal);
    Perform_1D_PPM_Time_Centered_Reconstruction(i, j, k, di, dj, dk, u_idx, prim, p_L,    p_R,    QPRES,  dt, dx_normal);
    Perform_1D_PPM_Time_Centered_Reconstruction(i, j, k, di, dj, dk, u_idx, prim, u_L,    u_R,    u_idx,  dt, dx_normal);
    Perform_1D_PPM_Time_Centered_Reconstruction(i, j, k, di, dj, dk, u_idx, prim, v_L,    v_R,    v_idx,  dt, dx_normal);
    Perform_1D_PPM_Time_Centered_Reconstruction(i, j, k, di, dj, dk, u_idx, prim, w_L,    w_R,    w_idx,  dt, dx_normal);
    Perform_1D_PPM_Time_Centered_Reconstruction(i, j, k, di, dj, dk, u_idx, prim, Y_L,    Y_R,    QY,     dt, dx_normal);

    // =========================================================================
    // 2. Fetch Central Cell Indicators for Non-Ideal Acoustic Corrections
    // =========================================================================
    Real rho_c = prim(i, j, k, QRHO);
    Real GAME  = prim(i, j, k, QGAME);

    Real cdtdx1 = 0.5 * dt / dx_trans1;
    Real cdtdx2 = 0.5 * dt / dx_trans2;

    // =========================================================================
    // 3. Complete Transverse Corrections Across All Prim Components
    // =========================================================================
    
    // --- DIRECTION 1 CORRECTIONS
    // Spatial differences from Transverse Flux Array 1 (Cell i is between Face i and i+1)
    Real dF1_rho = cdtdx1 * (trans_flux_1(i+t1_di, j+t1_dj, k+t1_dk, URHO ) - trans_flux_1(i, j, k, URHO ));
    Real dF1_um  = cdtdx1 * (trans_flux_1(i+t1_di, j+t1_dj, k+t1_dk, u_idx) - trans_flux_1(i, j, k, u_idx));
	Real dF1_vm  = cdtdx1 * (trans_flux_1(i+t1_di, j+t1_dj, k+t1_dk, v_idx) - trans_flux_1(i, j, k, v_idx));
    Real dF1_wm  = cdtdx1 * (trans_flux_1(i+t1_di, j+t1_dj, k+t1_dk, w_idx) - trans_flux_1(i, j, k, w_idx));
    Real dF1_Y   = cdtdx1 * (trans_flux_1(i+t1_di, j+t1_dj, k+t1_dk, UY)    - trans_flux_1(i, j, k, UY   ));

    // Spatial differences from Transverse Flux Array 2
    Real dF2_rho = cdtdx2 * (trans_flux_2(i+t2_di, j+t2_dj, k+t2_dk, URHO ) - trans_flux_2(i, j, k, URHO ));
    Real dF2_um  = cdtdx2 * (trans_flux_2(i+t2_di, j+t2_dj, k+t2_dk, u_idx) - trans_flux_2(i, j, k, u_idx));
	Real dF2_vm  = cdtdx2 * (trans_flux_2(i+t2_di, j+t2_dj, k+t2_dk, v_idx) - trans_flux_2(i, j, k, v_idx));
    Real dF2_wm  = cdtdx2 * (trans_flux_2(i+t2_di, j+t2_dj, k+t2_dk, w_idx) - trans_flux_2(i, j, k, w_idx));
    Real dF2_UY  = cdtdx2 * (trans_flux_2(i+t2_di, j+t2_dj, k+t2_dk, UY   ) - trans_flux_2(i, j, k, UY   ));

	// Cell i-1 is on the left side of Face i
	Real rho_L_cell = prim(i-di,j-dj,k-dk, QRHO );
    Real p_L_cell   = prim(i-di,j-dj,k-dk, QPRES);
    Real u_L_cell   = prim(i-di,j-dj,k-dk, u_idx);
    Real v_L_cell   = prim(i-di,j-dj,k-dk, v_idx);
    Real w_L_cell   = prim(i-di,j-dj,k-dk, w_idx);
    Real Y_L_cell   = prim(i-di,j-dj,k-dk, QY   );

    amrex::Real c_L_cell   = compute_pr_speed_of_sound(CAMR::fluid, rho_L_cell, p_L_cell, Y_L_cell); 
    amrex::Real gamma4_L   = (rho_L_cell * c_L_cell * c_L_cell) / p_L_cell;

    // --- TRANSVERSE FLUX DIVERGENCE ACROSS CELL (i-1, j) ---
    Real dF1_rho_L = trans_flux_1(i-di+t1_di, j-dj+t1_dj, k-dk+t1_dk, URHO)  - trans_flux_1(i-di, j-dj, k-dk, URHO );
    Real dF1_mx_L  = trans_flux_1(i-di+t1_di, j-dj+t1_dj, k-dk+t1_dk, u_idx) - trans_flux_1(i-di, j-dj, k-dk, u_idx);
    Real dF1_my_L  = trans_flux_1(i-di+t1_di, j-dj+t1_dj, k-dk+t1_dk, v_idx) - trans_flux_1(i-di, j-dj, k-dk, v_idx);
    Real dF1_mz_L  = trans_flux_1(i-di+t1_di, j-dj+t1_dj, k-dk+t1_dk, w_idx) - trans_flux_1(i-di, j-dj, k-dk, w_idx);
    Real dF1_rE_L  = trans_flux_1(i-di+t1_di, j-dj+t1_dj, k-dk+t1_dk, UEDEN) - trans_flux_1(i-di, j-dj, k-dk, UEDEN);
    Real dF1_Y_L   = trans_flux_1(i-di+t1_di, j-dj+t1_dj, k-dk+t1_dk, UY)    - trans_flux_1(i-di, j-dj, k-dk, UY   );

    Real dF2_rho_L = trans_flux_2(i-di+t2_di, j-dj+t2_dj, k-dk+t2_dk, URHO)  - trans_flux_2(i-di, j-dj, k-dk, URHO );
    Real dF2_mx_L  = trans_flux_2(i-di+t2_di, j-dj+t2_dj, k-dk+t2_dk, u_idx) - trans_flux_2(i-di, j-dj, k-dk, u_idx);
    Real dF2_my_L  = trans_flux_2(i-di+t2_di, j-dj+t2_dj, k-dk+t2_dk, v_idx) - trans_flux_2(i-di, j-dj, k-dk, v_idx);
    Real dF2_mz_L  = trans_flux_2(i-di+t2_di, j-dj+t2_dj, k-dk+t2_dk, w_idx) - trans_flux_2(i-di, j-dj, k-dk, w_idx);
    Real dF2_rE_L  = trans_flux_2(i-di+t2_di, j-dj+t2_dj, k-dk+t2_dk, UEDEN) - trans_flux_2(i-di, j-dj, k-dk, UEDEN);
    Real dF2_Y_L   = trans_flux_2(i-di+t2_di, j-dj+t2_dj, k-dk+t2_dk, UY)    - trans_flux_2(i-di, j-dj, k-dk, UY   );

	// --- TRANSFORM CONSERVATIVE DIVERGENCES TO PRIMITIVE CORRECTIONS ---
    // Left State Primitive Corrections via Jacobian-equivalent split form
    amrex::Real drho_dt_L = cdtdx1 * dF1_rho_L + cdtdx2 * dF2_rho_L;
    amrex::Real du_dt_L   = cdtdx1 * (dF1_mx_L - u_L_cell * dF1_rho_L) + cdtdx2 * (dF2_mx_L - u_L_cell * dF2_rho_L);
    amrex::Real dv_dt_L   = cdtdx1 * (dF1_my_L - v_L_cell * dF1_rho_L) + cdtdx2 * (dF2_my_L - v_L_cell * dF2_rho_L);
    amrex::Real dw_dt_L   = cdtdx1 * (dF1_mz_L - w_L_cell * dF1_rho_L) + cdtdx2 * (dF2_mz_L - w_L_cell * dF2_rho_L);
    amrex::Real dY_dt_L   = cdtdx1 * (dF1_Y_L  - Y_L_cell * dF1_rho_L) + cdtdx2 * (dF2_Y_L  - Y_L_cell * dF2_rho_L);

	// 3D Pressure Transformation including multi-directional kinetic coupling
    amrex::Real ke_flux_corr_y_L = u_L_cell * (dF1_mx_L - u_L_cell * dF1_rho_L) + v_L_cell * (dF1_my_L - v_L_cell * dF1_rho_L) + w_L_cell * (dF1_mz_L - w_L_cell * dF1_rho_L);
    amrex::Real ke_flux_corr_z_L = u_L_cell * (dF2_mx_L - u_L_cell * dF2_rho_L) + v_L_cell * (dF2_my_L - v_L_cell * dF2_rho_L) + w_L_cell * (dF2_mz_L - w_L_cell * dF2_rho_L);
    
    amrex::Real dp_dt_y_L = (gamma4_L - 1.0) * (dF1_rE_L - ke_flux_corr_y_L - p_L_cell * (dF1_my_L - v_L_cell * dF1_rho_L) / rho_L_cell);
    amrex::Real dp_dt_z_L = (gamma4_L - 1.0) * (dF2_rE_L - ke_flux_corr_z_L - p_L_cell * (dF2_mz_L - w_L_cell * dF2_rho_L) / rho_L_cell);
    amrex::Real dp_dt_L   = cdtdx1 * dp_dt_y_L + cdtdx2 * dp_dt_z_L;

	// =========================================================================
    // 3. APPLY 3D CTU TRANSVERSE CORRECTIONS TO INTERFACE PRIMITIVES
    // =========================================================================
    // Left side face update
	rho_L -= drho_dt_L;
	u_L -= du_dt_L / rho_L_cell;
	v_L -= dv_dt_L / rho_L_cell;
	w_L -= dw_dt_L / rho_L_cell;
	p_L -= dp_dt_L;
	Y_L -= dY_dt_L / rho_L_cell;

	// Cell i is on the right side of Face i
	Real rho_R_cell = prim(i   ,j   ,k   , QRHO );
    Real p_R_cell   = prim(i   ,j   ,k   , QPRES);
    Real u_R_cell   = prim(i   ,j   ,k   , u_idx);
    Real v_R_cell   = prim(i   ,j   ,k   , v_idx);
    Real w_R_cell   = prim(i   ,j   ,k   , w_idx);
    Real Y_R_cell   = prim(i   ,j   ,k   , QY   );

    amrex::Real c_R_cell   = compute_pr_speed_of_sound(CAMR::fluid, rho_L_cell, p_L_cell, Y_L_cell); // Can I use the pressure that is computed during this call?
    amrex::Real gamma4_R   = (rho_L_cell * c_L_cell * c_L_cell) / p_L_cell;

    // --- TRANSVERSE FLUX DIVERGENCE ACROSS CELL (i, j) ---
    Real dF1_rho_R = trans_flux_1(i   +t1_di, j   +t1_dj, k   +t1_dk, URHO)  - trans_flux_1(i   , j   , k   , URHO );
    Real dF1_mx_R  = trans_flux_1(i   +t1_di, j   +t1_dj, k   +t1_dk, u_idx) - trans_flux_1(i   , j   , k   , u_idx);
    Real dF1_my_R  = trans_flux_1(i   +t1_di, j   +t1_dj, k   +t1_dk, v_idx) - trans_flux_1(i   , j   , k   , v_idx);
    Real dF1_mz_R  = trans_flux_1(i   +t1_di, j   +t1_dj, k   +t1_dk, w_idx) - trans_flux_1(i   , j   , k   , w_idx);
    Real dF1_rE_R  = trans_flux_1(i   +t1_di, j   +t1_dj, k   +t1_dk, UEDEN) - trans_flux_1(i   , j   , k   , UEDEN);
    Real dF1_Y_R   = trans_flux_1(i   +t1_di, j   +t1_dj, k   +t1_dk, UY)    - trans_flux_1(i   , j   , k   , UY   );

    Real dF2_rho_R = trans_flux_2(i   +t2_di, j   +t2_dj, k   +t2_dk, URHO)  - trans_flux_2(i   , j   , k   , URHO );
    Real dF2_mx_R  = trans_flux_2(i   +t2_di, j   +t2_dj, k   +t2_dk, u_idx) - trans_flux_2(i   , j   , k   , u_idx);
    Real dF2_my_R  = trans_flux_2(i   +t2_di, j   +t2_dj, k   +t2_dk, v_idx) - trans_flux_2(i   , j   , k   , v_idx);
    Real dF2_mz_R  = trans_flux_2(i   +t2_di, j   +t2_dj, k   +t2_dk, w_idx) - trans_flux_2(i   , j   , k   , w_idx);
    Real dF2_rE_R  = trans_flux_2(i   +t2_di, j   +t2_dj, k   +t2_dk, UEDEN) - trans_flux_2(i   , j   , k   , UEDEN);
    Real dF2_Y_R   = trans_flux_2(i   +t2_di, j   +t2_dj, k   +t2_dk, UY)    - trans_flux_2(i   , j   , k   , UY   );

	// --- TRANSFORM CONSERVATIVE DIVERGENCES TO PRIMITIVE CORRECTIONS ---
    // Left State Primitive Corrections via Jacobian-equivalent split form
    amrex::Real drho_dt_R = cdtdx1 * dF1_rho_R + cdtdx2 * dF2_rho_R;
    amrex::Real du_dt_R   = cdtdx1 * (dF1_mx_R - u_R_cell * dF1_rho_R) + cdtdx2 * (dF2_mx_R - u_R_cell * dF2_rho_R);
    amrex::Real dv_dt_R   = cdtdx1 * (dF1_my_R - v_R_cell * dF1_rho_R) + cdtdx2 * (dF2_my_R - v_R_cell * dF2_rho_R);
    amrex::Real dw_dt_R   = cdtdx1 * (dF1_mz_R - w_R_cell * dF1_rho_R) + cdtdx2 * (dF2_mz_R - w_R_cell * dF2_rho_R);
    amrex::Real dY_dt_R   = cdtdx1 * (dF1_Y_R  - Y_R_cell * dF1_rho_R) + cdtdx2 * (dF2_Y_R  - Y_R_cell * dF2_rho_R);

	// 3D Pressure Transformation including multi-directional kinetic coupling
    amrex::Real ke_flux_corr_y_R = u_R_cell * (dF1_mx_R - u_R_cell * dF1_rho_R) + v_R_cell * (dF1_my_R - v_R_cell * dF1_rho_R) + w_R_cell * (dF1_mz_R - w_R_cell * dF1_rho_R);
    amrex::Real ke_flux_corr_z_R = u_R_cell * (dF2_mx_R - u_R_cell * dF2_rho_R) + v_R_cell * (dF2_my_R - v_R_cell * dF2_rho_R) + w_R_cell * (dF2_mz_R - w_R_cell * dF2_rho_R);
    
    amrex::Real dp_dt_y_R = (gamma4_R - 1.0) * (dF1_rE_R - ke_flux_corr_y_R - p_R_cell * (dF1_my_R - v_R_cell * dF1_rho_R) / rho_R_cell);
    amrex::Real dp_dt_z_R = (gamma4_R - 1.0) * (dF2_rE_R - ke_flux_corr_z_R - p_R_cell * (dF2_mz_R - w_R_cell * dF2_rho_R) / rho_R_cell);
    amrex::Real dp_dt_R   = cdtdx1 * dp_dt_y_R + cdtdx2 * dp_dt_z_R;

	// =========================================================================
    // 3. APPLY 3D CTU TRANSVERSE CORRECTIONS TO INTERFACE PRIMITIVES
    // =========================================================================
    // Left side face update
	rho_R -= drho_dt_R;
	u_R -= du_dt_R / rho_R_cell;
	v_R -= dv_dt_R / rho_R_cell;
	w_R -= dw_dt_R / rho_R_cell;
	p_R -= dp_dt_R;
	Y_R -= dY_dt_R / rho_R_cell;

	
    // =========================================================================
    // 4. Inject Phase-Aware Artificial Viscosity (Normal direction only)
    // =========================================================================
    Real q_visc_L = compute_phase_aware_viscosity_x(i   , j   , k   , u_idx, di, dj, dk, prim, dx_normal);
    Real q_visc_R = compute_phase_aware_viscosity_x(i+di, j+dj, k+dk, u_idx, di, dj, dk, prim, dx_normal);

    p_L += q_visc_L;
    p_R += q_visc_R;

    // =========================================================================
    // 5. Submit Fully Traced & Corrected States to exact Riemann Solver
    // =========================================================================
	amrex::Array1D<amrex::Real, 0, 5> W_L{rho_L, u_L, v_L, w_L, p_L, Y_L};
	amrex::Array1D<amrex::Real, 0, 5> W_R{rho_R, u_R, v_R, w_R, p_R, Y_R};
	auto riemann = solve_exact_pr_hem_riemann(CAMR::fluid, W_L, W_R);

    // Add Artificial Viscosity total energy heating dissipation matching 
    Real rho_face = riemann.state[QRHO];
    Real u_face   = riemann.state[QU];
    Real q_interface = 0.5 * (q_visc_L + q_visc_R);
    riemann.flux[UEDEN] += q_interface * u_face;

    // Output mapped riemann state and fluxes to final array indices
    gd_state(i, j, k, QRHO)   = riemann.state[QRHO];
    gd_state(i, j, k, u_idx)  = riemann.state[QU];
    gd_state(i, j, k, v_idx)  = riemann.state[QV];
    gd_state(i, j, k, w_idx)  = riemann.state[QW];   
    gd_state(i, j, k, QREINT) = riemann.state[QREINT];
    gd_state(i, j, k, QPRES)  = riemann.state[QPRES];
    gd_state(i, j, k, QY)     = riemann.state[QY];

    gd_flux(i, j, k, URHO)  = riemann.flux[URHO];
    gd_flux(i, j, k, u_idx) = riemann.flux[UMX];   
    gd_flux(i, j, k, v_idx) = riemann.flux[UMY];   
    gd_flux(i, j, k, w_idx) = riemann.flux[UMZ];   
    gd_flux(i, j, k, UEDEN) = riemann.flux[UEDEN];
    gd_flux(i, j, k, UY)    = riemann.flux[UY];
}
#endif
