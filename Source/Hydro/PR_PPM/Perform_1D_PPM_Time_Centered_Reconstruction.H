#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <cmath>
#include <algorithm>

/**
 * Helper function: Lightweight van Leer slope limiter to ensure 
 * monotonicity of cell-centered slopes before characteristic tracing.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_van_leer_slope(amrex::Real q_l, amrex::Real q_c, amrex::Real q_r)
{
    amrex::Real dl = q_c - q_l;
    amrex::Real dr = q_r - q_c;
    if (dl * dr <= 0.0) return 0.0;
    return (2.0 * dl * dr) / (dl + dr);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Perform_1D_PPM_Time_Centered_Reconstruction(
    int i, int j, int k,
    int di, int dj, int dk,                          // Direction strides (e.g., 1,0,0 for X)
    int u_idx,                                       // Normal velocity component index
    amrex::Array4<const amrex::Real> const& prim,
    amrex::Real& state_L, amrex::Real& state_R,      // Outputs: Time-centered edges for "comp"
    int comp, amrex::Real dt, amrex::Real dx)
{
#if 0
    using namespace amrex;

    // 1. Identify multi-dimensional neighbor indices along the tracking axis
    int i_m2 = i - 2*di; int j_m2 = j - 2*dj; int k_m2 = k - 2*dk;
    int i_m1 = i - 1*di; int j_m1 = j - 1*dj; int k_m1 = k - 1*dk;
    int i_p1 = i + 1*di; int j_p1 = j + 1*dj; int k_p1 = k + 1*dk;
    int i_p2 = i + 2*di; int j_p2 = j + 2*dj; int k_p2 = k + 2*dk;

    // 2. Fetch baseline cell-centered primitives and cached thermodynamic variables
    Real rho  = prim(i, j, k, PrimComp::RHO);
    Real u_c  = prim(i, j, k, u_idx);
    Real P    = prim(i, j, k, PrimComp::PRES);
    Real cs   = prim(i, j, k, PrimComp::CS);
    Real GAME = prim(i, j, k, PrimComp::GAME);
    Real q_c  = prim(i, j, k, comp);

    // 3. Compute 4th-order spatial edge states inside cell (i)
    Real q_L_spatial = (7.0/12.0) * (q_c + prim(i_m1,j_m1,k_m1, comp)) - 
                       (1.0/12.0) * (prim(i_p1,j_p1,k_p1, comp) + prim(i_m2,j_m2,k_m2, comp));
              
    Real q_R_spatial = (7.0/12.0) * (q_c + prim(i_p1,j_p1,k_p1, comp)) - 
                       (1.0/12.0) * (prim(i_m1,j_m1,k_m1, comp) + prim(i_p2,j_p2,k_p2, comp));

    // Monotonize baseline spatial edge values to local cell extrema
    if ((q_R_spatial - q_c)*(q_c - q_L_spatial) <= 0.0) {
        q_L_spatial = q_c; q_R_spatial = q_c;
    } else if (std::abs(q_R_spatial - q_c) > 2.0 * std::abs(q_L_spatial - q_c)) {
        q_R_spatial = 3.0 * q_c - 2.0 * q_L_spatial;
    } else if (std::abs(q_L_spatial - q_c) > 2.0 * std::abs(q_R_spatial - q_c)) {
        q_L_spatial = 3.0 * q_c - 2.0 * q_R_spatial;
    }

    // 4. Compute limited cell slopes for the characteristic expansion tracking
    Real dP  = compute_van_leer_slope(prim(i_m1,j_m1,k_m1, PrimComp::PRES), P, prim(i_p1,j_p1,k_p1, PrimComp::PRES));
    Real du  = compute_van_leer_slope(prim(i_m1,j_m1,k_m1, u_idx),          u_c, prim(i_p1,j_p1,k_p1, u_idx));
    Real dq  = compute_van_leer_slope(prim(i_m1,j_m1,k_m1, comp),          q_c, prim(i_p1,j_p1,k_p1, comp));

    // Dimensionless wave parameters (Courant numbers for each characteristic wave path)
    Real sigma_minus = (u_c - cs) * dt / dx;
    Real sigma_0     = u_c * dt / dx;
    Real sigma_plus  = (u_c + cs) * dt / dx;

    // =========================================================================
    // 5. CHARACTERISTIC WAVE CLIPPING LOGIC FOR TIME-CENTERING
    // =========================================================================
    Real sum_L = 0.0;
    Real sum_R = 0.0;

    if (comp == PrimComp::PRES) {
        Real alpha_minus = 0.5 * (dP / (rho * cs) - du);
        Real alpha_plus  = 0.5 * (dP / (rho * cs) + du);

        if (sigma_minus > 0.0) sum_L += alpha_minus * (1.0 - sigma_minus) * (rho * cs);
        if (sigma_plus  > 0.0) sum_L += alpha_plus  * (1.0 - sigma_plus)  * (rho * cs);

        if (sigma_minus < 0.0) sum_R += alpha_minus * (-1.0 - sigma_minus) * (rho * cs);
        if (sigma_plus  < 0.0) sum_R += alpha_plus  * (-1.0 - sigma_plus)  * (rho * cs);
    } 
    else if (comp == u_idx) {
        Real alpha_minus = 0.5 * (dP / (rho * cs) - du);
        Real alpha_plus  = 0.5 * (dP / (rho * cs) + du);

        if (sigma_minus > 0.0) sum_L -= alpha_minus * (1.0 - sigma_minus);
        if (sigma_plus  > 0.0) sum_L += alpha_plus  * (1.0 - sigma_plus);

        if (sigma_minus < 0.0) sum_R -= alpha_minus * (-1.0 - sigma_minus);
        if (sigma_plus  < 0.0) sum_R += alpha_plus  * (-1.0 - sigma_plus);
    } 
    else {
        // --- ADVECTION SCALAR ZONE (Density, Tangential Velocities, GAME, and cs) ---
        // These fields are decoupled from acoustic waves and travel strictly at u_c
        
        // Left Edge Time-Centering: Wave entering from the left means u_c > 0
        if (sigma_0 > 0.0) {
            sum_L += dq * (1.0 - sigma_0);
        }
        
        // Right Edge Time-Centering: Wave entering from the right means u_c < 0
        if (sigma_0 < 0.0) {
            sum_R += dq * (-1.0 - sigma_0);
        }

        // Density has an added correction term from acoustic pressure wave expansion
        if (comp == PrimComp::RHO) {
            Real alpha_minus = 0.5 * (dP / (rho * cs) - du);
            Real alpha_plus  = 0.5 * (dP / (rho * cs) + du);

            if (sigma_minus > 0.0) sum_L += alpha_minus * (1.0 - sigma_minus) * (1.0 / (cs * cs));
            if (sigma_plus  > 0.0) sum_L += alpha_plus  * (1.0 - sigma_plus)  * (1.0 / (cs * cs));
            
            if (sigma_minus < 0.0) sum_R += alpha_minus * (-1.0 - sigma_minus) * (1.0 / (cs * cs));
            if (sigma_plus  < 0.0) sum_R += alpha_plus  * (-1.0 - sigma_plus)  * (1.0 / (cs * cs));
        }
    }

    // =========================================================================
    // 6. Final Assembly & Extraction of Time-Centered Interface Values
    // =========================================================================
    state_L = q_L_spatial + 0.5 * sum_L;
    state_R = q_R_spatial + 0.5 * sum_R;

    // Enforce physical bounds for thermodynamic interpolations
    if (comp == PrimComp::RHO || comp == PrimComp::PRES || comp == PrimComp::CS) {
        if (state_L <= 1e-10) state_L = 1e-10;
        if (state_R <= 1e-10) state_R = 1e-10;
    }
    if (comp == PrimComp::GAME) {
        if (state_L < 1.001) state_L = 1.001;
        if (state_R < 1.001) state_R = 1.001;
    }
#endif
}

