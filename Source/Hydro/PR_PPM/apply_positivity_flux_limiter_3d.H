AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_positivity_flux_limiter_3d(
    int i, int j, int k,
    amrex::Array4<amrex::Real const> const& U_old, // Old conservative states
    amrex::Array4<amrex::Real> const& fx,         // X-fluxes (to be modified)
    amrex::Array4<amrex::Real> const& fy,         // Y-fluxes (to be modified)
    amrex::Array4<amrex::Real> const& fz,         // Z-fluxes (to be modified)
    PRParameters const& params,
    amrex::Real const dt, amrex::Real const dx, amrex::Real const dy, amrex::Real const dz)
{
    using namespace StateIdx3D;
    
    amrex::Real const rho_floor = 1e-5;
    amrex::Real const e_floor   = params.R * 10.0;

    // Predict cell-centered evolution using current second-order fluxes
    amrex::Real const rho_old = U_old(i, j, k, RHO);
    amrex::Real const E_old   = U_old(i, j, k, ENG);

    amrex::Real const dF_rho = (fx(i+1, j, k, RHO) - fx(i, j, k, RHO)) / dx +
                               (fy(i, j+1, k, RHO) - fy(i, j, k, RHO)) / dy +
                               (fz(i, j, k+1, RHO) - fz(i, j, k, RHO)) / dz;

    amrex::Real const dF_eng = (fx(i+1, j, k, ENG) - fx(i, j, k, ENG)) / dx +
                               (fy(i, j+1, k, ENG) - fy(i, j, k, ENG)) / dy +
                               (fz(i, j, k+1, ENG) - fz(i, j, k, ENG)) / dz;

    amrex::Real const rho_predicted = rho_old - dt * dF_rho;
    amrex::Real const E_predicted   = E_old   - dt * dF_eng;

    // Extrapolate projected kinetic fields to calculate predicted internal energy
    amrex::Real const mx_pred = U_old(i, j, k, MX) - dt * ((fx(i+1, j, k, MX) - fx(i, j, k, MX))/dx + (fy(i, j+1, k, MX) - fy(i, j, k, MX))/dy + (fz(i, j, k+1, MX) - fz(i, j, k, MX))/dz);
    amrex::Real const my_pred = U_old(i, j, k, MY) - dt * ((fx(i+1, j, k, MY) - fx(i, j, k, MX))/dx + (fy(i, j+1, k, MY) - fy(i, j, k, MY))/dy + (fz(i, j, k+1, MY) - fz(i, j, k, MY))/dz);
    amrex::Real const mz_pred = U_old(i, j, k, MZ) - dt * ((fx(i+1, j, k, MZ) - fx(i, j, k, MX))/dx + (fy(i, j+1, k, MZ) - fy(i, j, k, MY))/dy + (fz(i, j, k+1, MZ) - fz(i, j, k, MZ))/dz);
    
    amrex::Real const ke_predicted = 0.5 * (mx_pred*mx_pred + my_pred*my_pred + mz_pred*mz_pred) / amrex::max(rho_floor, rho_predicted);
    amrex::Real const e_predicted  = (E_predicted - ke_predicted) / amrex::max(rho_floor, rho_predicted);

    // Determine if the update violates physical positivity bounds
    if (rho_predicted < rho_floor || e_predicted < e_floor) {
        // Compute an adaptive scaling factor theta to choke the flux divergences back to safety
        amrex::Real theta = 1.0;
        
        if (rho_predicted < rho_floor) {
            theta = amrex::min(theta, (rho_old - rho_floor) / amrex::max(1e-12, dt * dF_rho));
        }
        if (e_predicted < e_floor) {
            amrex::Real const current_e = (E_old - 0.5*(U_old(i,j,k,MX)*U_old(i,j,k,MX)+U_old(i,j,k,MY)*U_old(i,j,k,MY)+U_old(i,j,k,MZ)*U_old(i,j,k,MZ))/rho_old)/rho_old;
            theta = amrex::min(theta, (current_e - e_floor) / amrex::max(1e-12, current_e - e_predicted));
        }

        // Smoothly clamp theta parameter
        theta = amrex::clamp(theta, 0.0, 1.0);

        // Scale high-order correction portions of the fluxes surrounding cell (i,j,k)
        for (int n = 0; n <= QNS; ++n) {
            fx(i,   j,   k,   n) = theta * fx(i,   j,   k,   n);
            fx(i+1, j,   k,   n) = theta * fx(i+1, j,   k,   n);
            fy(i,   j,   k,   n) = theta * fy(i,   j,   k,   n);
            fy(i,   j+1, k,   n) = theta * fy(i,   j+1, k,   n);
            fz(i,   j,   k,   n) = theta * fz(i,   j,   k,   n);
            fz(i,   j,   k+1, n) = theta * fz(i,   j,   k+1, n);
        }
    }
}

