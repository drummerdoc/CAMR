#include "AMReX_PROB_AMR_F.H"
#include "AMReX_ParmParse.H"
#include "CAMR.H"
#include "prob.H"
#include "EOS.H"

extern "C" {
    void amrex_probinit(const int* /*init*/,
                        const int* /*name*/,
                        const int* /*namelen*/,
                        const amrex_real* /*problo*/,
                        const amrex_real* /*probhi*/)
    {
        amrex::ParmParse pp("prob");

        pp.query("p0",       CAMR::h_prob_parm->p0);
        pp.query("T0",       CAMR::h_prob_parm->T0);
        pp.query("u0",       CAMR::h_prob_parm->u0);
        pp.query("p_amb",    CAMR::h_prob_parm->p_amb);
        pp.query("p_res",    CAMR::h_prob_parm->p_res);
        pp.query("T_res",    CAMR::h_prob_parm->T_res);
        pp.query("gap_yc",    CAMR::h_prob_parm->gap_yc);
        pp.query("gap_half",  CAMR::h_prob_parm->gap_half);
        pp.query("gap_taper", CAMR::h_prob_parm->gap_taper);
        pp.query("ramp_time", CAMR::h_prob_parm->ramp_time);
        pp.query("res_alpha1", CAMR::h_prob_parm->res_alpha1);
        pp.query("res_dome_taper", CAMR::h_prob_parm->res_dome_taper);
        pp.query("res_u",      CAMR::h_prob_parm->res_u);
        pp.query("sym_ylo",    CAMR::h_prob_parm->sym_ylo);
        pp.query("amb_vapor",  CAMR::h_prob_parm->amb_vapor);

        // Two-phase saturated reservoir (rung 4): precompute the saturated
        // liquid/vapor per-phase (rho,e) at T_res once, so the per-cell
        // bcnormal fill is cheap.
        if (CAMR::h_prob_parm->res_alpha1 < amrex::Real(1.0)) {
            amrex::Real rhoL, eL, rhoV, eV;
            EOS::co2_sat_LV(CAMR::h_prob_parm->T_res, rhoL, eL, rhoV, eV);
            CAMR::h_prob_parm->res_rhoL = rhoL;
            CAMR::h_prob_parm->res_eL   = eL;
            CAMR::h_prob_parm->res_rhoV = rhoV;
            CAMR::h_prob_parm->res_eV   = eV;
            amrex::Print() << "  CO2_PipeBreak: saturated two-phase reservoir @ T_res="
                           << CAMR::h_prob_parm->T_res << " K, alpha1_liq="
                           << CAMR::h_prob_parm->res_alpha1
                           << "  rhoL=" << rhoL << " rhoV=" << rhoV
                           << " eL=" << eL << " eV=" << eV << std::endl;
        }

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         CAMR::h_prob_parm, CAMR::h_prob_parm+1,
                         CAMR::d_prob_parm);
    }
}
