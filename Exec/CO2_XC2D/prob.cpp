#include "AMReX_PROB_AMR_F.H"
#include "AMReX_ParmParse.H"
#include "CAMR.H"
#include "prob.H"

extern "C" {
    void amrex_probinit(const int* /*init*/,
                        const int* /*name*/,
                        const int* /*namelen*/,
                        const amrex_real* /*problo*/,
                        const amrex_real* /*probhi*/)
    {
        amrex::ParmParse pp("prob");

        pp.query("p_L",      CAMR::h_prob_parm->p_L);
        pp.query("T_L",      CAMR::h_prob_parm->T_L);
        pp.query("u_L",      CAMR::h_prob_parm->u_L);

        pp.query("p_R",      CAMR::h_prob_parm->p_R);
        pp.query("T_R",      CAMR::h_prob_parm->T_R);
        pp.query("u_R",      CAMR::h_prob_parm->u_R);

        pp.query("x_diaph",  CAMR::h_prob_parm->x_diaph);

        CAMR::h_prob_parm->p_amb_lo = CAMR::h_prob_parm->p_L;
        CAMR::h_prob_parm->p_amb_hi = CAMR::h_prob_parm->p_R;
        pp.query("p_amb_lo", CAMR::h_prob_parm->p_amb_lo);
        pp.query("p_amb_hi", CAMR::h_prob_parm->p_amb_hi);

        CAMR::h_prob_parm->p_amb = CAMR::h_prob_parm->p_amb_hi;
        pp.query("p_amb", CAMR::h_prob_parm->p_amb);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         CAMR::h_prob_parm, CAMR::h_prob_parm+1,
                         CAMR::d_prob_parm);
    }
}
