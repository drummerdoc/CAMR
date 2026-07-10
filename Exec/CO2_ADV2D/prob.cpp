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

        pp.query("P0", CAMR::h_prob_parm->P0);
        pp.query("T1", CAMR::h_prob_parm->T1);
        pp.query("T2", CAMR::h_prob_parm->T2);

        pp.query("u0", CAMR::h_prob_parm->u0);
        pp.query("v0", CAMR::h_prob_parm->v0);

        pp.query("alpha_mean", CAMR::h_prob_parm->alpha_mean);
        pp.query("alpha_amp",  CAMR::h_prob_parm->alpha_amp);
        pp.query("k_x",        CAMR::h_prob_parm->k_x);
        pp.query("k_y",        CAMR::h_prob_parm->k_y);
        pp.query("alpha_pow",  CAMR::h_prob_parm->alpha_pow);

        CAMR::h_prob_parm->p_amb = CAMR::h_prob_parm->P0;

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         CAMR::h_prob_parm, CAMR::h_prob_parm+1,
                         CAMR::d_prob_parm);
    }
}
