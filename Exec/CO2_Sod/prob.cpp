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

        pp.query("type",  CAMR::h_prob_parm->type);
        pp.query("p_l",   CAMR::h_prob_parm->p_l);
        pp.query("p_r",   CAMR::h_prob_parm->p_r);
        pp.query("T_l",   CAMR::h_prob_parm->T_l);
        pp.query("T_r",   CAMR::h_prob_parm->T_r);
        pp.query("u_l",   CAMR::h_prob_parm->u_l);
        pp.query("u_r",   CAMR::h_prob_parm->u_r);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         CAMR::h_prob_parm, CAMR::h_prob_parm+1,
                         CAMR::d_prob_parm);
    }
}
