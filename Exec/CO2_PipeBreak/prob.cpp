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

        pp.query("p0",       CAMR::h_prob_parm->p0);
        pp.query("T0",       CAMR::h_prob_parm->T0);
        pp.query("u0",       CAMR::h_prob_parm->u0);
        pp.query("p_amb",    CAMR::h_prob_parm->p_amb);
        pp.query("p_res",    CAMR::h_prob_parm->p_res);
        pp.query("T_res",    CAMR::h_prob_parm->T_res);
        pp.query("gap_yc",   CAMR::h_prob_parm->gap_yc);
        pp.query("gap_half", CAMR::h_prob_parm->gap_half);
        pp.query("ramp_time",    CAMR::h_prob_parm->ramp_time);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         CAMR::h_prob_parm, CAMR::h_prob_parm+1,
                         CAMR::d_prob_parm);
    }
}
