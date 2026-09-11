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
        ProbParmDevice& p = *CAMR::h_prob_parm;
        pp.query("p0",    p.p0);
        pp.query("T0",    p.T0);
        pp.query("u0",    p.u0);
        const bool has_pamb = pp.query("p_amb", p.p_amb);
        pp.query("ic",    p.ic);
        pp.query("amp",   p.amp);
        pp.query("xc",    p.xc);
        pp.query("yc",    p.yc);
        pp.query("w",     p.w);
        pp.query("two_phase", p.two_phase);
        pp.query("a1_0",  p.a1_0);
        pp.query("alpha_trace", p.alpha_trace);
        if (p.two_phase != 0) {
            EOS::co2_sat_LV(p.T0, p.rhoL, p.eL, p.rhoV, p.eV);
            // The mixture's own pressure: the PR coexistence pressure of the
            // saturated liquid state (EOS::Psat is the reference-fluid
            // correlation and differs from it, 41.95 vs 41.59 bar at 280 K).
            amrex::Real Y[NUM_SPECIES]; Y[0] = 1.0; for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = 0.0;
            EOS::REY2P_phase(p.rhoL, p.eL, Y, hem::Phase3::Liquid, p.p0);
            if (!has_pamb) p.p_amb = p.p0;   // far field on the dome unless the deck says otherwise
            amrex::Print() << "  CO2_NSCBC: saturated mixture at T0=" << p.T0
                           << " K, P_LV=" << p.p0/1.0e5 << " bar (EOS::Psat " << EOS::Psat(p.T0)/1.0e5 << ")" << " bar, alpha_1=" << p.a1_0
                           << " (rhoL=" << p.rhoL << ", rhoV=" << p.rhoV << ")\n";
        }
        if (p.ic < 0 || p.ic > 4) amrex::Abort("prob.ic: 0..4");
        if (p.ic == 4 && p.two_phase == 0) amrex::Abort("prob.ic=4 needs prob.two_phase=1");
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         CAMR::h_prob_parm, CAMR::h_prob_parm+1,
                         CAMR::d_prob_parm);
    }
}
