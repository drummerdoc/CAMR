// GERG Gibbs-driving-force probe. Same I/O contract as probe_pr.
// (rho,e) -> T by bisection on the phase-locked branch, then read State.g.
#include "gerg_co2_guard.H"
#include <cstdio>
#include <cmath>
static double T_from_e(double rho,double e,gerg::Phase ph)
{
    double lo=180.0, hi=420.0;
    for (int i=0;i<200;++i) {
        double m=0.5*(lo+hi);
        double em=gerg::state_TR_phase(m,rho,ph).e;
        if (!std::isfinite(em)) return -1.0;
        if (em>e) hi=m; else lo=m;
    }
    return 0.5*(lo+hi);
}
int main()
{
    double r1,e1,r2,e2;
    while (scanf("%lf %lf %lf %lf",&r1,&e1,&r2,&e2)==4) {
        double T1=T_from_e(r1,e1,gerg::Phase::Liquid);
        double T2=T_from_e(r2,e2,gerg::Phase::Vapor);
        gerg::State s1=gerg::state_TR_phase(T1,r1,gerg::Phase::Liquid);
        gerg::State s2=gerg::state_TR_phase(T2,r2,gerg::Phase::Vapor);
        int ok = T1>0&&T2>0&&std::isfinite(s1.g)&&std::isfinite(s2.g);
        printf("%.8e %.8e %.8e %.8e %.8e %d\n",T1,T2,s1.g,s2.g,s1.g-s2.g,ok);
    }
    return 0;
}
