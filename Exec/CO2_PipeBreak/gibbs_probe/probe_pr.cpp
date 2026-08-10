// PR Gibbs-driving-force probe. stdin: "rho1 e1 rho2 e2" per line (SI).
// stdout: "T1 T2 g1 g2 g1-g2 ok"   g = h - T*s  [J/kg]
#include "hem_pr_state.H"
#include <cstdio>
#include <cmath>
int main()
{
    const hem::PRFluid CO2 = hem::PRFluid::make(
        304.13, 7.3773e6, 0.22394, 0.04401,
        2.35681300e+00, 8.98412990e-03, -7.12206320e-06,
        2.45730071e-09, -1.42885166e-13, 0.0);
    double r1,e1,r2,e2;
    while (scanf("%lf %lf %lf %lf",&r1,&e1,&r2,&e2)==4) {
        hem::State s1 = hem::state_from_rho_e_phase(CO2,r1,e1,hem::Phase3::Liquid);
        hem::State s2 = hem::state_from_rho_e_phase(CO2,r2,e2,hem::Phase3::Vapor);
        double g1 = s1.h - s1.T*s1.s, g2 = s2.h - s2.T*s2.s;
        int ok = std::isfinite(g1)&&std::isfinite(g2)&&s1.T>0&&s2.T>0;
        printf("%.8e %.8e %.8e %.8e %.8e %d\n",s1.T,s2.T,g1,g2,g1-g2,ok);
    }
    return 0;
}
