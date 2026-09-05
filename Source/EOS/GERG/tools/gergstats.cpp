// Bulk classification of PS per-phase states against the guarded GERG surface.
// stdin rows: rho e a1 rho1 e1 rho2 e2
// stdout: one row per cell: P1 P2 Pmix T1 T2 c1 c2 cmix flag_ext1 flag_ext2 flag_Tclamp1 flag_Tclamp2
#include <cstdio>
#include <cmath>
#include "../gerg_co2_guard.H"

static double Tbis_phase(double rho, double e, gerg::Phase ph)
{ double lo=gerg::T_trip, hi=gerg::T_max;
  for (int i=0;i<64;++i){double m=0.5*(lo+hi); if (gerg::state_TR_phase(m,rho,ph).e>e) hi=m; else lo=m;}
  return 0.5*(lo+hi); }

int main()
{
    double rho,e,a1,r1,e1,r2,e2;
    while (scanf("%lf %lf %lf %lf %lf %lf %lf",&rho,&e,&a1,&r1,&e1,&r2,&e2)==7) {
        const double T1=Tbis_phase(r1,e1,gerg::Phase::Liquid);
        const gerg::State s1=gerg::state_TR_phase(T1,r1,gerg::Phase::Liquid);
        const double T2=Tbis_phase(r2,e2,gerg::Phase::Vapor);
        const gerg::State s2=gerg::state_TR_phase(T2,r2,gerg::Phase::Vapor);
        const double Pmix=a1*s1.P+(1.0-a1)*s2.P;
        const double inv = a1/(r1*s1.c*s1.c) + (1.0-a1)/(r2*s2.c*s2.c);
        const double cmix = 1.0/std::sqrt(rho*inv);
        const int ext1 = (r1 < gerg::rhoL_edge(T1)) ? 1 : 0;
        const int ext2 = (r2 > gerg::rhoV_edge(T2)) ? 1 : 0;
        const int cl1 = (T1 <= gerg::T_trip*(1.0+1e-9)) ? 1 : 0;
        const int cl2 = (T2 <= gerg::T_trip*(1.0+1e-9)) ? 1 : 0;
        const int hi1 = (T1 >= gerg::T_max*(1.0-1e-9)) ? 1 : 0;
        const int hi2 = (T2 >= gerg::T_max*(1.0-1e-9)) ? 1 : 0;
        // g_T_from_e_phase (GERGTab) REJECTS table service for IN-DOME branch
        // states and falls through to the 64-iteration g_T_bisect_phase.
        const int dome1 = (T1 < gerg::SAT_T_HI &&
                           r1 < gerg::rhoL_sat(T1) && r1 > gerg::rhoV_sat(T1)) ? 1 : 0;
        const int dome2 = (T2 < gerg::SAT_T_HI &&
                           r2 < gerg::rhoL_sat(T2) && r2 > gerg::rhoV_sat(T2)) ? 1 : 0;
        printf("%.6e %.6e %.6e %.4f %.4f %.3f %.3f %.3f %d %d %d %d %d %d %d %d\n",
               s1.P, s2.P, Pmix, T1, T2, s1.c, s2.c, cmix, ext1, ext2, cl1, cl2, hi1, hi2,
               dome1, dome2);
    }
    return 0;
}
