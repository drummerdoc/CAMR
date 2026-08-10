// Read "rho e a1 rho1 e1 rho2 e2" rows on stdin; emit the GERG guarded
// evaluations the PS pressure path would produce.
#include <cstdio>
#include <cmath>
#include "gerg_co2_guard.H"

static double Tbis_phase(double rho, double e, gerg::Phase ph)
{ double lo=gerg::T_trip, hi=gerg::T_max;
  for (int i=0;i<64;++i){double m=0.5*(lo+hi); if (gerg::state_TR_phase(m,rho,ph).e>e) hi=m; else lo=m;}
  return 0.5*(lo+hi); }
static double Tbis_auto(double rho, double e)
{ double lo=gerg::T_trip, hi=gerg::T_max;
  for (int i=0;i<64;++i){double m=0.5*(lo+hi); if (gerg::state_TR_auto(m,rho).e>e) hi=m; else lo=m;}
  return 0.5*(lo+hi); }

int main()
{
    double rho,e,a1,r1,e1,r2,e2;
    printf("%8s %10s %7s %8s %10s %8s %10s | %10s %8s %8s | %10s %8s %8s | %10s %10s %8s %8s\n",
      "rho","e","a1","rho1","e1","rho2","e2",
      "P1[bar]","T1","c1","P2[bar]","T2","c2","Pmix[bar]","Pauto[bar]","Tauto","cauto");
    while (scanf("%lf %lf %lf %lf %lf %lf %lf",&rho,&e,&a1,&r1,&e1,&r2,&e2)==7) {
        const double T1=Tbis_phase(r1,e1,gerg::Phase::Liquid);
        const gerg::State s1=gerg::state_TR_phase(T1,r1,gerg::Phase::Liquid);
        const double T2=Tbis_phase(r2,e2,gerg::Phase::Vapor);
        const gerg::State s2=gerg::state_TR_phase(T2,r2,gerg::Phase::Vapor);
        const double Ta=Tbis_auto(rho,e);
        const gerg::State sa=gerg::state_TR_auto(Ta,rho);
        const double Pmix=a1*s1.P+(1.0-a1)*s2.P;
        printf("%8.2f %10.4g %7.5f %8.2f %10.4g %8.3f %10.4g | %10.4f %8.2f %8.1f | %10.4f %8.2f %8.1f | %10.4f %10.4f %8.2f %8.1f\n",
          rho,e,a1,r1,e1,r2,e2, s1.P/1e5,T1,s1.c, s2.P/1e5,T2,s2.c, Pmix/1e5, sa.P/1e5,Ta,sa.c);
    }
    return 0;
}
