#include <cstdio>
#include <cmath>
#include "gerg_co2_guard.H"
static double Tbis(double rho,double e,gerg::Phase ph){double lo=gerg::T_trip,hi=gerg::T_max;
 for(int i=0;i<64;++i){double m=.5*(lo+hi); if(gerg::state_TR_phase(m,rho,ph).e>e) hi=m; else lo=m;} return .5*(lo+hi);}
int main(int argc,char**argv){
 double T=(argc>1)?atof(argv[1]):230.0;
 printf("# VAPOR branch at fixed T=%.2f, edge rhoV_edge=%.4f, rhoV_sat=%.3f\n",T,gerg::rhoV_edge(T),gerg::rhoV_sat(T));
 printf("%10s %14s %10s %14s %10s\n","rho","P[bar]","c[m/s]","dPdrho_T","in_ext");
 double re=gerg::rhoV_edge(T);
 for(double f=0.90;f<=1.10001;f+=0.01){
   double rho=re*f; gerg::State s=gerg::state_TR_phase(T,rho,gerg::Phase::Vapor);
   printf("%10.4f %14.6f %10.2f %14.1f %10d\n",rho,s.P/1e5,s.c,s.dPdrho_T,(rho>re));
 }
 printf("\n# same, but through the (rho,e) inversion actually used by the solver (fixed e)\n");
 double e=gerg::state_TR_phase(T,re,gerg::Phase::Vapor).e;
 printf("# e = %.6e  (energy of the edge state)\n",e);
 printf("%10s %10s %14s %10s %14s\n","rho","T_inv","P[bar]","c[m/s]","dPdrho_T");
 for(double f=0.90;f<=1.10001;f+=0.01){
   double rho=re*f; double Ti=Tbis(rho,e,gerg::Phase::Vapor);
   gerg::State s=gerg::state_TR_phase(Ti,rho,gerg::Phase::Vapor);
   printf("%10.4f %10.3f %14.6f %10.2f %14.1f\n",rho,Ti,s.P/1e5,s.c,s.dPdrho_T);
 }
 return 0;}
