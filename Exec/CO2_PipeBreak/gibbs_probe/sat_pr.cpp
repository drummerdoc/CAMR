#include "hem_pr_state.H"
#include <cstdio>
#include <cmath>
int main(){
    const hem::PRFluid CO2 = hem::PRFluid::make(304.13,7.3773e6,0.22394,0.04401,
        2.35681300e+00,8.98412990e-03,-7.12206320e-06,2.45730071e-09,-1.42885166e-13,0.0);
    printf("# T Psat[bar] rhoL rhoV gL gV gL-gV\n");
    for (double T=220.0; T<=300.0; T+=5.0) {
        hem::State L=hem::satStateL(CO2,T), V=hem::satStateV(CO2,T);
        double gL=L.h-L.T*L.s, gV=V.h-V.T*V.s;
        printf("%6.1f %9.4f %9.3f %9.4f %14.6e %14.6e %13.5e\n",
               T,0.5*(L.P+V.P)/1e5,L.rho,V.rho,gL,gV,gL-gV);
    }
    return 0;
}
