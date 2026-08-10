#include "gerg_co2_guard.H"
#include <cstdio>
int main(){
    printf("# T Psat[bar] rhoL rhoV gL gV gL-gV\n");
    for (double T=220.0; T<=300.0; T+=5.0) {
        double rL=gerg::rhoL_sat(T), rV=gerg::rhoV_sat(T);
        gerg::State L=gerg::state_TR_phase(T,rL,gerg::Phase::Liquid);
        gerg::State V=gerg::state_TR_phase(T,rV,gerg::Phase::Vapor);
        printf("%6.1f %9.4f %9.3f %9.4f %14.6e %14.6e %13.5e\n",
               T,gerg::Psat_spl(T)/1e5,rL,rV,L.g,V.g,L.g-V.g);
    }
    return 0;
}
