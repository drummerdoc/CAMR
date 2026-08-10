#include "gerg_co2_guard.H"
#include <cstdio>
#include <cmath>
int main(){
    double r1,T1,r2,T2;
    while(scanf("%lf %lf %lf %lf",&r1,&T1,&r2,&T2)==4){
        gerg::State a=gerg::state_TR_phase(T1,r1,gerg::Phase::Liquid);
        gerg::State b=gerg::state_TR_phase(T2,r2,gerg::Phase::Vapor);
        printf("%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %d\n",
               a.g,b.g,a.h,b.h,a.s,b.s,a.e,b.e,
               (std::isfinite(a.g)&&std::isfinite(b.g))?1:0);
    }
    return 0;
}
