#include "hem_pr_state.H"
#include <cstdio>
#include <cmath>
int main(){
    const hem::PRFluid CO2 = hem::PRFluid::make(304.13,7.3773e6,0.22394,0.04401,
        2.35681300e+00,8.98412990e-03,-7.12206320e-06,2.45730071e-09,-1.42885166e-13,0.0);
    double r1,T1,r2,T2;
    while(scanf("%lf %lf %lf %lf",&r1,&T1,&r2,&T2)==4){
        hem::State a=hem::state_from_T_v(CO2,T1,hem::v_PR_from_rho(CO2,r1));
        hem::State b=hem::state_from_T_v(CO2,T2,hem::v_PR_from_rho(CO2,r2));
        double g1=a.h-T1*a.s, g2=b.h-T2*b.s;
        printf("%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %d\n",
               g1,g2,a.h,b.h,a.s,b.s,a.e,b.e,
               (std::isfinite(g1)&&std::isfinite(g2))?1:0);
    }
    return 0;
}
