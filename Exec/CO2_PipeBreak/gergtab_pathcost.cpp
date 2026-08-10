// Measure which EOS path the per-phase branch inversion takes for REAL cell
// states, and what the table-seeded polish would cost/achieve.
// stdin rows: rho e a1 rho1 e1 rho2 e2   (same format as gergstats)
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "gerg_co2_guard.H"
#include "gergtab_bicubic.H"

static long nev;   // guarded-eval counter
static gerg::State ev(double T, double rho, gerg::Phase ph)
{ ++nev; return gerg::state_TR_phase(T, rho, ph); }

static double bisect(double rho, double e, gerg::Phase ph)
{ double lo=gerg::T_trip, hi=gerg::T_max;
  for (int i=0;i<64;++i){double m=.5*(lo+hi); if (ev(m,rho,ph).e>e) hi=m; else lo=m;}
  return .5*(lo+hi); }

// g_polish_phase with a configurable step count / clamp
static bool polish(double rho, double e, gerg::Phase ph, double T,
                   int nstep, double clamp, double& Tout)
{
    for (int it=0; it<nstep; ++it) {
        const gerg::State s0 = ev(T, rho, ph);
        const gerg::State s1 = ev(T+0.1, rho, ph);
        const double slope = std::fmax((s1.e-s0.e)/0.1, 100.0);
        double dT = (e - s0.e)/slope;
        if (dT >  clamp) dT =  clamp;
        if (dT < -clamp) dT = -clamp;
        T += dT;
        if (T < gerg::T_trip) T = gerg::T_trip;
        if (T > gerg::T_max)  T = gerg::T_max;
    }
    Tout = T;
    return std::fabs(ev(T,rho,ph).e - e) <= 1.0e-6*std::fmax(std::fabs(e), 1.0e4);
}

struct Acc { long n=0, ok=0, dome=0, cur_bisect=0, acc2=0, acc4=0, acc6=0;
             double ev_cur=0, ev_new2=0, ev_new4=0, ev_new6=0, maxdT=0; };

static void one(double rho, double e, gerg::Phase ph, Acc& a)
{
    ++a.n;
    double Tt; const bool ok = (ph==gerg::Phase::Liquid)
        ? gergtab::T_liq(rho,e,Tt) : gergtab::T_vap(rho,e,Tt);
    if (ok) ++a.ok;
    const bool indome = ok && (Tt < gerg::SAT_T_HI &&
                       rho < gerg::rhoL_sat(Tt) && rho > gerg::rhoV_sat(Tt));
    if (indome) ++a.dome;

    nev=0; const double Tref = bisect(rho,e,ph); const double ev_bis = nev;

    // CURRENT code path
    if (ok && !indome) { nev=0; double To; polish(rho,e,ph,Tt,2,2.0,To);
                         a.ev_cur += nev; }
    else { a.ev_cur += ev_bis; ++a.cur_bisect; }

    // PROPOSED: always seed from the table when the lookup succeeds
    if (ok) {
        if (std::fabs(Tt-Tref) > a.maxdT) a.maxdT = std::fabs(Tt-Tref);
        for (int k=0;k<3;++k) {
            const int ns = (k==0)?2:((k==1)?4:6);
            nev=0; double To; const bool good = polish(rho,e,ph,Tt,ns,2.0,To);
            double cost = nev + (good?0.0:ev_bis);
            if (k==0){ a.ev_new2+=cost; if(good) ++a.acc2; }
            if (k==1){ a.ev_new4+=cost; if(good) ++a.acc4; }
            if (k==2){ a.ev_new6+=cost; if(good) ++a.acc6; }
        }
    } else { a.ev_new2+=ev_bis; a.ev_new4+=ev_bis; a.ev_new6+=ev_bis; }
}

int main()
{
    double rho,e,a1,r1,e1,r2,e2;
    Acc L, V;
    while (scanf("%lf %lf %lf %lf %lf %lf %lf",&rho,&e,&a1,&r1,&e1,&r2,&e2)==7) {
        one(r1,e1,gerg::Phase::Liquid,L);
        one(r2,e2,gerg::Phase::Vapor ,V);
    }
    const char* nm[2] = {"LIQUID","VAPOR "};
    Acc* aa[2] = {&L,&V};
    for (int k=0;k<2;++k) {
        Acc& a = *aa[k]; if (!a.n) continue;
        printf("%s n=%ld  table_ok=%.1f%%  in-dome(rejected)=%.1f%%  "
               "CURRENT bisect-path=%.1f%%\n", nm[k], a.n,
               100.0*a.ok/a.n, 100.0*a.dome/a.n, 100.0*a.cur_bisect/a.n);
        printf("        mean guarded evals/call: CURRENT %6.1f | polish2 %5.1f "
               "(accept %.1f%%) | polish4 %5.1f (accept %.1f%%) | polish6 %5.1f (accept %.1f%%)\n",
               a.ev_cur/a.n, a.ev_new2/a.n, 100.0*a.acc2/std::fmax(a.ok,1.0),
               a.ev_new4/a.n, 100.0*a.acc4/std::fmax(a.ok,1.0),
               a.ev_new6/a.n, 100.0*a.acc6/std::fmax(a.ok,1.0));
        printf("        speedup vs CURRENT: polish2 %.1fx  polish4 %.1fx  polish6 %.1fx"
               "   max |T_seed - T_exact| = %.3f K\n",
               a.ev_cur/std::fmax(a.ev_new2,1.0), a.ev_cur/std::fmax(a.ev_new4,1.0),
               a.ev_cur/std::fmax(a.ev_new6,1.0), a.maxdT);
    }
    return 0;
}
