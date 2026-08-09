// =====================================================================
//  gen_gergtab.cpp -- GERGTab table generator (self-contained: computes
//  AND emits gergtab_table_params.H + gergtab_table_data.cpp).
//
//  TABULATION TARGET (deliberate, differs from PRTab's "physical
//  surface with holes"): the table stores the GERG backend's OWN
//  fixed-count-bisection T(rho,e) inversion on each guarded surface
//  (auto + branchL + branchV) -- i.e., exactly the function the
//  analytic fallback computes (Source/EOS/GERG/EOS.H g_T_from_e_*).
//  That function is TOTAL, so there are no holes to fill, and
//  out-of-domain runtime queries fall back to the identical function
//  (no accuracy cliff).  In far-extension zones where e(T) is
//  non-monotone the bisection root can jump; the table then smooths a
//  jump the fallback also has -- localized, trace-phase territory.
//
//  Build/run (host g++, no AMReX):
//      g++ -O2 -std=c++17 gen_gergtab.cpp -o gen_gergtab && ./gen_gergtab [N]
//  Emits into the parent directory (Source/EOS/GERGTab/).
//  Also prints an off-grid table-vs-bisection accuracy sweep.
// =====================================================================
#include "../../GERG/gerg_co2_guard.H"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

static double T_from_e(double rho, double e, int mode)   // mode 0=auto,1=L,2=V
{
    double lo = gerg::T_trip, hi = gerg::T_max;
    for (int it = 0; it < 64; ++it) {                      // MUST match GERG/EOS.H
        const double m = 0.5*(lo + hi);
        double em;
        if (mode == 0) em = gerg::state_TR_auto(m, rho).e;
        else em = gerg::state_TR_phase(m, rho,
                    (mode==1) ? gerg::Phase::Liquid : gerg::Phase::Vapor).e;
        if (em > e) hi = m; else lo = m;
    }
    return 0.5*(lo + hi);
}

int main(int argc, char** argv)
{
    const int N = (argc > 1) ? std::atoi(argv[1]) : 256;
    const double LR0 = -3.0, LR1 = std::log10(1500.0);
    // e range: span of the AUTO (equilibrium, bounded) surface + margin
    double e_lo = 1e300, e_hi = -1e300;
    for (int i = 0; i <= 48; ++i) {
        const double r = std::pow(10.0, LR0 + (LR1-LR0)*i/48.0);
        for (double T : {gerg::T_trip, gerg::T_max}) {
            const double e = gerg::state_TR_auto(T, r).e;
            if (e < e_lo) e_lo = e; if (e > e_hi) e_hi = e;
        }
    }
    const double mar = 0.02*(e_hi - e_lo);
    e_lo -= mar; e_hi += mar;
    std::printf("grid %dx%d  logr [%g, %g]  e [%g, %g]\n", N, N, LR0, LR1, e_lo, e_hi);

    std::vector<double> TA(N*N), TL(N*N), TV(N*N);
    for (int i = 0; i < N; ++i) {
        const double r = std::pow(10.0, LR0 + (LR1-LR0)*i/(N-1));
        for (int j = 0; j < N; ++j) {
            const double e = e_lo + (e_hi-e_lo)*j/(N-1);
            TA[i*N+j] = T_from_e(r, e, 0);
            TL[i*N+j] = T_from_e(r, e, 1);
            TV[i*N+j] = T_from_e(r, e, 2);
        }
    }

    // ---- VALIDITY MASKS: the guarded surfaces are only C1 at the
    //  extension seams / dome boundary; bicubic degrades to ~h^2 there
    //  (measured: relT 1e-3 -> relP 0.43 in the stiff liquid).  Flag
    //  curvature spikes (2nd differences of T well above the smooth-
    //  region level), DILATE by the 4x4 stencil footprint, and let the
    //  runtime fall back to the identical analytic bisection there.
    auto mask_of = [&](const std::vector<double>& A) {
        std::vector<unsigned char> bad(N*N, 0), dil(N*N, 0);
        // curvature threshold: 8x the 90th percentile of |D2| (robust)
        std::vector<double> d2s; d2s.reserve(2*N*N);
        auto D2e = [&](int i,int j){ return std::fabs(A[i*N+j-1]-2*A[i*N+j]+A[i*N+j+1]); };
        auto D2r = [&](int i,int j){ return std::fabs(A[(i-1)*N+j]-2*A[i*N+j]+A[(i+1)*N+j]); };
        for (int i = 1; i < N-1; ++i) for (int j = 1; j < N-1; ++j) {
            d2s.push_back(D2e(i,j)); d2s.push_back(D2r(i,j));
        }
        std::vector<double> tmp = d2s;
        std::nth_element(tmp.begin(), tmp.begin()+size_t(0.90*tmp.size()), tmp.end());
        const double thresh = 8.0*std::max(tmp[size_t(0.90*tmp.size())], 1e-6);
        for (int i = 1; i < N-1; ++i) for (int j = 1; j < N-1; ++j)
            if (D2e(i,j) > thresh || D2r(i,j) > thresh) bad[i*N+j] = 1;
        // dilate by 2 in each direction (stencil reach) + domain edges
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
            bool b = (i < 2 || i > N-3 || j < 2 || j > N-3);
            for (int di = -2; di <= 2 && !b; ++di)
                for (int dj = -2; dj <= 2 && !b; ++dj) {
                    const int ii = i+di, jj = j+dj;
                    if (ii>=0 && ii<N && jj>=0 && jj<N && bad[ii*N+jj]) b = true;
                }
            dil[i*N+j] = b ? 0 : 1;   // 1 = table USABLE at this stencil base
        }
        long ok = 0; for (auto v : dil) ok += v;
        std::printf("  mask coverage %.1f%% (thresh %.3g K)\n", 100.0*ok/(N*N), thresh);
        return dil;
    };
    std::printf("validity masks:\n");
    const std::vector<unsigned char> OKA = mask_of(TA), OKL = mask_of(TL), OKV = mask_of(TV);

    // ---- emit params ----
    std::FILE* f = std::fopen("../gergtab_table_params.H", "w");
    std::fprintf(f,
        "#ifndef GERGTAB_TABLE_PARAMS_H\n#define GERGTAB_TABLE_PARAMS_H\n"
        "// GENERATED by tools/gen_gergtab.cpp -- do not edit.\n"
        "namespace gergtab {\n"
        "constexpr int    NLR = %d;\nconstexpr int    NE  = %d;\n"
        "constexpr double LR0 = %.17g;\nconstexpr double LR1 = %.17g;\n"
        "constexpr double E0  = %.17g;\nconstexpr double E1  = %.17g;\n"
        "extern const double GT_TA[NLR*NE];\n"
        "extern const double GT_TL[NLR*NE];\n"
        "extern const double GT_TV[NLR*NE];\n"
        "extern const unsigned char GT_OKA[NLR*NE];\n"
        "extern const unsigned char GT_OKL[NLR*NE];\n"
        "extern const unsigned char GT_OKV[NLR*NE];\n"
        "}\n#endif\n", N, N, LR0, LR1, e_lo, e_hi);
    std::fclose(f);

    // ---- emit data ----
    f = std::fopen("../gergtab_table_data.cpp", "w");
    std::fprintf(f, "// GENERATED by tools/gen_gergtab.cpp -- do not edit.\n"
                    "#include \"gergtab_table_params.H\"\nnamespace gergtab {\n");
    auto emit = [&](const char* nm, const std::vector<double>& a) {
        std::fprintf(f, "extern const double %s[NLR*NE] = {", nm);
        for (size_t k = 0; k < a.size(); ++k)
            std::fprintf(f, "%s%.17g", k ? "," : "", a[k]);
        std::fprintf(f, "};\n");
    };
    emit("GT_TA", TA); emit("GT_TL", TL); emit("GT_TV", TV);
    auto emitm = [&](const char* nm, const std::vector<unsigned char>& a) {
        std::fprintf(f, "extern const unsigned char %s[NLR*NE] = {", nm);
        for (size_t k = 0; k < a.size(); ++k)
            std::fprintf(f, "%s%d", k ? "," : "", int(a[k]));
        std::fprintf(f, "};\n");
    };
    emitm("GT_OKA", OKA); emitm("GT_OKL", OKL); emitm("GT_OKV", OKV);
    std::fprintf(f, "}\n");
    std::fclose(f);
    std::printf("wrote ../gergtab_table_params.H, ../gergtab_table_data.cpp\n");

    // ---- off-grid accuracy sweep (bicubic done at runtime; here report
    //      bilinear-midpoint bound as a cheap proxy + exact at nodes) ----
    double worst = 0, wr = 0, we = 0; int mode_w = 0;
    for (int i = 0; i < 200; ++i) {
        const double lr = LR0 + (LR1-LR0)*(i+0.37)/200.0;
        const double r = std::pow(10.0, lr);
        for (int j = 0; j < 200; ++j) {
            const double e = e_lo + (e_hi-e_lo)*(j+0.61)/200.0;
            for (int m = 0; m < 3; ++m) {
                // bilinear from the table
                const double x = (lr-LR0)/(LR1-LR0)*(N-1), y = (e-e_lo)/(e_hi-e_lo)*(N-1);
                const int i0 = int(x), j0 = int(y);
                const double tx = x-i0, ty = y-j0;
                const std::vector<double>& A = (m==0)?TA:((m==1)?TL:TV);
                const double Tb = (1-tx)*(1-ty)*A[i0*N+j0] + tx*(1-ty)*A[(i0+1)*N+j0]
                                + (1-tx)*ty*A[i0*N+j0+1]   + tx*ty*A[(i0+1)*N+j0+1];
                const double Te = T_from_e(r, e, m);
                const double rel = std::fabs(Tb-Te)/Te;
                if (rel > worst) { worst = rel; wr = r; we = e; mode_w = m; }
            }
        }
    }
    std::printf("off-grid BILINEAR proxy: worst relT %.3e at rho %.4g e %.4g (mode %d)\n"
                "(runtime bicubic is strictly better away from kinks)\n",
                worst, wr, we, mode_w);
    return 0;
}
