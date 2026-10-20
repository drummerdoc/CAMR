// =====================================================================
//  gpu_bench.cu — EOS-strategy microbenchmark on the harvested deployment
//  states: analytic PR inversion, fixed-count bisection, bicubic table,
//  and MLP inference (FP64 / FP32, with and without Newton polish), one
//  thread per state, as a CAMR ParallelFor would run them.
//
//  Reads data/deploy.csv (phase, rho, e, T, ... ok) and tiles the valid
//  rows to N states.  Prints ns/state and Gstates/s per kernel, the max
//  |T - T_ref| against the branch-locked inversion, and the effect of
//  warp divergence (harvest order vs sorted by branch and regime).
//
//  Build (needs the generated PRTab tables for the table kernels):  ./build_gpu_bench.sh
//  Run:  ./gpu_bench [N=2097152] [data/deploy.csv]
// =====================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cuda_runtime.h>

// --- CAMR's PR machinery, device-callable copies made by build_gpu_bench.sh ---
#include "hem_pr_state.H"

// --- the PRTab bicubic, with the tables passed as device pointers ----------
#include "prtab_table_params.H"
namespace prtab {   // the generated host arrays (prtab_table_data.cpp); copied to the device at start-up
  extern const double TBL_TL[NLR*NE], TBL_TV[NLR*NE];
  extern const double PTL_F[NP*NP], PTL_FX[NP*NP], PTL_FY[NP*NP], PTL_FXY[NP*NP];
  extern const double PTV_F[NP*NP], PTV_FX[NP*NP], PTV_FY[NP*NP], PTV_FXY[NP*NP];
}
__device__ __forceinline__ double d_cr(double p0,double p1,double p2,double p3,double t){
  double a=-0.5*p0+1.5*p1-1.5*p2+0.5*p3, b=p0-2.5*p1+2.0*p2-0.5*p3, c=-0.5*p0+0.5*p2;
  return ((a*t+b)*t+c)*t+p1;
}
__device__ double d_bicubic(const double* __restrict__ TT, double rho, double e){
  using namespace prtab;
  double lr=log10(rho<1e-3?1e-3:rho);
  double fu=(lr-LR0)/(LR1-LR0)*(NLR-1), fv=(e-E0)/(E1-E0)*(NE-1);
  if(fu<1.0)fu=1.0; if(fu>NLR-2)fu=NLR-2; if(fv<1.0)fv=1.0; if(fv>NE-2)fv=NE-2;
  int iu=(int)fu, iv=(int)fv;
  if(iu>NLR-3)iu=NLR-3; if(iv>NE-3)iv=NE-3;
  double tu=fu-iu, tv=fv-iv, rT[4];
  for(int m=-1;m<=2;++m){ int r=iu+m; const double* row=TT+r*NE+iv;
    rT[m+1]=d_cr(__ldg(row-1),__ldg(row),__ldg(row+1),__ldg(row+2),tv); }
  return d_cr(rT[0],rT[1],rT[2],rT[3],tu);
}
__device__ __forceinline__ double d_h1(double f0,double f1,double d0,double d1,double t){
  double t2=t*t, t3=t2*t;
  return f0*(2*t3-3*t2+1)+f1*(3*t2-2*t3)+d0*(t3-2*t2+t)+d1*(t3-t2);
}
__device__ double d_patchT(const double* F,const double* FX,const double* FY,const double* FXY,double rho,double e){
  using namespace prtab;
  double lr=log10(rho<1e-3?1e-3:rho);
  double fu=(lr-PLR0)/(PLR1-PLR0)*(NP-1), fv=(e-PE0)/(PE1-PE0)*(NP-1);
  if(fu<0)fu=0; if(fu>NP-1)fu=NP-1; if(fv<0)fv=0; if(fv>NP-1)fv=NP-1;
  int i=(int)fu, j=(int)fv; if(i>NP-2)i=NP-2; if(j>NP-2)j=NP-2;
  double u=fu-i, v=fv-j; int a=i*NP+j, b=(i+1)*NP+j, c=i*NP+j+1, d=(i+1)*NP+j+1;
  double fj =d_h1(F[a],F[b],FX[a],FX[b],u), gj =d_h1(FY[a],FY[b],FXY[a],FXY[b],u);
  double fj1=d_h1(F[c],F[d],FX[c],FX[d],u), gj1=d_h1(FY[c],FY[d],FXY[c],FXY[d],u);
  return d_h1(fj,fj1,gj,gj1,v);
}
struct Tables { const double *TL,*TV,*PLF,*PLFX,*PLFY,*PLFXY,*PVF,*PVFX,*PVFY,*PVFXY; };
__device__ double d_table_T(const Tables& tb, int ph, double rho, double e){
  using namespace prtab;
  double T=d_bicubic(ph==1?tb.TL:tb.TV, rho, e);
  double lr=log10(rho<1e-3?1e-3:rho);
  if(HAVE_PATCH && lr>=PLR0 && lr<=PLR1 && e>=PE0 && e<=PE1)
    T=(ph==1)?d_patchT(tb.PLF,tb.PLFX,tb.PLFY,tb.PLFXY,rho,e):d_patchT(tb.PVF,tb.PVFX,tb.PVFY,tb.PVFXY,rho,e);
  return T;
}

// --- the fluid --------------------------------------------------------------
__constant__ hem::PRFluid c_fluid;

__device__ __forceinline__ hem::Phase3 d_br(int ph){ return ph==1?hem::Phase3::Liquid:hem::Phase3::Vapor; }
__device__ __forceinline__ double d_e_of_T(double T, double rho, int ph){
  return hem::state_from_T_v(c_fluid, T, c_fluid.M/rho, d_br(ph)).e;
}

// K1: CAMR's branch-locked inversion (bracketed Illinois, data-dependent trip count)
__global__ void k_pr_branch(int n, const double* rho, const double* e, const int* ph, double* T){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
  hem::State s; bool ok=hem::state_from_rho_e_phase_try(c_fluid, rho[i], e[i], d_br(ph[i]), s);
  T[i]= ok ? s.T : -1.0;
}
// K2: fixed-count bisection on e(T) (the GERG design, divergence-free)
__global__ void k_pr_bisect(int n, int iters, const double* rho, const double* e, const int* ph, double* T){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
  double lo=1.0, hi=5000.0, et=e[i], r=rho[i]; int p=ph[i];
  for(int k=0;k<iters;++k){ double m=0.5*(lo+hi); double em=d_e_of_T(m,r,p); if(em<et) lo=m; else hi=m; }
  T[i]=0.5*(lo+hi);
}
// K3: bicubic table lookup (PRTab design, table-terminal)
__global__ void k_table(int n, Tables tb, const double* rho, const double* e, const int* ph, double* T){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
  T[i]=d_table_T(tb, ph[i], rho[i], e[i]);
}
// Newton polish on the analytic surface from a seed (GERGTab design): k steps,
// derivative by one-sided difference, 2 e(T) evaluations per step.
__device__ __forceinline__ double d_polish(double T0, double rho, double et, int p, int steps){
  double T=T0;
  for(int k=0;k<steps;++k){ double h=1e-3*T; double e0=d_e_of_T(T,rho,p), e1=d_e_of_T(T+h,rho,p);
    double de=(e1-e0)/h; if(de>0.0) T-= (e0-et)/de; }
  return T;
}
__global__ void k_table_newton(int n, int steps, Tables tb, const double* rho, const double* e, const int* ph, double* T){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
  double T0=d_table_T(tb, ph[i], rho[i], e[i]);
  T[i]=d_polish(T0, rho[i], e[i], ph[i], steps);
}
// K4/K5: MLP of the trained size (6 inputs, W wide, D hidden layers, SiLU), weights in global memory.
template<typename R, int W, int D>
__global__ void k_mlp(int n, const R* __restrict__ wts, const double* rho, const double* e, const int* ph, double* T){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
  R x[6]; { double lr=log10(rho[i]); x[0]=(R)(lr/3.0); x[1]=(R)(e[i]/1e5);
    // stand-ins for the loci features (same arithmetic cost class as the real ones: a few flops each)
    x[2]=(R)(e[i]/1e5-0.5*lr); x[3]=(R)(e[i]/1e5-0.3*lr); x[4]=(R)sqrt(fabs(1.0-lr/3.0)); x[5]=(R)(ph[i]==1?1:-1); }
  R a[W], b[W]; const R* w=wts;
  // layer 0: 6 -> W
  for(int o=0;o<W;++o){ R s=w[o*7+6]; for(int k=0;k<6;++k) s+=w[o*7+k]*x[k]; a[o]=s/(R(1)+exp(-s)); }
  w+=W*7;
  for(int l=1;l<D;++l){
    for(int o=0;o<W;++o){ R s=w[o*(W+1)+W]; for(int k=0;k<W;++k) s+=w[o*(W+1)+k]*a[k]; b[o]=s/(R(1)+exp(-s)); }
    w+=W*(W+1);
    for(int k=0;k<W;++k) a[k]=b[k];
  }
  R s=w[W]; for(int k=0;k<W;++k) s+=w[k]*a[k];
  T[i]=(double)s*100.0;
}
template<typename R, int W, int D>
__global__ void k_mlp_newton(int n, int steps, const R* __restrict__ wts, const double* rho, const double* e, const int* ph, double* T){
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
  R x[6]; { double lr=log10(rho[i]); x[0]=(R)(lr/3.0); x[1]=(R)(e[i]/1e5);
    x[2]=(R)(e[i]/1e5-0.5*lr); x[3]=(R)(e[i]/1e5-0.3*lr); x[4]=(R)sqrt(fabs(1.0-lr/3.0)); x[5]=(R)(ph[i]==1?1:-1); }
  R a[W], b[W]; const R* w=wts;
  for(int o=0;o<W;++o){ R s=w[o*7+6]; for(int k=0;k<6;++k) s+=w[o*7+k]*x[k]; a[o]=s/(R(1)+exp(-s)); }
  w+=W*7;
  for(int l=1;l<D;++l){
    for(int o=0;o<W;++o){ R s=w[o*(W+1)+W]; for(int k=0;k<W;++k) s+=w[o*(W+1)+k]*a[k]; b[o]=s/(R(1)+exp(-s)); }
    w+=W*(W+1);
    for(int k=0;k<W;++k) a[k]=b[k];
  }
  R s=w[W]; for(int k=0;k<W;++k) s+=w[k]*a[k];
  // the net is untrained here: seed the polish from the reference-sized guess 300 K + the net's
  // (tiny, random) output so the Newton cost is measured on a realistic iterate, not the accuracy
  double T0=300.0+ (double)s*1e-3;
  T[i]=d_polish(T0, rho[i], e[i], ph[i], steps);
}

__global__ void k_probe(double rho, double e, int ph, double* out){
  hem::State s; bool ok=hem::state_from_rho_e_phase_try(c_fluid, rho, e, d_br(ph), s);
  out[0]=ok?1.0:0.0; out[1]=s.T; out[2]=d_e_of_T(280.0, rho, ph); out[3]=c_fluid.M; out[4]=0.0;
}
// --- host --------------------------------------------------------------------
#define CK(x) do{cudaError_t _e=(x); if(_e!=cudaSuccess){fprintf(stderr,"CUDA %s @%d: %s\n",#x,__LINE__,cudaGetErrorString(_e)); exit(1);} }while(0)

struct Row { int ph; double rho, e, T; };

template<typename F> double time_kernel(F launch, int reps=5){
  cudaEvent_t a,b; CK(cudaEventCreate(&a)); CK(cudaEventCreate(&b));
  launch(); CK(cudaGetLastError()); CK(cudaDeviceSynchronize());            // warm-up
  std::vector<float> ms;
  for(int r=0;r<reps;++r){ CK(cudaEventRecord(a)); launch(); CK(cudaEventRecord(b)); CK(cudaEventSynchronize(b));
    float m; CK(cudaEventElapsedTime(&m,a,b)); ms.push_back(m); }
  std::sort(ms.begin(),ms.end()); return ms[ms.size()/2];
}

int main(int argc, char** argv){
  int N = argc>1 ? atoi(argv[1]) : (1<<21);
  const char* path = argc>2 ? argv[2] : "data/deploy.csv";
  // ---- load the valid deployment rows
  std::vector<Row> rows; { FILE* f=fopen(path,"r"); if(!f){perror(path); return 1;} char line[4096]; fgets(line,sizeof line,f);
    while(fgets(line,sizeof line,f)){ int ph, ok; double rho,e,T,P,c; char rest[2048];
      if(sscanf(line,"%d,%lf,%lf,%lf,%lf,%lf,%d,%[^\n]",&ph,&rho,&e,&T,&P,&c,&ok,rest)>=7 && ok==1 && T>0) rows.push_back({ph,rho,e,T}); }
    fclose(f); }
  printf("# %zu valid deployment states from %s, tiled to N = %d\n", rows.size(), path, N);
  std::vector<double> h_rho(N), h_e(N), h_T(N); std::vector<int> h_ph(N);
  for(int i=0;i<N;++i){ const Row& r=rows[i%rows.size()]; h_rho[i]=r.rho; h_e[i]=r.e; h_ph[i]=r.ph; h_T[i]=r.T; }
  // sorted copy (by branch, then by the harvested T as a regime proxy) for the divergence test
  std::vector<int> order(N); for(int i=0;i<N;++i) order[i]=i;
  std::sort(order.begin(),order.end(),[&](int a,int b){ return h_ph[a]!=h_ph[b] ? h_ph[a]<h_ph[b] : h_T[a]<h_T[b]; });
  std::vector<double> s_rho(N), s_e(N); std::vector<int> s_ph(N);
  for(int i=0;i<N;++i){ s_rho[i]=h_rho[order[i]]; s_e[i]=h_e[order[i]]; s_ph[i]=h_ph[order[i]]; }

  double *d_rho,*d_e,*d_T,*d_T2,*d_srho,*d_se; int *d_ph,*d_sph;
  CK(cudaMalloc(&d_rho,N*8)); CK(cudaMalloc(&d_e,N*8)); CK(cudaMalloc(&d_T,N*8)); CK(cudaMalloc(&d_T2,N*8)); CK(cudaMalloc(&d_ph,N*4));
  CK(cudaMalloc(&d_srho,N*8)); CK(cudaMalloc(&d_se,N*8)); CK(cudaMalloc(&d_sph,N*4));
  CK(cudaMemcpy(d_rho,h_rho.data(),N*8,cudaMemcpyHostToDevice)); CK(cudaMemcpy(d_e,h_e.data(),N*8,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(d_ph,h_ph.data(),N*4,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(d_srho,s_rho.data(),N*8,cudaMemcpyHostToDevice)); CK(cudaMemcpy(d_se,s_e.data(),N*8,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(d_sph,s_ph.data(),N*4,cudaMemcpyHostToDevice));
  const hem::PRFluid CO2 = hem::PRFluid::make(304.13, 7.3773e6, 0.22394, 0.04401,
      2.35681300e+00, 8.98412990e-03, -7.12206320e-06, 2.45730071e-09, -1.42885166e-13, 0.0);
  CK(cudaMemcpyToSymbol(c_fluid,&CO2,sizeof CO2));
  // tables -> device
  auto up=[&](const double* src,size_t n){ double* d; CK(cudaMalloc(&d,n*8)); CK(cudaMemcpy(d,src,n*8,cudaMemcpyHostToDevice)); return (const double*)d; };
  using namespace prtab;
  Tables tb{ up(TBL_TL,(size_t)NLR*NE), up(TBL_TV,(size_t)NLR*NE),
             up(PTL_F,(size_t)NP*NP), up(PTL_FX,(size_t)NP*NP), up(PTL_FY,(size_t)NP*NP), up(PTL_FXY,(size_t)NP*NP),
             up(PTV_F,(size_t)NP*NP), up(PTV_FX,(size_t)NP*NP), up(PTV_FY,(size_t)NP*NP), up(PTV_FXY,(size_t)NP*NP) };
  printf("# tables: 2 branch tables of %d x %d doubles (%.1f MB) + 2 patches of %d^2 x 4 (%.2f MB); HAVE_PATCH=%d\n",
         NLR,NE,2.0*NLR*NE*8/1e6, NP, 8.0*NP*NP*8/1e6, (int)HAVE_PATCH);
  // random weights for the nets
  constexpr int W=128, D=4; size_t nw = W*7 + (size_t)(D-1)*W*(W+1) + (W+1);
  std::vector<double> hw(nw); srand(1); for(auto& v:hw) v=(rand()/(double)RAND_MAX-0.5)*0.2;
  std::vector<float> hwf(hw.begin(),hw.end());
  double* d_w64; float* d_w32; CK(cudaMalloc(&d_w64,nw*8)); CK(cudaMalloc(&d_w32,nw*4));
  CK(cudaMemcpy(d_w64,hw.data(),nw*8,cudaMemcpyHostToDevice)); CK(cudaMemcpy(d_w32,hwf.data(),nw*4,cudaMemcpyHostToDevice));
  printf("# net: 6 -> %d x %d -> 1, %zu weights (%.0f KB FP64, %.0f KB FP32)\n", W, D, nw, nw*8/1e3, nw*4/1e3);

  { double* d; CK(cudaMalloc(&d,5*8)); k_probe<<<1,1>>>(122.6251543, 88039.19965, 2, d); CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    double h[5]; CK(cudaMemcpy(h,d,5*8,cudaMemcpyDeviceToHost));
    printf("# probe: vapour (rho 122.63, e 88039) -> ok %.0f T %.6f K (host PR: 280.000)  e(280K) %.6e (host: 8.8039e4)  M %.5f\n", h[0],h[1],h[2],h[3]); }
  const int B=128, G=(N+B-1)/B;
  auto report=[&](const char* name, double ms, bool acc, bool fp32=false){
    double maxd=0, rms=0; int cnt=0;
    if(acc){ std::vector<double> t(N), ref(N); CK(cudaMemcpy(t.data(),d_T,N*8,cudaMemcpyDeviceToHost)); CK(cudaMemcpy(ref.data(),d_T2,N*8,cudaMemcpyDeviceToHost));
      for(int i=0;i<N;++i){ if(ref[i]>0){ double d=fabs(t[i]-ref[i]); maxd=std::max(maxd,d); rms+=d*d; ++cnt; } } rms=sqrt(rms/std::max(cnt,1)); }
    printf("%-34s %8.2f ms  %7.2f ns/state  %7.3f Gstate/s", name, ms, ms*1e6/N, N/ms/1e6);
    if(acc) printf("   vs branch inversion: rms %.2e K  max %.2e K", rms, maxd);
    printf("\n");
  };
  // reference: K1 on harvest order -> d_T2
  double ms1 = time_kernel([&]{ k_pr_branch<<<G,B>>>(N,d_rho,d_e,d_ph,d_T2); });
  report("K1 PR branch inversion (harvest order)", ms1, false);
  double ms1s= time_kernel([&]{ k_pr_branch<<<G,B>>>(N,d_srho,d_se,d_sph,d_T); });
  report("K1 PR branch inversion (sorted)", ms1s, false);
  for(int it : {32, 64}){
    char nm[64]; snprintf(nm,sizeof nm,"K2 fixed %d-step bisection", it);
    report(nm, time_kernel([&]{ k_pr_bisect<<<G,B>>>(N,it,d_rho,d_e,d_ph,d_T); }), true);
  }
  report("K3 bicubic table (+ dome patch)", time_kernel([&]{ k_table<<<G,B>>>(N,tb,d_rho,d_e,d_ph,d_T); }), true);
  report("K3s bicubic table (sorted)", time_kernel([&]{ k_table<<<G,B>>>(N,tb,d_srho,d_se,d_sph,d_T); }), false);
  if (argc > 3) {   // dump per-state table / seed+2Newton / reference T for the first tile (regime analysis offline)
    std::vector<double> t3(N), t3n(N), ref(N);
    k_table<<<G,B>>>(N,tb,d_rho,d_e,d_ph,d_T); CK(cudaDeviceSynchronize()); CK(cudaMemcpy(t3.data(),d_T,N*8,cudaMemcpyDeviceToHost));
    k_table_newton<<<G,B>>>(N,2,tb,d_rho,d_e,d_ph,d_T); CK(cudaDeviceSynchronize()); CK(cudaMemcpy(t3n.data(),d_T,N*8,cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(ref.data(),d_T2,N*8,cudaMemcpyDeviceToHost));
    FILE* o=fopen(argv[3],"w"); fprintf(o,"phase,rho,e,T_ref,T_table,T_table_newton2\n");
    for(size_t i=0;i<rows.size() && (int)i<N;++i) fprintf(o,"%d,%.9e,%.9e,%.9e,%.9e,%.9e\n",h_ph[i],h_rho[i],h_e[i],ref[i],t3[i],t3n[i]);
    fclose(o); printf("# dumped %zu states to %s\n", std::min(rows.size(),(size_t)N), argv[3]);
  }
  for(int st : {1, 2}){
    char nm[64]; snprintf(nm,sizeof nm,"K3+N table seed + %d Newton", st);
    report(nm, time_kernel([&]{ k_table_newton<<<G,B>>>(N,st,tb,d_rho,d_e,d_ph,d_T); }), true);
  }
  report("K4 MLP 128x4 FP64", time_kernel([&]{ k_mlp<double,W,D><<<G,B>>>(N,d_w64,d_rho,d_e,d_ph,d_T); }), false);
  report("K5 MLP 128x4 FP32", time_kernel([&]{ k_mlp<float,W,D><<<G,B>>>(N,d_w32,d_rho,d_e,d_ph,d_T); }), false);
  report("K5+N MLP FP32 seed + 2 Newton (cost)", time_kernel([&]{ k_mlp_newton<float,W,D><<<G,B>>>(N,2,d_w32,d_rho,d_e,d_ph,d_T); }), false);
  // smaller nets for the size scan
  report("K4b MLP 64x3 FP64", time_kernel([&]{ k_mlp<double,64,3><<<G,B>>>(N,d_w64,d_rho,d_e,d_ph,d_T); }), false);
  report("K5b MLP 64x3 FP32", time_kernel([&]{ k_mlp<float,64,3><<<G,B>>>(N,d_w32,d_rho,d_e,d_ph,d_T); }), false);
  report("K4c MLP 32x3 FP64", time_kernel([&]{ k_mlp<double,32,3><<<G,B>>>(N,d_w64,d_rho,d_e,d_ph,d_T); }), false);
  return 0;
}
