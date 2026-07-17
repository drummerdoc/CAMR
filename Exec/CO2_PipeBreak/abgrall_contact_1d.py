#!/usr/bin/env python3
"""
abgrall_contact_1d.py — minimal 1-D two-phase HLLC harness that reproduces the
CAMR PS "stiff-contact per-phase-energy" instability (the satjet slug ring) in
isolation, and is the development/verification vehicle for task #64 (Abgrall
contact per-phase-pressure reset).  See ../../LEARNINGS.md.

Model: 6-eq PS mixture slots [alpha1, m1=a1*rho1, m2=a2*rho2, rhou, UE1, UE2],
per-phase STIFFENED-GAS EOS (P_k=(g-1)rho e - g Pinf; e=(P+g Pinf)/((g-1)rho);
c^2=g(P+Pinf)/rho).  Liquid is stiff (large Pinf) to mimic CO2 near the dome —
the Abgrall pathology is EOS-agnostic, so SG is a faithful, invertible proxy for
the PR mechanism; the fix developed here ports directly to PS_HLLC.H.

Interface: HLLC single-valued Godunov state replicating PS_HLLC::ps_star_state
  E_k_star = E_k + (S_M - u)(S_M + P_mix/q_k),  q_k = rho_k (S_K - u)   [the BUG]
Conserved slots (m1,m2,rhou,UE1,UE2) updated by flux div; alpha1 upwinded by S_M.

CANONICAL TEST: uniform (P,u), only alpha1 (and per-phase densities) jump at a
contact, advected at u>0.  A pressure-oscillation-free scheme keeps P EXACTLY
flat.  Without the reset, stiff liquid rings P and it grows (Abgrall).  With
--reset, the contact carries no per-phase pressure jump and P stays flat.

  python3 abgrall_contact_1d.py            # baseline
  python3 abgrall_contact_1d.py --reset    # reset branch (WIP, see STATUS)
  python3 abgrall_contact_1d.py --steps 400 --n 200

STATUS (WIP — do not trust yet; findings recorded in ../../LEARNINGS.md):
  * With the MONOTONE stiffened-gas EOS, the pure advected uniform-(p,u) contact
    oscillation is BOUNDED and small (~1-2% of P0) and does NOT run away over
    600 steps.  Dilute liquid (a1=0.05) makes it SMALLER, not larger.
  * => this minimal test does NOT reproduce the strong CAMR slug ring (~63%).
    Conclusion: the CAMR ring is NOT the textbook advected-contact Abgrall
    oscillation; it needs the DYNAMIC per-phase state inconsistency (m1, alpha1,
    UE1 evolving under the split conserved-mass / deposited-energy / frozen-alpha
    operators) + continuous inlet forcing.  #64 should target OPERATOR
    CONSISTENCY between per-phase mass and per-phase energy transport, not only a
    contact-pressure reset.
  * The --reset branch currently makes it WORSE (bug: likely the hand-rolled
    non-conservative alpha update and/or the reset ke/energy bookkeeping).  Debug
    the alpha update first (replace with clean upwind of d_t a + u d_x a = 0).
  * NEXT: either (a) drive this harness (inlet BC + let rho1=m1/a1 evolve) to
    reproduce the ring, or (b) build the richer test directly in CAMR.
"""
import argparse
import numpy as np

# --- per-phase stiffened-gas EOS -------------------------------------------
# liquid (phase 1): stiff, c1~350 m/s at P0,rho1;  vapor (phase 2): ~ideal
G1, PINF1 = 1.19, 7.8e7
G2, PINF2 = 1.40, 0.0

def P_of(rho, e, g, pinf):   return (g-1.0)*rho*e - g*pinf
def e_of(rho, P, g, pinf):   return (P + g*pinf) / ((g-1.0)*rho)
def c2_of(rho, P, g, pinf):  return g*(P+pinf)/rho

def ctoprim(Q):
    a1 = np.clip(Q[0], 1e-9, 1-1e-9); a2 = 1.0-a1
    m1, m2, rhou, UE1, UE2 = Q[1], Q[2], Q[3], Q[4], Q[5]
    rho = m1+m2; u = rhou/rho
    rho1 = m1/a1; rho2 = m2/a2
    ke = 0.5*u*u
    e1 = UE1/m1 - ke; e2 = UE2/m2 - ke
    P1 = P_of(rho1, e1, G1, PINF1); P2 = P_of(rho2, e2, G2, PINF2)
    Pmix = a1*P1 + a2*P2
    Y1 = m1/rho; Y2 = m2/rho
    cf2 = Y1*c2_of(rho1, P1, G1, PINF1) + Y2*c2_of(rho2, P2, G2, PINF2)  # Wallis frozen
    return a1, a2, rho, u, rho1, rho2, e1, e2, P1, P2, Pmix, np.sqrt(np.maximum(cf2, 1.0))

def phys_flux_state(Q):
    a1, a2, rho, u, rho1, rho2, e1, e2, P1, P2, Pmix, c = ctoprim_scalar(Q)
    F = np.zeros(6)
    F[0] = a1*u                          # alpha (handled non-cons; kept for ref)
    F[1] = Q[1]*u; F[2] = Q[2]*u
    F[3] = Q[3]*u + Pmix
    F[4] = (Q[4] + a1*Pmix)*u
    F[5] = (Q[5] + a2*Pmix)*u
    return F

def ctoprim_scalar(Q):
    a1 = min(max(Q[0], 1e-9), 1-1e-9); a2 = 1.0-a1
    m1, m2, rhou, UE1, UE2 = Q[1], Q[2], Q[3], Q[4], Q[5]
    rho = m1+m2; u = rhou/rho; rho1 = m1/a1; rho2 = m2/a2
    ke = 0.5*u*u; e1 = UE1/m1-ke; e2 = UE2/m2-ke
    P1 = P_of(rho1, e1, G1, PINF1); P2 = P_of(rho2, e2, G2, PINF2)
    Pmix = a1*P1+a2*P2
    Y1 = m1/rho; Y2 = m2/rho
    cf2 = Y1*c2_of(rho1, P1, G1, PINF1)+Y2*c2_of(rho2, P2, G2, PINF2)
    return a1, a2, rho, u, rho1, rho2, e1, e2, P1, P2, Pmix, np.sqrt(max(cf2, 1.0))

def hllc_face(QL, QR, reset, dtol_a=1e-3, dtol_p=0.02):
    a1L,a2L,rhoL,uL,r1L,r2L,e1L,e2L,P1L,P2L,PmL,cL = ctoprim_scalar(QL)
    a1R,a2R,rhoR,uR,r1R,r2R,e1R,e2R,P1R,P2R,PmR,cR = ctoprim_scalar(QR)
    SL = min(uL-cL, uR-cR); SR = max(uL+cR, uR+cR)
    SM = ((PmR-PmL) + rhoL*uL*(SL-uL) - rhoR*uR*(SR-uR)) / (rhoL*(SL-uL) - rhoR*(SR-uR))
    if SL >= 0.0:  return phys_flux_state(QL)
    if SR <= 0.0:  return phys_flux_state(QR)
    # star side by contact sign
    if SM >= 0.0:
        Q,K = QL,'L'; SK=SL; a1K,a2K=a1L,a2L; uK=uL; r1K,r2K=r1L,r2L
        e1K,e2K=e1L,e2L; P1K,P2K,PmK=P1L,P2L,PmL; m1K,m2K=QL[1],QL[2]
        dA=abs(a1L-a1R); dP=abs(PmL-PmR)/max(PmL,1.0)
    else:
        Q,K = QR,'R'; SK=SR; a1K,a2K=a1R,a2R; uK=uR; r1K,r2K=r1R,r2R
        e1K,e2K=e1R,e2R; P1K,P2K,PmK=P1R,P2R,PmR; m1K,m2K=QR[1],QR[2]
        dA=abs(a1L-a1R); dP=abs(PmL-PmR)/max(PmR,1.0)
    rK = (SK-uK)/(SK-SM)
    m1s = m1K*rK; m2s = m2K*rK; rho_s = m1s+m2s
    q1 = r1K*(SK-uK); q2 = r2K*(SK-uK)
    ke_s = 0.5*SM*SM
    # per-phase specific TOTAL energy at star (the standard PS/HLLC formula = the bug)
    E1s = e1K + ke_s + (SM-uK)*(SM + PmK/q1) - 0.5*uK*uK   # careful: e1K is specific INTERNAL
    E2s = e2K + ke_s + (SM-uK)*(SM + PmK/q2) - 0.5*uK*uK
    # NOTE: E_kK(total) = e_kK + 0.5 uK^2 ; formula: E_ks = E_kK + (SM-uK)(SM+Pm/q_k)
    E1s = (e1K+0.5*uK*uK) + (SM-uK)*(SM + PmK/q1)
    E2s = (e2K+0.5*uK*uK) + (SM-uK)*(SM + PmK/q2)
    if reset and dA > dtol_a and dP < dtol_p:
        # ---- Abgrall contact reset (#64): pin per-phase pressure to upwind ----
        rho1_s = m1s/a1K; rho2_s = m2s/a2K       # alpha advects with contact
        e1s = e_of(rho1_s, P1K, G1, PINF1)       # EOS inversion at upwind per-phase P
        e2s = e_of(rho2_s, P2K, G2, PINF2)
        E1s = e1s + ke_s; E2s = e2s + ke_s
    Qs = np.array([a1K, m1s, m2s, rho_s*SM, m1s*E1s, m2s*E2s])
    return phys_flux_state(Qs)

def run(n=100, steps=200, cfl=0.4, reset=False):
    L = 1.0; dx = L/n; x = (np.arange(n)+0.5)*dx
    P0, u0 = 4.0e6, 50.0
    # uniform P,u ; alpha1 (and per-phase densities fixed by P0) jump at x=0.5
    a1 = np.where(x < 0.5, 0.9, 0.1)
    rho1 = 800.0*np.ones(n); rho2 = 100.0*np.ones(n)  # set by (P0,T); uniform
    e1 = e_of(rho1, P0, G1, PINF1); e2 = e_of(rho2, P0, G2, PINF2)
    a2 = 1-a1; m1 = a1*rho1; m2 = a2*rho2; rho = m1+m2
    ke = 0.5*u0*u0
    Q = np.zeros((6, n))
    Q[0]=a1; Q[1]=m1; Q[2]=m2; Q[3]=rho*u0; Q[4]=m1*(e1+ke); Q[5]=m2*(e2+ke)
    def wrap(i): return i % n
    for s in range(steps):
        # max wave speed
        cmax = 0.0
        for i in range(n):
            *_, u, _,_,_,_,_,_,_,c = ctoprim_scalar(Q[:,i]); cmax=max(cmax, abs(u)+c)
        dt = cfl*dx/cmax
        F = np.zeros((6, n))       # flux at face i-1/2
        SMf = np.zeros(n)
        for i in range(n):
            QL = Q[:, wrap(i-1)]; QR = Q[:, i]
            F[:, i] = hllc_face(QL, QR, reset)
            # face velocity for non-cons alpha (recompute SM cheaply via ctoprim avg)
            _,_,rL,uL,_,_,_,_,_,_,PmL,cL = ctoprim_scalar(QL)
            _,_,rR,uR,_,_,_,_,_,_,PmR,cR = ctoprim_scalar(QR)
            SL=min(uL-cL,uR-cR); SR=max(uL+cR,uR+cR)
            SMf[i]=((PmR-PmL)+rL*uL*(SL-uL)-rR*uR*(SR-uR))/(rL*(SL-uL)-rR*(SR-uR))
        Qn = Q.copy()
        for i in range(n):
            ip = wrap(i+1)
            # conserved slots
            for k in (1,2,3,4,5):
                Qn[k, i] = Q[k, i] - dt/dx*(F[k, ip] - F[k, i])
            # alpha1 non-conservative upwind by face velocity
            uf_l = SMf[i]; uf_r = SMf[ip]
            aup_l = Q[0, wrap(i-1)] if uf_l > 0 else Q[0, i]
            aup_r = Q[0, i]         if uf_r > 0 else Q[0, ip]
            Qn[0, i] = Q[0, i] - dt/dx*(uf_r*aup_r - uf_l*aup_l) + dt/dx*Q[0,i]*(uf_r-uf_l)
        Q = Qn
    # report P profile flatness
    P = np.array([ctoprim_scalar(Q[:, i])[10] for i in range(n)])
    Prel = (P.max()-P.min())/P0
    print("reset=%-5s  n=%d steps=%d :  P range = [%.4f, %.4f] bar   spread/P0 = %.3e" % (
        str(reset), n, steps, P.min()/1e5, P.max()/1e5, Prel))
    return Prel

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=100)
    ap.add_argument('--steps', type=int, default=200)
    ap.add_argument('--cfl', type=float, default=0.4)
    ap.add_argument('--reset', action='store_true')
    a = ap.parse_args()
    run(n=a.n, steps=a.steps, cfl=a.cfl, reset=a.reset)
