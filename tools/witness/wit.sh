#!/bin/bash
# wit.sh <case> <tag> <reftag>  (case: N1 N1b N4 FF)
CAMR_HOME=${CAMR_HOME:-$(cd "$(dirname "$0")/../.." && pwd)}
c=$1; t=$2; r=$3; N=$CAMR_HOME/Exec/CO2_NSCBC; S=$CAMR_HOME/Exec/CO2_RiemannSuite
export CHUNK=3000; E=$N/exe_keep/nc_$t.ex; I=$N/inputs.base; R=$N/runs/audit
case $c in
N1)  $N/runchunk.sh $R N1_$t $E $I stop_time=4.0e-3 amr.plot_per=5.0e-4 amr.n_cell=256 prob.amp=1.0 prob.ic=1; $CAMR_HOME/tools/witness/cmpruns.sh $R N1_${r}_plt_ $R N1_${t}_plt_ | tail -3;;
N1b) $N/runchunk.sh $R N1b_$t $E $I stop_time=6.0e-3 amr.plot_per=5.0e-4 amr.n_cell=256 prob.amp=1.0 prob.ic=1 prob.two_phase=1 prob.a1_0=0.05; $CAMR_HOME/tools/witness/cmpruns.sh $R N1b_${r}_plt_ $R N1b_${t}_plt_ | tail -3;;
N4)  $N/runchunk.sh $R N4vent_$t $E $I stop_time=12.0e-3 amr.plot_per=1.0e-3 amr.n_cell=256 prob.ic=0 prob.p0=2.2e6 CAMR.ps_bc_nscbc_type_xlo=vent CAMR.ps_bc_nscbc_type_xhi=vent; $CAMR_HOME/tools/witness/cmpruns.sh $R N4vent_${r}_plt_ $R N4vent_${t}_plt_ | tail -3;;
FF)  cd $S; mkdir -p flashv_$t; $S/exe_keep/nc_$t.ex inputs CAMR.cfl=0.25 prob.alpha_trace=1.0e-6 prob.u_L=0 prob.u_R=0 stop_time=2.5e-4 amr.plot_int=-1 amr.plot_per=2.0e-5 amr.v=0 CAMR.v=1 CAMR.ps_diag_mass=1 geometry.prob_lo=0.0 geometry.prob_hi=1.0 amr.n_cell=256 prob.x_diaph=0.5 prob.phase_L=1 prob.T_L=320 prob.p_L=1.0e7 prob.x_qual_L=0.0 prob.phase_R=1 prob.T_R=320 prob.p_R=1.0e7 prob.x_qual_R=0.0 CAMR.lo_bc=SlipWall CAMR.hi_bc=Inflow prob.p_amb=1.0e5 CAMR.ps_bc_use_nscbc=1 CAMR.ps_bc_nscbc_type_xhi=vent CAMR.ps_bc_nscbc_probe=1 amr.plot_file=flashv_$t/v2s_ > flashv_$t/v2s.log 2>&1; echo FF rc=$?
     $CAMR_HOME/tools/witness/cmpruns.sh flashv_$r v2s_ flashv_$t v2s_ | tail -3; grep "PS-BCRX" flashv_$r/v2s.log | cmp - <(grep "PS-BCRX" flashv_$t/v2s.log) 2>/dev/null; ;;
esac
