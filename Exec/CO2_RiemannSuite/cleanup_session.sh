#!/usr/bin/env bash
# Remove this session's scratch (Riemann battery + #86/#88 diagnostics).
# SAFE by construction: only touches session-scratch prefixes; PRESERVES the
# c1_* regression baselines and all kept scripts / inputs / docs / PDFs.
# Run on the HOST:  bash Exec/CO2_RiemannSuite/cleanup_session.sh
set -u
cd "$(dirname "$0")" || exit 1   # -> Exec/CO2_RiemannSuite

echo "[1/3] scratch plotfile/run dirs (NOT the c1_ baselines) ..."
rm -rf plt_* rgr_* fs_* mt_std* mt_camr* rebase_B9_* b9band_* b7t_* b7c_* dchk_* \
       ic_A1_* off_A1_* chk0_A1_* exact_A1_* verifA1_* gw0* g_m0* std_ref __pycache__ \
       ex_wp_* ex_hllc_* s2_* f0c_* f0x_* f1f_* f1t_* f1g_* tmp_build_s*

echo "[2/3] scratch logs / temp images / one-off scripts ..."
rm -f *.log tdiag_*.txt _a1-*.png _b4-*.png _b7-*.png _pg-*.png \
      a1_adjudication.png decomp_jointtarget.png decomp_profiles.png mt_B2_compare.png \
      ab.py ab1.py ab_ripple.py mkplot.py mt_compare.py

echo "[3/3] CO2_PipeBreak scratch logs ..."
( cd ../CO2_PipeBreak 2>/dev/null && rm -f gate.log p3*.log pk.log )

echo "done. PRESERVED: c1_* baselines, inputs.decomp, full_suite.py, exact_suite.py,"
echo "  _stage2.py, ps_plotfile.py, err_vs_analytic.py, montage_ac.py,"
echo "  measure_ripple.py, *.pdf.  (run_ac_suite.py and a1_adjudication.py are"
echo "  deleted from the tree -- list refreshed 2026-08-24.)"
echo "NOTE: do NOT 'git clean -fdX' in this repo -- the c1_* baselines are gitignored"
echo "      and would be deleted."
