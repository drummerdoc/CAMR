#!/usr/bin/env bash
# Remove the scratch that the harnesses in this directory write.
# Touches only the known run-directory prefixes and their leftovers; every
# tracked file (scripts, decks, refs/, characterization/) is left alone.
# Run on the host:  bash Exec/CO2_RiemannSuite/cleanup_session.sh
set -u
cd "$(dirname "$0")" || exit 1   # -> Exec/CO2_RiemannSuite

echo "[1/2] harness run directories ..."
#   vc_ vc4_ vcrp0_ vcrp1_ vp_ vz_ v4_ vcb12_   verify_canonical.py
#   ex_wp_                                    exact_suite.py
#   ch_                                       characterize.py --defaults
#   fs_                                       full_suite.py / characterize.py (matched config)
#   hem_                                      hem_limit.py
#   cvg* conv_data                            gen_convergence.py
#   std_ref                                   full_suite.py standalone CSV mirror
#   flash_red_runs                            flashing_front.py
#   *.old.*                                   AMReX renames of overwritten plotfiles
rm -rf vc_* vc4_* vcrp0_* vcrp1_* vp_* vz_* v4_* vcb12_* ex_wp_* ch_* fs_* hem_* \
       cvg* conv_data *.old.* std_ref flash_red_runs __pycache__

echo "[2/2] scratch logs ..."
rm -f *.log

echo "done."
echo "NOTE: do not 'git clean -fdX' here without looking: the generated EOS tables"
echo "      and the local build directory are gitignored and would go with it."
