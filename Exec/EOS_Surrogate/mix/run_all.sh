#!/bin/bash
# Trace the CO2 + N2 phase envelopes, validate them, draw figure 7.  About one minute.
set -e
cd "$(dirname "$0")"
PY=${PY:-../../CO2_PipeBreak/yt_env/bin/python}
# the venv is untracked and lives in the main checkout; fall back to it from a worktree
[ -x "$PY" ] || PY=/home/marcusd/src/CAMR/Exec/CO2_PipeBreak/yt_env/bin/python
$PY envelope.py            # envelopes_z*.csv, envelopes_PT.csv, critical_points.csv
$PY validate.py | tee validate_output.txt   # gate + CoolProp table (+ coolprop_points.csv)
$PY fig7.py                # ../figures/fig7_phase_envelopes.png
