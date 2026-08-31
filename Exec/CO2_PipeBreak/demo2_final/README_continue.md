# Continuing the demo2 run from this checkpoint

`chk_sj2_00500` is the AMReX checkpoint at coarse step 500, t = 2.4753e-3 s
(the run's final steps 501-505 to t = 2.5e-3 are not checkpointed; the restart
recomputes them identically in ~2 minutes, so nothing is lost).

From `Exec/CO2_PipeBreak/`, with the 2-D exe built as usual:

    ./CAMR2d.gnu.TPROF.PS.PR.ex inputs.satjet_demo2 \
        amr.restart=demo2_final/chk_sj2_00500 \
        stop_time=<NEW_END_TIME> \
        max_step=1000000 \
        amr.check_int=50 amr.plot_int=10

Notes:
- stop_time: the deck says 2.5e-3; set your new end time (e.g. 1.0e-2).
- max_step: the deck caps at 20000 total level-0 steps - raise it for a long run.
- Plotfiles/checkpoints continue the plt_sj2_/chk_sj2_ numbering in the
  directory you launch from; each plt is ~8 MB (every 10 coarse steps) and
  each chk ~8-10 MB (every 50), so budget disk for a long run or raise the
  intervals.
- Pace on one core at t~2.5 ms was ~2 coarse steps/9 min and slows as the
  refined region grows; MPI (CAMR2d.llvm.TPROF.MPI.PS.PR.ex) is the practical
  choice for a long extension.
- plt_sj2_00505 here is the final state at exactly t = 2.5e-3 (for analysis);
  run3_partial.log is the full concatenated log of the original run.
