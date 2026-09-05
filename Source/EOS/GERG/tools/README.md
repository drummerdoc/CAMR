# GERG host tools

Two standalone host programs over the guarded GERG-2008 CO2 surface (`../gerg_co2_guard.H`); no AMReX.

- `gergstats.cpp`: reads `rho e a1 rho1 e1 rho2 e2` rows on stdin and prints per-phase P, T, c, the mixture P and c, and the extension/clamp/in-dome flags for each cell (bulk classification of solver states against the guarded surface).
- `gerg_edge_cliff.cpp`: scans the vapour branch across `rhoV_edge(T)` at a fixed T (argument, default 230 K), first at fixed T and then through the fixed-e inversion the solver uses, printing P, c and dP/drho|T.

Compile with the host compiler: `g++ -O2 -std=c++17 -o gergstats gergstats.cpp` (same for `gerg_edge_cliff.cpp`).
