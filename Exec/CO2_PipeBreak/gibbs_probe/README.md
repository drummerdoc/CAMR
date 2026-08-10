# Gibbs driving-force probes (PR vs GERG)

Standalone host tools, no AMReX, no run required.

    E=../../../Source/EOS
    g++ -O2 -std=c++17 -DHEM_NO_AMREX -I$E/PR   sat_pr.cpp   -o sat_pr
    g++ -O2 -std=c++17                -I$E/GERG sat_gerg.cpp -o sat_gerg
    g++ -O2 -std=c++17 -DHEM_NO_AMREX -I$E/PR   rt_pr.cpp    -o rt_pr
    g++ -O2 -std=c++17                -I$E/GERG rt_gerg.cpp  -o rt_gerg

`sat_*`  : calibration. Prints g_L - g_V on each EOS's own saturation curve,
           220-300 K. Must be ~0 -- it is, to 1e-9 relative (PR) and 1e-10
           (GERG), so neither Gibbs implementation is buggy.

`rt_*`   : stdin "rho1 T1 rho2 T2", stdout "g1 g2 h1 h2 s1 s2 e1 e2 ok".
           State is specified by (rho, T) NOT (rho, e): internal energy
           carries an arbitrary reference that differs between the two
           EOSs by ~204 kJ/kg, so (rho, e) does NOT name the same physical
           state in both.

`probe_*`: same but keyed on (rho, e). RETAINED ONLY as a cautionary
           example -- do not use it to compare the two EOSs.
