CAMR_HOME      ?= ../..
AMREX_HOME     ?= ../../../amrex

# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
#TINY_PROFILE = TRUE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
#DEBUG = TRUE

#TEST = TRUE
#USE_ASSERTION = TRUE

# CAMR
USE_EB = FALSE

# Phase 1 wiring:  RealFluidCO2 is a γ=1.4 clone of GammaLaw so this
# case is numerically identical to Exec/Sod/ and exercises only the
# new build-system routing.  Phase 2 replaces the shim with a real
# Peng-Robinson CO2 EOS; the Exec inputs will then move to CO2 units.
Eos_Model := RealFluidCO2

# GNU Make
Bpack := ./Make.package
Blocs := .

include $(CAMR_HOME)/Exec/Make.CAMR
