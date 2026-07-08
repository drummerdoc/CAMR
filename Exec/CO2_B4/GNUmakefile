CAMR_HOME      ?= ../..
AMREX_HOME     ?= ../../../amrex

USE_PS_HYDRO = TRUE

# AMReX
# DIM is conditionally-assigned so a `make DIM=1` on the command line
# builds the 1D variant (produces CAMR1d.*.PS.ex alongside 2D binaries
# in this same directory).  The B4-Cross-critical Riemann problem is
# physically 1D — a diaphragm-split tube — so DIM=1 gives a ~32×
# wall-clock speedup at the same x-resolution for scheme-development
# iteration on the R-star validation.
DIM ?= 2
COMP = gnu
PRECISION = DOUBLE

# Profiling
USE_PS_HYDRO=TRUE
PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
TINY_PROFILE = TRUE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE

# CAMR
USE_EB = FALSE
Eos_Model := RealFluidCO2

# GNU Make
Bpack := ./Make.package
Blocs := .

include $(CAMR_HOME)/Exec/Make.CAMR
