CAMR_HOME      ?= ../..
AMREX_HOME     ?= ../../../amrex

USE_PS_HYDRO = TRUE

# AMReX
# DIM is conditionally-assigned so a `make DIM=1` on the command line
# builds the 1D variant (produces CAMR1d.*.PS.ex alongside 2D binaries
# in this same directory).  The T-Blowdown case is physically 1D — a
# closed pipe with vent at one end — so DIM=1 gives a ~32× wall-clock
# speedup at the same x-resolution for scheme-development iteration.
DIM ?= 2
COMP = gnu
PRECISION = DOUBLE

# Profiling
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
Eos_Model := PR

# GNU Make
Bpack := ./Make.package
Blocs := .

include $(CAMR_HOME)/Exec/Make.CAMR
