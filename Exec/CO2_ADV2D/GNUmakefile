# CO2_XC2D — 2D diagonal cross-critical CO2 Riemann problem.
# Copy of CO2_B4's build config (2D, gnu, PS hydro, PR EOS).
CAMR_HOME      ?= ../..
AMREX_HOME     ?= ../../../amrex

USE_PS_HYDRO = TRUE

DIM = 2
COMP = gnu
PRECISION = DOUBLE

PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
TINY_PROFILE = TRUE

USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

DEBUG = FALSE

USE_EB = FALSE

Eos_Model := PR

Bpack := ./Make.package
Blocs := .

include $(CAMR_HOME)/Exec/Make.CAMR
