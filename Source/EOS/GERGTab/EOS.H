// =====================================================================
//  GERGTab: the GERG-2008 backend with a bicubic table seed for the
//  (rho,e) -> T inversion.  All physics lives in EOS/GERG/EOS.H; the
//  USE_GERGTAB_EOS define (set by Make.CAMR for Eos_Model := GERGTab)
//  enables the table hooks there, and this directory supplies
//  gergtab_bicubic.H and the generated tables.  Dials (documented at
//  their accessors in GERG/EOS.H): CAMR.eos_table, CAMR.eos_table_auto,
//  CAMR.eos_diag.
// =====================================================================
#include "EOS/GERG/EOS.H"
