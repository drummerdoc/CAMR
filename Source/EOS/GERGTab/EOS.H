// =====================================================================
//  GERGTab = the GERG-2008 backend with the bicubic table fast path.
//  All physics lives in EOS/GERG/EOS.H; the USE_GERGTAB_EOS define
//  (set by Make.CAMR for Eos_Model := GERGTab) enables the table hooks
//  there, and this dir supplies gergtab_bicubic.H + generated tables.
//  Runtime gates: CAMR.eos_table, CAMR.eos_table_auto (default 1;
//  deprecated aliases eos_mlp / eos_mlp_auto honored).
// =====================================================================
#include "EOS/GERG/EOS.H"
