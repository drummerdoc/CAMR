#include <AMReX_ParmParse.H>

#include "CAMR.H"
#include "Tagging.H"

void
CAMR::read_tagging_params()
{
  amrex::ParmParse pp("tagging");

  pp.query("denerr", tagging_parm->denerr);
  pp.query("max_denerr_lev", tagging_parm->max_denerr_lev);
  pp.query("dengrad", tagging_parm->dengrad);
  pp.query("max_dengrad_lev", tagging_parm->max_dengrad_lev);

  pp.query("presserr", tagging_parm->presserr);
  pp.query("max_presserr_lev", tagging_parm->max_presserr_lev);
  pp.query("pressgrad", tagging_parm->pressgrad);
  pp.query("max_pressgrad_lev", tagging_parm->max_pressgrad_lev);

  pp.query("velerr", tagging_parm->velerr);
  pp.query("max_velerr_lev", tagging_parm->max_velerr_lev);
  pp.query("velgrad", tagging_parm->velgrad);
  pp.query("max_velgrad_lev", tagging_parm->max_velgrad_lev);

  pp.query("temperr", tagging_parm->temperr);
  pp.query("max_temperr_lev", tagging_parm->max_temperr_lev);
  pp.query("tempgrad", tagging_parm->tempgrad);
  pp.query("max_tempgrad_lev", tagging_parm->max_tempgrad_lev);

  //  AUDIT 2026-08-24 B14: tagging.alphaerr / tagging.max_alphaerr_lev were
  //  read into TaggingParm fields NO tagging routine consumed — a deck
  //  setting them got no α refinement and no warning.  Removed rather than
  //  wired: the LIVE route for α-based refinement is the errtagging path,
  //      amr.refinement_indicators = <name>
  //      amr.<name>.field_name     = ps_alpha1
  //  (+ value_greater / adjacent_difference_greater etc., CAMR_error.cpp),
  //  which the two-phase decks already use.
}
