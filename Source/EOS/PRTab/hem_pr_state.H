// PRTab forwarder — the real device-inline PR machinery lives in
// PR/.  Forwarding (rather than keeping a second physical copy)
// guarantees a SINGLE definition, so it cannot clash with
// (formerly PS_relax_device.H's explicit "EOS/PR/hem_pr_state.H" include;
//  that file was deleted 2026-08-17, S1)
// (both resolve to the same file; the include guard then dedupes).
// Resolved via the -I on Source/.
#include "EOS/PR/hem_pr_state.H"
