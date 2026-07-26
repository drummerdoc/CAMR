// PRTab forwarder — the real device-inline PR machinery lives in
// PR/.  Forwarding (rather than keeping a second physical copy)
// guarantees a SINGLE definition, so it cannot clash with
// PS_relax_device.H's explicit "EOS/PR/hem_pr_state.H" include
// (both resolve to the same file; the include guard then dedupes).
// Resolved via the -I on Source/.
#include "EOS/PR/hem_pr_state.H"
