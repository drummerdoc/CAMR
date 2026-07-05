#include "IndexDefines.H"

void
init_pass_map(PassMap* pmap)
{
  int curMapIndx = 0;
  for (int i = 0; i < NUM_ADV; ++i) {
    pmap->upassMap[curMapIndx] = i + UFA;
    pmap->qpassMap[curMapIndx] = i + QFA;
    curMapIndx++;
  }
  for (int i = 0; i < NUM_SPECIES; ++i) {
    pmap->upassMap[curMapIndx] = i + UFS;
    pmap->qpassMap[curMapIndx] = i + QFS;
    curMapIndx++;
  }
  for (int i = 0; i < NUM_AUX; ++i) {
    pmap->upassMap[curMapIndx] = i + UFX;
    pmap->qpassMap[curMapIndx] = i + QFX;
    curMapIndx++;
  }
#ifdef USE_PS_HYDRO
  // Pelanti-Shyue extra components (α₁, α_k ρ_k, α_k ρ_k E_k) are
  // marked passive so Godunov/MOL passively advect them along the
  // contact until PS_umeth (Phase 4c) takes over.
  for (int i = 0; i < NUM_PS; ++i) {
    pmap->upassMap[curMapIndx] = i + UPS;
    pmap->qpassMap[curMapIndx] = i + QPS;
    curMapIndx++;
  }
#endif
}
