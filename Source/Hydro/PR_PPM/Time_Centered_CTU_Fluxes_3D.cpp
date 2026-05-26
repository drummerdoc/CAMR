#include "Time_Centered_CTU_Fluxes_3D.H"

#include "Compute_Face_Flux_3D_CTU.H"
#include "IndexDefines.H"

using namespace amrex;

void Time_Centered_CTU_Fluxes_3D(FArrayBox& fab_prim, const Box& bx,
								 GpuArray<FArrayBox, AMREX_SPACEDIM>& flux,
								 Real dt, const GpuArray<double, AMREX_SPACEDIM>& dx)
{
  const Box bxg1 = grow(bx, 1);
  const Box bxg2 = grow(bx, 2);
  auto const& prim = fab_prim.array();
  // -------------------------------------------------------------------------
  // STAGE 1: Pure 1D Flux Arrays (No transverse influence)
  // -------------------------------------------------------------------------
  const Box xflxbx = surroundingNodes(grow(bxg2, 0, -1), 0);
  FArrayBox flux_x_1d(xflxbx,NTHERM); flux_x_1d.setVal<RunOn::Host>(0.0); auto const& fx_1d = flux_x_1d.array();
  FArrayBox prim_x_1d(xflxbx,NTHERM); auto const& px_1d = prim_x_1d.array();
  ParallelFor(xflxbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	Compute_Face_Flux_3D_CTU(i,j,k, 1,0,0, 0,0,0, 0,0,0, UMX, UMY, UMZ, prim, fx_1d, px_1d, fx_1d, px_1d, fx_1d, px_1d, dt, dx);
  });

  const Box yflxbx = surroundingNodes(grow(bxg2, 1, -1), 1);
  FArrayBox flux_y_1d(yflxbx,NTHERM);  flux_y_1d.setVal<RunOn::Host>(0.0); auto const& fy_1d = flux_y_1d.array();
  FArrayBox prim_y_1d(yflxbx,NTHERM); auto const& py_1d = prim_y_1d.array();
  ParallelFor(yflxbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	Compute_Face_Flux_3D_CTU(i,j,k, 0,1,0, 0,0,0, 0,0,0, UMX, UMY, UMZ, prim, fy_1d, fy_1d, py_1d, fy_1d, py_1d, py_1d, dt, dx);
  });

  const Box zflxbx = surroundingNodes(grow(bxg2, 2, -1), 2);
  FArrayBox flux_z_1d(zflxbx,NTHERM);  flux_z_1d.setVal<RunOn::Host>(0.0); auto const& fz_1d = flux_z_1d.array();
  FArrayBox prim_z_1d(zflxbx,NTHERM); auto const& pz_1d = prim_z_1d.array();
  ParallelFor(zflxbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	Compute_Face_Flux_3D_CTU(i,j,k, 0,0,1, 0,0,0, 0,0,0, UMX, UMY, UMZ, prim, fz_1d, pz_1d, fz_1d, pz_1d, fz_1d, pz_1d, dt, dx);
  });

  // -------------------------------------------------------------------------
  // STAGE 2: Intermediate Cross-Flux Generation
  // -------------------------------------------------------------------------
  const Box xcflxbx = surroundingNodes(grow(bxg1, 0, -1), 0);
  FArrayBox flux_xy(xcflxbx,NTHERM);  auto const& f_xy = flux_xy.array();
  FArrayBox flux_xz(xcflxbx,NTHERM);  auto const& f_xz = flux_xz.array();
  FArrayBox prim_xy(xcflxbx,QTHERM);  auto const& p_xy = prim_xy.array();
  FArrayBox prim_xz(xcflxbx,QTHERM);  auto const& p_xz = prim_xy.array();
  ParallelFor(xcflxbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	// X-edges corrected by Y or Z 1D Flux cross gradients
	Compute_Face_Flux_3D_CTU(i,j,k, 1,0,0, 0,1,0, 0,0,0, UMX, UMY, UMZ, prim, fy_1d, py_1d, fx_1d, py_1d, f_xy, p_xy, dt, dx);
	Compute_Face_Flux_3D_CTU(i,j,k, 1,0,0, 0,0,1, 0,0,0, UMX, UMZ, UMY, prim, fz_1d, pz_1d, fx_1d, py_1d, f_xz, p_xz, dt, dx);
  });

  const Box ycflxbx = surroundingNodes(grow(bxg1, 1, -1), 1);
  FArrayBox flux_yx(ycflxbx,NTHERM);  auto const& f_yx = flux_yx.array();
  FArrayBox flux_yz(ycflxbx,NTHERM);  auto const& f_yz = flux_yz.array();
  FArrayBox prim_yx(ycflxbx,QTHERM);  auto const& p_yx = prim_yx.array();
  FArrayBox prim_yz(ycflxbx,QTHERM);  auto const& p_yz = prim_yz.array();
  ParallelFor(ycflxbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	// Y-edges corrected by X or Z 1D Flux cross gradients
	Compute_Face_Flux_3D_CTU(i,j,k, 0,1,0, 1,0,0, 0,0,0, UMY, UMX, UMZ, prim, fx_1d, px_1d, fy_1d, py_1d, f_yx, p_yx, dt, dx);
	Compute_Face_Flux_3D_CTU(i,j,k, 0,1,0, 0,0,1, 0,0,0, UMY, UMZ, UMX, prim, fz_1d, pz_1d, fy_1d, py_1d, f_yz, p_yz, dt, dx);
  });

  const Box zcflxbx = surroundingNodes(grow(bxg1, 2, -1), 2);
  FArrayBox flux_zx(zcflxbx,NTHERM);  auto const& f_zx = flux_zx.array();
  FArrayBox flux_zy(zcflxbx,NTHERM);  auto const& f_zy = flux_zy.array();
  FArrayBox prim_zx(zcflxbx,QTHERM);  auto const& p_zx = prim_zx.array();
  FArrayBox prim_zy(zcflxbx,QTHERM);  auto const& p_zy = prim_zy.array();
  ParallelFor(zcflxbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	// Z-edges corrected by X or Y 1D Flux cross gradients
	Compute_Face_Flux_3D_CTU(i,j,k, 0,0,1, 1,0,0, 0,0,0, UMZ, UMX, UMY, prim, fx_1d, px_1d, fz_1d, pz_1d, f_zx, p_zx, dt, dx);
	Compute_Face_Flux_3D_CTU(i,j,k, 0,0,1, 0,1,0, 0,0,0, UMZ, UMY, UMX, prim, fy_1d, py_1d, fz_1d, pz_1d, f_zy, p_zy, dt, dx);
  });

  // -------------------------------------------------------------------------
  // STAGE 3: Final 3D Unsplit Flux Construction & Grid Update
  // -------------------------------------------------------------------------
  const Box xbx = surroundingNodes(bx, 0); auto const& final_f_x = flux[0].array();
  FArrayBox final_prim_x(xbx,NTHERM);  auto const& final_p_x = final_prim_x.array();
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	Compute_Face_Flux_3D_CTU(i,j,k, 1,0,0, 0,1,0, 0,0,1, UMX, UMY, UMZ, prim, f_zy, p_zy, f_yz, p_yz, final_f_x, final_p_x, dt, dx);
  });

  const Box ybx = surroundingNodes(bx, 1); auto const& final_f_y = flux[1].array();
  FArrayBox final_prim_y(ybx,NTHERM);  auto const& final_p_y = final_prim_y.array();
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	Compute_Face_Flux_3D_CTU(i,j,k, 0,1,0, 1,0,0, 0,0,1, UMY, UMX, UMZ, prim, f_zx, p_zx, f_xz, p_xz, final_f_y, final_p_y, dt, dx);
  });

  const Box zbx = surroundingNodes(bx, 2); auto const& final_f_z = flux[2].array();
  FArrayBox final_prim_z(zbx,NTHERM);  auto const& final_p_z = final_prim_z.array();
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	Compute_Face_Flux_3D_CTU(i,j,k, 0,0,1, 1,0,0, 0,1,0, UMZ, UMX, UMY, prim, f_yx, p_yx, f_xy, p_xy, final_f_z, final_p_z, dt, dx);
  });
}
