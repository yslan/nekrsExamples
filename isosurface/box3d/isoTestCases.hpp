// =========================================================================== //
//  isoTestCases.hpp -- named analytic test fields for the isoSurface plugin.
//
//  Each maker builds the nodal field on THIS rank's mesh and returns, together
//  with it, the isolevel and a closed-form distance-to-the-level-set function.
//  That last part is what makes these tests rather than demos: the error is
//  measured against the analytic surface, not against another run.
//
//      auto tc = isoTest::sphere(mesh);              // named constructor
//      iso = std::make_unique<isoSurface>(mesh, tc.field(), prm(tc));
//      ...
//      auto err = isoTest::check(tc, iso->roots(0));  // max distance, root count
//
//  What each case is for (report section references):
//
//    plane     x + 2y - z = 0.5      degree 1 in x.  On a STRAIGHT element the
//                                    interpolant is exact and the roots must hit
//                                    machine precision; on a curved one the map
//                                    makes f(x(r)) degree N and it is still
//                                    exact.  The baseline sanity check.
//    sphere    |x-x0|^2 = R^2        curvature and a closed surface.  On a
//                                    curved element f(x(r)) is NOT a polynomial,
//                                    so the error floor is the interpolation
//                                    error -- report section 2.  The test
//                                    measures that floor and compares.
//    poly5     prod (r-r_i) = 0      degree 5 in the REFERENCE coordinate, so
//                                    exact for N >= 5 and degree-DEFICIENT for
//                                    N > 5.  Exercises deflation (section 3.3):
//                                    the one case where plain route A fails.
//    tangency  (|x-x0|^2-R^2)^2 = 0  |grad f| -> 0 on the surface.  Must be
//                                    FLAGGED (bit0), not silently returned:
//                                    sqrt(eps) is the best any method can do
//                                    (section 2.5).
//    tgv       Q-criterion of TGV    a REAL flow field, not a manufactured one.
//                                    Eight vortex tubes; genuine topology.
// =========================================================================== //
#if !defined(iso_test_cases_hpp_)
#define iso_test_cases_hpp_

#include "isoSurface.hpp"

// <algorithm> <functional> <limits> <memory> <string> come from v1.hpp ->
// platform.hpp -> nekrsSys.hpp.  <cmath> does not: nekrsSys.hpp includes the C
// header <math.h>, not <cmath>, so std::sqrt/std::fabs are not guaranteed.
#include <cmath>

namespace isoTest
{

struct TestCase {
  std::string name;
  std::string description;
  dfloat level = 0;
  std::vector<dfloat> f;  // nodal values on this rank
  int expectedRoots = -1; // -1 = no analytic expectation
  bool expectTangency = false;

  // distance from a physical point to the TRUE level set; NaN if unavailable
  std::function<double(double, double, double)> dist;

  // ---- companion field, for the milestone-5 interpolation check -----------
  // A case may carry an extra field defined on the SAME nodes, together with its
  // closed form.  Interpolating `aux` onto the surface and comparing against
  // `auxExact` there is a direct test of the sampling path -- and, for a smooth
  // non-polynomial like the TGV velocity, one whose error must sit at the
  // interpolation floor rather than at machine precision.
  //   aux      nodal values, component-major: [c * Nlocal + i]
  //   auxComp  number of components
  //   auxExact (x,y,z,c) -> the exact value of component c
  std::string auxName;
  int auxComp = 0;
  std::vector<dfloat> aux;
  std::function<double(double, double, double, int)> auxExact;

  isoSurface::field field() const
  {
    // host field: the caller passes f directly to the vector overload
    return isoSurface::field{name, {}};
  }
};

namespace detail
{
inline void nodes(mesh_t *mesh, std::vector<dfloat> &X, std::vector<dfloat> &Y, std::vector<dfloat> &Z)
{
  auto [hx, hy, hz] = mesh->xyzHost();
  X = hx;
  Y = hy;
  Z = hz;
}

// Midpoint of the field's GLOBAL range.  Used as the default isolevel so a case
// always has something to find, whatever mesh it is handed.
inline void bbox(const std::vector<dfloat> &X,
                 const std::vector<dfloat> &Y,
                 const std::vector<dfloat> &Z,
                 dfloat *lo,
                 dfloat *hi)
{
  lo[0] = hi[0] = X.empty() ? 0 : X[0];
  lo[1] = hi[1] = Y.empty() ? 0 : Y[0];
  lo[2] = hi[2] = Z.empty() ? 0 : Z[0];
  for (size_t i = 0; i < X.size(); ++i) {
    lo[0] = std::min(lo[0], X[i]);
    hi[0] = std::max(hi[0], X[i]);
    lo[1] = std::min(lo[1], Y[i]);
    hi[1] = std::max(hi[1], Y[i]);
    lo[2] = std::min(lo[2], Z[i]);
    hi[2] = std::max(hi[2], Z[i]);
  }
  MPI_Allreduce(MPI_IN_PLACE, lo, 3, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, hi, 3, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm());
}

inline dfloat midLevel(const std::vector<dfloat> &f)
{
  dfloat lo = f.empty() ? 0 : f[0], hi = lo;
  for (const auto v : f) {
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  MPI_Allreduce(MPI_IN_PLACE, &lo, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &hi, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm());
  return 0.5 * (lo + hi);
}

} // namespace detail

// --------------------------------------------------------------------------- //
// `level` defaults to NaN, meaning "pick the midpoint of the field's global
// range".  A hard-coded level is a trap: on a mesh that does not straddle it the
// case silently finds zero roots and looks like a pass.  Pass an explicit level
// only when you know the mesh.
inline TestCase plane(mesh_t *mesh, dfloat level = std::numeric_limits<dfloat>::quiet_NaN())
{
  TestCase tc;
  tc.name = "plane";
  tc.description = "x + 2y - z (degree 1: the interpolant is exact)";
  std::vector<dfloat> X, Y, Z;
  detail::nodes(mesh, X, Y, Z);
  tc.f.resize(mesh->Nlocal);
  for (dlong i = 0; i < mesh->Nlocal; ++i) {
    tc.f[i] = X[i] + 2 * Y[i] - Z[i];
  }
  tc.level = std::isnan(level) ? detail::midLevel(tc.f) : level;
  const double c = static_cast<double>(tc.level);
  tc.dist = [c](double x, double y, double z) { return std::fabs(x + 2 * y - z - c) / std::sqrt(6.0); };
  return tc;
}

inline TestCase sphere(mesh_t *mesh, dfloat R = 0.8, dfloat cx = 4.5, dfloat cy = 1.5, dfloat cz = 2.25)
{
  TestCase tc;
  tc.name = "sphere";
  tc.description = "|x-x0|^2 = R^2 (curved: the floor is the interpolation error)";
  tc.level = R * R;
  std::vector<dfloat> X, Y, Z;
  detail::nodes(mesh, X, Y, Z);
  tc.f.resize(mesh->Nlocal);
  for (dlong i = 0; i < mesh->Nlocal; ++i) {
    const dfloat dx = X[i] - cx, dy = Y[i] - cy, dz = Z[i] - cz;
    tc.f[i] = dx * dx + dy * dy + dz * dz;
  }
  const double x0 = cx, y0 = cy, z0 = cz, rr = R;
  tc.dist = [x0, y0, z0, rr](double x, double y, double z) {
    return std::fabs(std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0)) - rr);
  };
  return tc;
}

// Degree-5 polynomial in the REFERENCE coordinate r.  For N > 5 the top Legendre
// modes vanish, which is exactly the degree-deficient case: without deflation
// the comrade matrix divides by roundoff (report 3.3).  Every fibre in r has the
// same five roots, which makes the expected count exact and checkable.
inline TestCase poly5(mesh_t *mesh)
{
  TestCase tc;
  tc.name = "poly5";
  tc.description = "prod(r - r_i), degree 5 in r: degree deflation";
  tc.level = 0.0;
  const int Nq = mesh->Nq;
  const double ri[5] = {-0.9, -0.42, 0.07, 0.37, 0.88};
  tc.f.resize(mesh->Nlocal);
  for (dlong e = 0; e < mesh->Nelements; ++e) {
    for (int c = 0; c < Nq; ++c) {
      for (int b = 0; b < Nq; ++b) {
        for (int a = 0; a < Nq; ++a) {
          const double r = static_cast<double>(mesh->r[(c * Nq + b) * Nq + a]);
          double p = 1.0;
          for (int q = 0; q < 5; ++q) {
            p *= (r - ri[q]);
          }
          tc.f[e * mesh->Np + (c * Nq + b) * Nq + a] = static_cast<dfloat>(p);
        }
      }
    }
  }
  // five roots on every r-fibre, none on the s- or t-fibres
  tc.expectedRoots = 5 * Nq * Nq * static_cast<int>(mesh->Nelements);
  return tc; // no dist(): the level set is in reference space
}

inline TestCase tangency(mesh_t *mesh, dfloat R = 0.8, dfloat cx = 4.5, dfloat cy = 1.5, dfloat cz = 2.25)
{
  TestCase tc;
  tc.name = "tangency";
  tc.description = "(|x-x0|^2 - R^2)^2 = 0: |grad f| -> 0, must be FLAGGED";
  tc.level = 0.0;
  tc.expectTangency = true;
  std::vector<dfloat> X, Y, Z;
  detail::nodes(mesh, X, Y, Z);
  tc.f.resize(mesh->Nlocal);
  for (dlong i = 0; i < mesh->Nlocal; ++i) {
    const dfloat dx = X[i] - cx, dy = Y[i] - cy, dz = Z[i] - cz;
    const dfloat q = dx * dx + dy * dy + dz * dz - R * R;
    tc.f[i] = q * q;
  }
  const double x0 = cx, y0 = cy, z0 = cz, rr = R;
  tc.dist = [x0, y0, z0, rr](double x, double y, double z) {
    return std::fabs(std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0)) - rr);
  };
  return tc;
}

// Q-criterion of the Taylor-Green vortex, evaluated analytically from the
// standard initial condition.  The domain is mapped to [0,2pi]^3 so one mesh
// spans the whole TGV cell regardless of its physical extent.
inline TestCase tgv(mesh_t *mesh, dfloat level = 0.15)
{
  TestCase tc;
  tc.name = "tgv";
  tc.description = "Q-criterion of the Taylor-Green vortex (a real flow field)";
  tc.level = level;

  std::vector<dfloat> X, Y, Z;
  detail::nodes(mesh, X, Y, Z);

  // map the mesh bounding box onto [0,2pi]^3
  dfloat lo[3] = {X[0], Y[0], Z[0]}, hi[3] = {X[0], Y[0], Z[0]};
  for (dlong i = 0; i < mesh->Nlocal; ++i) {
    lo[0] = std::min(lo[0], X[i]);
    hi[0] = std::max(hi[0], X[i]);
    lo[1] = std::min(lo[1], Y[i]);
    hi[1] = std::max(hi[1], Y[i]);
    lo[2] = std::min(lo[2], Z[i]);
    hi[2] = std::max(hi[2], Z[i]);
  }
  MPI_Allreduce(MPI_IN_PLACE, lo, 3, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, hi, 3, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm());

  const double twopi = 2.0 * M_PI;
  // the box->[0,2pi]^3 map, captured so auxExact can reproduce it exactly
  const double ox = lo[0], oy = lo[1], oz = lo[2];
  const double sxk = twopi / std::max<double>(hi[0] - lo[0], 1e-30);
  const double syk = twopi / std::max<double>(hi[1] - lo[1], 1e-30);
  const double szk = twopi / std::max<double>(hi[2] - lo[2], 1e-30);

  // The velocity is ALREADY computed below to build Q -- it is the gradient
  // tensor's source -- so keeping it costs one store and turns this case into an
  // analytic check on interpolation (plan D).
  tc.auxName = "velocity";
  tc.auxComp = 3;
  tc.aux.resize(static_cast<size_t>(3) * mesh->Nlocal);
  tc.auxExact = [ox, oy, oz, sxk, syk, szk](double px, double py, double pz, int c) {
    const double x = sxk * (px - ox), y = syk * (py - oy), z = szk * (pz - oz);
    if (c == 0) {
      return std::sin(x) * std::cos(y) * std::cos(z);
    }
    if (c == 1) {
      return -std::cos(x) * std::sin(y) * std::cos(z);
    }
    return 0.0;
  };

  tc.f.resize(mesh->Nlocal);
  for (dlong i = 0; i < mesh->Nlocal; ++i) {
    const double x = twopi * (X[i] - lo[0]) / std::max<double>(hi[0] - lo[0], 1e-30);
    const double y = twopi * (Y[i] - lo[1]) / std::max<double>(hi[1] - lo[1], 1e-30);
    const double z = twopi * (Z[i] - lo[2]) / std::max<double>(hi[2] - lo[2], 1e-30);
    // u =  sin x cos y cos z, v = -cos x sin y cos z, w = 0
    const double cx = std::cos(x), sx = std::sin(x);
    const double cy = std::cos(y), sy = std::sin(y);
    const double cz = std::cos(z), sz = std::sin(z);
    tc.aux[i] = static_cast<dfloat>(sx * cy * cz);
    tc.aux[mesh->Nlocal + i] = static_cast<dfloat>(-cx * sy * cz);
    tc.aux[2 * mesh->Nlocal + i] = static_cast<dfloat>(0.0);
    double g[3][3] = {{cx * cy * cz, -sx * sy * cz, -sx * cy * sz},
                      {sx * sy * cz, -cx * cy * cz, cx * sy * sz},
                      {0.0, 0.0, 0.0}};
    double Q = 0.0;
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        const double S = 0.5 * (g[a][b] + g[b][a]);
        const double O = 0.5 * (g[a][b] - g[b][a]);
        Q += 0.5 * (O * O - S * S);
      }
    }
    tc.f[i] = static_cast<dfloat>(Q);
  }
  return tc;
}

// --------------------------------------------------------------------------- //
//  check -- measure the error of an extracted RootSet against the analytic
//  surface, and report it on rank 0.  Returns the max distance (or NaN when the
//  case has no closed form, e.g. poly5 and tgv).
// --------------------------------------------------------------------------- //
struct CheckResult {
  double maxDist = 0.0;
  long long nRoots = 0;
  long long nTangency = 0;
  long long nOnFace = 0;
  long long nMultiEdge = 0;
};

inline CheckResult check(const TestCase &tc, const RootSet &rs, bool print = true)
{
  CheckResult cr;
  cr.nRoots = rs.nRoots;
  double worst = 0.0;
  for (dlong i = 0; i < rs.nRoots; ++i) {
    if (tc.dist) {
      worst = std::max(
          worst,
          tc.dist(static_cast<double>(rs.x[i]), static_cast<double>(rs.y[i]), static_cast<double>(rs.z[i])));
    }
    if (rs.flag[i] & RootSet::Tangency) {
      ++cr.nTangency;
    }
    if (rs.flag[i] & RootSet::OnFace) {
      ++cr.nOnFace;
    }
    if (rs.flag[i] & RootSet::MultiEdge) {
      ++cr.nMultiEdge;
    }
  }
  cr.maxDist = worst;

  MPI_Allreduce(MPI_IN_PLACE, &cr.nRoots, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nTangency, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nOnFace, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nMultiEdge, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.maxDist, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm());

  if (print && platform->comm.mpiRank() == 0) {
    printf("  isoTest[%-8s] roots=%-8lld", tc.name.c_str(), cr.nRoots);
    if (tc.dist) {
      printf(" maxDist=%.3e", cr.maxDist);
    } else {
      printf(" maxDist=   n/a  ");
    }
    if (tc.expectedRoots >= 0) {
      printf(" expected=%-8d %s", tc.expectedRoots, (cr.nRoots == tc.expectedRoots) ? "MATCH" : "MISMATCH");
    }
    if (cr.nTangency || tc.expectTangency) {
      printf(" tangency=%lld%s",
             cr.nTangency,
             (tc.expectTangency && cr.nTangency == 0) ? " <-- EXPECTED SOME" : "");
    }
    if (cr.nMultiEdge) {
      printf(" multiEdge=%lld", cr.nMultiEdge);
    }
    if (cr.nOnFace) {
      printf(" onFace=%lld", cr.nOnFace);
    }
    printf("\n");
  }
  return cr;
}

// --------------------------------------------------------------------------- //
//  checkSurface -- the milestone 2/3 error.
//
//  Measured INSIDE the triangles, not at their nodes.  Sampling only the nodes
//  would flatter every cell type equally: the corners are roots, hence exact by
//  construction, and a P2 cell's midside nodes are projected onto the surface,
//  so they are exact too.  What a flat triangulation actually gets wrong is the
//  INTERIOR, which is the whole reason milestone 3 exists -- so the samples are
//  taken at interior barycentric points, pushed through the cell's own shape
//  functions (the same ones VTK uses to render it).
// --------------------------------------------------------------------------- //
struct SurfaceCheckResult {
  double maxDist = 0.0; // worst in-triangle distance to the true surface
  double maxNode = 0.0; // worst distance at a NODE, for contrast
  long long nCell = 0;
  long long nPoint = 0;
  long long nCorner = 0;
  long long nSample = 0;
};

inline SurfaceCheckResult checkSurface(const TestCase &tc, const SurfaceSet &ss, bool print = true)
{
  SurfaceCheckResult cr;
  cr.nCell = ss.nCell;
  cr.nPoint = ss.nPoint;
  cr.nCorner = ss.nCorner;

  // interior barycentric samples; none on an edge or at a vertex
  static const double bc[7][3] = {{0.25, 0.25, 0.50},
                                  {0.50, 0.25, 0.25},
                                  {0.25, 0.50, 0.25},
                                  {1.0 / 3, 1.0 / 3, 1.0 / 3},
                                  {0.10, 0.45, 0.45},
                                  {0.45, 0.10, 0.45},
                                  {0.45, 0.45, 0.10}};
  const int npc = ss.nodesPerCell;

  if (tc.dist) {
    for (dlong i = 0; i < ss.nPoint; ++i) {
      cr.maxNode = std::max(cr.maxNode,
                            tc.dist(static_cast<double>(ss.px[i]),
                                    static_cast<double>(ss.py[i]),
                                    static_cast<double>(ss.pz[i])));
    }
    for (dlong c = 0; c < ss.nCell; ++c) {
      for (int q = 0; q < 7; ++q) {
        const double L1 = bc[q][0], L2 = bc[q][1], L3 = bc[q][2];
        double w[6];
        if (npc == 3) {
          w[0] = L1;
          w[1] = L2;
          w[2] = L3;
        } else { // VTK_QUADRATIC_TRIANGLE shape functions, VTK node order
          w[0] = L1 * (2 * L1 - 1);
          w[1] = L2 * (2 * L2 - 1);
          w[2] = L3 * (2 * L3 - 1);
          w[3] = 4 * L1 * L2;
          w[4] = 4 * L2 * L3;
          w[5] = 4 * L3 * L1;
        }
        double P[3] = {0, 0, 0};
        for (int m = 0; m < npc; ++m) {
          const dlong p = ss.conn[npc * c + m];
          P[0] += w[m] * static_cast<double>(ss.px[p]);
          P[1] += w[m] * static_cast<double>(ss.py[p]);
          P[2] += w[m] * static_cast<double>(ss.pz[p]);
        }
        cr.maxDist = std::max(cr.maxDist, tc.dist(P[0], P[1], P[2]));
        ++cr.nSample;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &cr.nCell, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nPoint, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nCorner, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nSample, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.maxDist, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.maxNode, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm());

  if (print && platform->comm.mpiRank() == 0) {
    printf("  isoSrf [%-8s] order=%d cells=%-7lld points=%-7lld (%lld corners + %lld mid",
           tc.name.c_str(),
           ss.cellOrder,
           cr.nCell,
           cr.nPoint,
           cr.nCorner,
           cr.nPoint - cr.nCorner);
    if (cr.nCorner) {
      printf(", %.2f/corner", double(cr.nPoint - cr.nCorner) / double(cr.nCorner));
    }
    printf(")");
    if (tc.dist) {
      printf(" inTri=%.3e node=%.3e", cr.maxDist, cr.maxNode);
    }
    printf("\n");
  }
  return cr;
}

// --------------------------------------------------------------------------- //
//  checkFields -- the milestone-5 error.
//
//  Compares the INTERPOLATED companion field against its closed form at every
//  surface point, corners and P2 midside nodes alike, and reports the two groups
//  separately: they exercise different paths (a corner inherits its (r,s,t) from
//  a root, a midside node gets it from the Newton projection), so a bug in one
//  is invisible if the two are averaged together.
//
//  For a smooth NON-polynomial like the TGV velocity the error must sit at the
//  interpolation floor of the discretisation, not at machine precision: the
//  quantity being measured is the interpolant, not the arithmetic.
// --------------------------------------------------------------------------- //
struct FieldCheckResult {
  double maxCorner = 0.0;
  double maxMid = 0.0;
  double refScale = 0.0; // max |exact|, so the error can be read as relative
  long long nCorner = 0;
  long long nMid = 0;
  bool ran = false;
};

inline FieldCheckResult
checkFields(const TestCase &tc, const SurfaceSet &ss, const FieldSet &fs, bool print = true)
{
  FieldCheckResult cr;
  if (!tc.auxExact || tc.auxComp <= 0) {
    return cr; // this case carries no companion field
  }

  int fi = -1;
  for (int f = 0; f < fs.nField(); ++f) {
    if (fs.name[f] == tc.auxName) {
      fi = f;
      break;
    }
  }
  if (fi < 0) {
    return cr; // the field was not registered for output
  }
  cr.ran = true;

  const int nc = fs.nComp[fi];
  for (dlong i = 0; i < ss.nPoint; ++i) {
    const double px = static_cast<double>(ss.px[i]);
    const double py = static_cast<double>(ss.py[i]);
    const double pz = static_cast<double>(ss.pz[i]);
    double e = 0.0;
    for (int c = 0; c < nc; ++c) {
      const double ex = tc.auxExact(px, py, pz, c);
      const double got = static_cast<double>(fs.val[fi][static_cast<size_t>(i) * nc + c]);
      e = std::max(e, std::fabs(got - ex));
      cr.refScale = std::max(cr.refScale, std::fabs(ex));
    }
    if (i < ss.nCorner) {
      cr.maxCorner = std::max(cr.maxCorner, e);
      ++cr.nCorner;
    } else {
      cr.maxMid = std::max(cr.maxMid, e);
      ++cr.nMid;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &cr.maxCorner, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.maxMid, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.refScale, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nCorner, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &cr.nMid, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());

  if (print && platform->comm.mpiRank() == 0) {
    printf("  isoFld [%-8s] %s(%d): corner %lld err=%.3e | midside %lld err=%.3e | scale %.3g\n",
           tc.name.c_str(),
           tc.auxName.c_str(),
           nc,
           cr.nCorner,
           cr.maxCorner,
           cr.nMid,
           cr.maxMid,
           cr.refScale);
  }
  return cr;
}

// --------------------------------------------------------------------------- //
//  writeCheckpoint -- dump one test scalar to its own checkpoint file, so the
//  field the plugin actually saw can be opened next to the extracted roots.
//
//  Uses nrs_t::writeToFile, which takes exactly the field-list type the plugin
//  already uses:
//      writeToFile(name, time, mesh, {{"name", {o_field}}, ...})
//  One file per case, stem = "iso" + the case name ("isoplane", "isosphere",
//  ...), so the .fld series is self-describing.
//
//  The device buffer is allocated per call and released on return: these are
//  one-shot diagnostics at tstep 0, not a per-step output stream.
// --------------------------------------------------------------------------- //
template <typename NrsT>
inline void writeCheckpoint(NrsT *nrs, mesh_t *mesh, const TestCase &tc, double time = 0.0)
{
  auto o_f = platform->device.malloc<dfloat>(mesh->Nlocal);
  o_f.copyFrom(tc.f.data(), mesh->Nlocal);

  // The VARIABLE name must be one the writer recognises.  iofldNek accepts only
  // {mesh, velocity, pressure, temperature} and scalarNN (iofldNek.cpp:5) and
  // aborts on anything else -- "plane" is not a legal field name.  tavg solves
  // the same problem by renaming to scalarNN (tavg.cpp:173).  Here each case
  // gets its OWN file, so there is exactly one field per file and
  // "temperature" is the natural slot; the case identity lives in the FILE
  // name, which is unconstrained.
  std::vector<std::tuple<std::string, std::vector<deviceMemory<dfloat>>>> list;
  list.push_back({"temperature", {deviceMemory<dfloat>(o_f)}});

  nrs->writeToFile("iso" + tc.name, time, mesh, list, /* enforceOutXYZ */ true);

  if (platform->comm.mpiRank() == 0) {
    printf("  isoTest[%-8s] checkpoint -> iso%s0.f*  (field stored as 'temperature')\n",
           tc.name.c_str(),
           tc.name.c_str());
  }
}

} // namespace isoTest

#endif // iso_test_cases_hpp_
