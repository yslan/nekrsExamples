// =========================================================================== //
//  isoSurface.hpp -- isosurfaces from a spectral-element field by global
//  polynomial rootfinding.  HOST half; the kernels are in isoSurface.okl.
//
//      #ifdef __okl__
//      #include "isoSurface.okl"     // the kernels
//      #endif
//      #include "isoSurface.hpp"     // this file -- NOT inside a guard
//
//  See README.md to use it and DEVELOPER.md to modify it: data structures, the
//  API, the eight kernels, what stays on the host, and the capacity bounds.
//
//  A COPY.  The source is plugin/v6_{host.hpp,device.okl}, which is where the
//  tests are; a fix made here is a fix no test will see.
// =========================================================================== //
#if !defined(nekrs_isosurface_hpp_)
#define nekrs_isosurface_hpp_

// The development header plugin/v1.hpp declares a class of the SAME name, so a
// translation unit including both gets a redefinition with no hint as to why.
// Say why.  (Only plugin/ ever does this; a user's case will not.)
#if defined(nekrs_isosurface_v1_hpp_)
#error "isoSurface.hpp and plugin/v1.hpp both define class isoSurface -- include one"
#endif

#include "platform.hpp"
#include "mesh.h"

#include <limits>

// platform.hpp -> nekrsSys.hpp already provides <algorithm> <cstdint> <memory>
// <string> <vector> <cmath>.  These are not in it (they resolve transitively
// today, but that is what breaks on a compiler bump):
#include <complex>
#include <cstddef>
#include <fstream>       // ofstream, in writeVtm()
#include <iomanip>       // setw/setfill, in outName()
#include <sstream>       // ostringstream, in outName() and the VTU header
#include <unordered_map> // the P2 edge dedup

// =========================================================================== //
//  INLINED from plugin/isoKernels.hpp (+ plugin/isoAberthReal.h).
//  Pure functions on POD arrays: no NekRS types, no allocation, no I/O.
// =========================================================================== //
// The real-pair Aberth, shared verbatim with v3_device.okl.  Included
// at namespace scope (its symbols are `static inline`, prefixed `iso*`) so the
// host and the device run the SAME source text rather than two transcriptions
// that can drift.  See plans/E_occa_port_v2.md and plugin/isoAberthReal.h.
// --- INLINED from plugin/isoAberthReal.h ---------------------------------
#if !defined(ISO_AB_REAL)
#define ISO_AB_REAL double
#endif

//  Value and derivative of sum_k a[k] P_k(z) at complex z = (zr, zi).
//
//  Same three-term recurrence as the complex overload in isoKernels.hpp, with
//  every product written out:  (x+iy)(u+iv) = (xu-yv) + i(xv+yu).
static inline void isoLegendreValDerC(const ISO_AB_REAL *a,
                                      int n,
                                      ISO_AB_REAL zr,
                                      ISO_AB_REAL zi,
                                      ISO_AB_REAL *valr,
                                      ISO_AB_REAL *vali,
                                      ISO_AB_REAL *derr,
                                      ISO_AB_REAL *deri)
{
  if (n <= 0) {
    *valr = a[0];
    *vali = 0.0;
    *derr = 0.0;
    *deri = 0.0;
    return;
  }

  // P0 = 1, P1 = z ; D0 = 0, D1 = 1
  ISO_AB_REAL P0r = 1.0, P0i = 0.0;
  ISO_AB_REAL P1r = zr, P1i = zi;
  ISO_AB_REAL D0r = 0.0, D0i = 0.0;
  ISO_AB_REAL D1r = 1.0, D1i = 0.0;

  ISO_AB_REAL vr = a[0] + a[1] * P1r;
  ISO_AB_REAL vi = a[1] * P1i;
  ISO_AB_REAL dr = a[1] * D1r;
  ISO_AB_REAL di = a[1] * D1i;

  for (int k = 1; k < n; ++k) {
    // NB the parentheses around each cast are REQUIRED, not style.  OCCA's
    // expression parser cannot apply a binary operator whose right operand is a
    // cast-to-parenthesised-expression: `a / (T)(x + 1)` is read as a division
    // followed by a call and fails with "Unable to apply operator".  Wrapping
    // the cast -- `a / ((T)(x + 1))` -- parses.  Verified against the real
    // parser with `occa translate`.  See notes/22.
    const ISO_AB_REAL kden = ((ISO_AB_REAL)(k + 1));
    const ISO_AB_REAL c0 = ((ISO_AB_REAL)(2 * k + 1)) / kden;
    const ISO_AB_REAL c1 = ((ISO_AB_REAL)k) / kden;

    // P2 = c0 * z * P1 - c1 * P0
    const ISO_AB_REAL zP1r = zr * P1r - zi * P1i;
    const ISO_AB_REAL zP1i = zr * P1i + zi * P1r;
    const ISO_AB_REAL P2r = c0 * zP1r - c1 * P0r;
    const ISO_AB_REAL P2i = c0 * zP1i - c1 * P0i;

    // D2 = c0 * (P1 + z * D1) - c1 * D0
    const ISO_AB_REAL zD1r = zr * D1r - zi * D1i;
    const ISO_AB_REAL zD1i = zr * D1i + zi * D1r;
    const ISO_AB_REAL D2r = c0 * (P1r + zD1r) - c1 * D0r;
    const ISO_AB_REAL D2i = c0 * (P1i + zD1i) - c1 * D0i;

    vr += a[k + 1] * P2r;
    vi += a[k + 1] * P2i;
    dr += a[k + 1] * D2r;
    di += a[k + 1] * D2i;

    P0r = P1r;
    P0i = P1i;
    P1r = P2r;
    P1i = P2i;
    D0r = D1r;
    D0i = D1i;
    D1r = D2r;
    D1i = D2i;
  }

  *valr = vr;
  *vali = vi;
  *derr = dr;
  *deri = di;
}

//  Initial guess: Chebyshev points pushed off the real axis by delta.
//
//  delta != 0 is REQUIRED -- a purely real start cannot reach a complex root
//  (measured: 1161 failures at delta = 0 over 4000 draws, 0 for delta in
//  [0.01, 0.2]).  See notes/14; the value itself is mid-plateau, not tuned.
static inline void isoAberthInitR(int n, ISO_AB_REAL delta, ISO_AB_REAL *zr, ISO_AB_REAL *zi)
{
  for (int i = 0; i < n; ++i) {
    const ISO_AB_REAL th = 3.14159265358979323846 * (((ISO_AB_REAL)i) + 0.5) / ((ISO_AB_REAL)n);
    zr[i] = cos(th);
    zi[i] = delta * sin(th);
  }
}

// --------------------------------------------------------------------------- //
//  The iteration itself.  FIXED niter, no convergence test, no early exit:
//  identical work per fibre, so a GPU warp never diverges (notes/10, notes/14).
//
//      w      = p(z_i) / p'(z_i)
//      S      = sum_{j != i} 1 / (z_i - z_j)      the repulsion term
//      z_i   -= w / (1 - w S)
//
//  *finalCorr (may be null) returns the largest final correction -- the safety
//  net that says whether niter was enough.
// --------------------------------------------------------------------------- //
static inline void isoAberthIterateR(const ISO_AB_REAL *a,
                                     int n,
                                     int niter,
                                     ISO_AB_REAL *zr,
                                     ISO_AB_REAL *zi,
                                     ISO_AB_REAL *finalCorr)
{
  ISO_AB_REAL corrMax = 0.0;
  for (int it = 0; it < niter; ++it) {
    corrMax = 0.0;
    for (int i = 0; i < n; ++i) {
      ISO_AB_REAL vr, vi, dr, di;
      isoLegendreValDerC(a, n, zr[i], zi[i], &vr, &vi, &dr, &di);

      // w = v / d   (complex divide, guarded exactly as the complex version:
      // |d| > 0, i.e. the denominator is not an exact zero)
      const ISO_AB_REAL d2 = dr * dr + di * di;
      ISO_AB_REAL wr = 0.0, wi = 0.0;
      if (d2 > 0.0) {
        wr = (vr * dr + vi * di) / d2;
        wi = (vi * dr - vr * di) / d2;
      }

      // S = sum_{j != i} 1 / (z_i - z_j)
      ISO_AB_REAL Sr = 0.0, Si = 0.0;
      for (int j = 0; j < n; ++j) {
        if (j == i) {
          continue;
        }
        const ISO_AB_REAL xr = zr[i] - zr[j];
        const ISO_AB_REAL xi = zi[i] - zi[j];
        const ISO_AB_REAL x2 = xr * xr + xi * xi;
        if (x2 > 0.0) {
          Sr += xr / x2;
          Si += -xi / x2;
        }
      }

      // den = 1 - w S
      const ISO_AB_REAL wSr = wr * Sr - wi * Si;
      const ISO_AB_REAL wSi = wr * Si + wi * Sr;
      const ISO_AB_REAL denr = 1.0 - wSr;
      const ISO_AB_REAL deni = -wSi;

      // corr = w / den, falling back to w when den underflows
      const ISO_AB_REAL den2 = denr * denr + deni * deni;
      ISO_AB_REAL cr, ci;
      if (den2 > 1e-300) {
        cr = (wr * denr + wi * deni) / den2;
        ci = (wi * denr - wr * deni) / den2;
      } else {
        cr = wr;
        ci = wi;
      }
      // a non-finite correction means the iteration has left the rails; freeze
      // this root rather than poison the whole fibre with a NaN
      if (!(cr == cr) || !(ci == ci) || (cr - cr) != 0.0 || (ci - ci) != 0.0) {
        cr = 0.0;
        ci = 0.0;
      }

      zr[i] -= cr;
      zi[i] -= ci;

      const ISO_AB_REAL ac = sqrt(cr * cr + ci * ci);
      corrMax = (ac > corrMax) ? ac : corrMax;
    }
  }
  if (finalCorr) {
    *finalCorr = corrMax;
  }
}

// --- end of the inlined isoAberthReal.h ----------------------------------

namespace isoDetail
{

//  Compile-time limits.  N <= 20 covers every spectral-element order in use and
//  lets every buffer live on the stack, which is what an OCCA port needs.
static constexpr int maxN = 20;
static constexpr int maxNq = maxN + 1;

//  K_GLL[N] -- the EXACT worst-case over/undershoot of a degree-N GLL
//  interpolant, as a fraction of the range of its own nodal values:
//
//      K_N = max_x sum_j max(L_j(x),0) - 1 = (Lambda_N - 1)/2
//
//  ATTAINABLE, not a bound, so culling on the plain nodal min/max is not
//  rigorous and can drop a thin crossing.  Indexed by ORDER: K_GLL[8] is N=8.
//  notes/15 derives it.
static constexpr double K_GLL[maxNq] = {0.000000, 0.000000, 0.125000, 0.250000, 0.317941, 0.389297, 0.436872,
                                        0.486184, 0.522820, 0.560482, 0.590271, 0.620736, 0.645834, 0.671410,
                                        0.693094, 0.715132, 0.734219, 0.753579, 0.770624, 0.787886, 0.803284};

// K_N for any order; the fit is accurate to < 0.008 over N = 4..20.
inline double gllOvershoot(int N)
{
  if (N < 0) {
    return 0.0;
  }
  if (N < maxNq) {
    return K_GLL[N];
  }
  return 0.2932 * std::log(static_cast<double>(N)) - 0.0821;
}

//  Element culling.  An element can contain a crossing of `level` only if
//
//      min_j f_j - theta*K_N*R  <=  level  <=  max_j f_j + theta*K_N*R,
//      R = max_j f_j - min_j f_j
//
//  theta = 1 is rigorous and the default; smaller culls more but trades away the
//  guarantee (notes/15).
template <typename T> inline bool elementMayCross(const T *f, int Np, double level, int N, double theta)
{
  if (Np <= 0) {
    return false;
  }
  double lo = static_cast<double>(f[0]);
  double hi = lo;
  for (int i = 1; i < Np; ++i) {
    const double v = static_cast<double>(f[i]);
    lo = (v < lo) ? v : lo;
    hi = (v > hi) ? v : hi;
  }
  const double margin = theta * gllOvershoot(N) * (hi - lo);
  return (level >= lo - margin) && (level <= hi + margin);
}

//  Legendre series: value and derivative by the three-term recurrence.
//
//      (k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}
//
//  run FORWARDS as an evaluation.  O(n) work, O(1) storage, no matrix -- this is
//  the kernel primitive of route C (report section 5.1).  Complex argument,
//  because Aberth iterates in the complex plane.
inline void legendreValDer(const double *a,
                           int n,
                           const std::complex<double> &z,
                           std::complex<double> &val,
                           std::complex<double> &der)
{
  if (n <= 0) {
    val = std::complex<double>(a[0], 0.0);
    der = std::complex<double>(0.0, 0.0);
    return;
  }
  std::complex<double> P0(1.0, 0.0), P1 = z;
  std::complex<double> D0(0.0, 0.0), D1(1.0, 0.0);
  val = a[0] * P0 + a[1] * P1;
  der = a[1] * D1;
  for (int k = 1; k < n; ++k) {
    const double c0 = static_cast<double>(2 * k + 1) / static_cast<double>(k + 1);
    const double c1 = static_cast<double>(k) / static_cast<double>(k + 1);
    const std::complex<double> P2 = c0 * z * P1 - c1 * P0;
    const std::complex<double> D2 = c0 * (P1 + z * D1) - c1 * D0;
    val += a[k + 1] * P2;
    der += a[k + 1] * D2;
    P0 = P1;
    P1 = P2;
    D0 = D1;
    D1 = D2;
  }
}

// Real-argument overload: used for the transverse derivative and for residuals.
inline void legendreValDer(const double *a, int n, double x, double &val, double &der)
{
  if (n <= 0) {
    val = a[0];
    der = 0.0;
    return;
  }
  double P0 = 1.0, P1 = x, D0 = 0.0, D1 = 1.0;
  val = a[0] + a[1] * x;
  der = a[1];
  for (int k = 1; k < n; ++k) {
    const double c0 = static_cast<double>(2 * k + 1) / static_cast<double>(k + 1);
    const double c1 = static_cast<double>(k) / static_cast<double>(k + 1);
    const double P2 = c0 * x * P1 - c1 * P0;
    const double D2 = c0 * (P1 + x * D1) - c1 * D0;
    val += a[k + 1] * P2;
    der += a[k + 1] * D2;
    P0 = P1;
    P1 = P2;
    D0 = D1;
    D1 = D2;
  }
}

//  Degree deflation (report section 3.3).  Drop trailing Legendre coefficients
//  below aFloor * max|a|, so the leading coefficient the comrade/Aberth
//  machinery relies on is genuine rather than roundoff.  Returns the true
//  degree.  aFloor = 1e-12 is the measured safe plateau: over 20,000 randomized
//  fibres it gave 0 root-count mismatches against route B, where the naive
//  choice gave 24.
inline int deflateDegree(const double *a, int N, double aFloor)
{
  double amax = 0.0;
  for (int k = 0; k <= N; ++k) {
    const double v = std::fabs(a[k]);
    amax = (v > amax) ? v : amax;
  }
  if (amax == 0.0) {
    return 0;
  }
  const double cut = aFloor * amax;
  int d = 0;
  for (int k = 0; k <= N; ++k) {
    if (std::fabs(a[k]) > cut) {
      d = k;
    }
  }
  return d;
}

// --------------------------------------------------------------------------- //
//  ROUTE C: Ehrlich-Aberth on the Legendre series.  All n roots at once,
//
//      z_i <- z_i - w_i / (1 - w_i * S_i),   w_i = p(z_i)/p'(z_i),
//      S_i = sum_{j != i} 1/(z_i - z_j)
//
//  Newton on p(x)/prod_{j!=i}(x-z_j): S_i is a mutual repulsion, which stops two
//  estimates collapsing onto the same root.  Cubic.
//
//  INITIALIZATION is what breaks this method.  Chebyshev points nudged off the
//  real axis, delta = 0.1: delta must be NONZERO or the iteration stays real and
//  cannot reach the complex roots, which are the majority.  The textbook
//  circle-of-Cauchy-radius start diverges here (notes/14).
//
//  niter is FIXED -- cubic, then it STALLS at machine precision, so overshooting
//  costs flops and never accuracy.  12 = 10 with margin.
// --------------------------------------------------------------------------- //
inline void aberthInit(int n, double delta, std::complex<double> *z)
{
  for (int i = 0; i < n; ++i) {
    const double th = M_PI * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
    z[i] = std::complex<double>(std::cos(th), delta * std::sin(th));
  }
}

inline void
aberthIterate(const double *a, int n, int niter, std::complex<double> *z, double *finalCorr = nullptr)
{
  double corrMax = 0.0;
  for (int it = 0; it < niter; ++it) {
    corrMax = 0.0;
    for (int i = 0; i < n; ++i) {
      std::complex<double> v, d;
      legendreValDer(a, n, z[i], v, d);
      const std::complex<double> w = (std::abs(d) > 0.0) ? v / d : std::complex<double>(0.0, 0.0);

      std::complex<double> S(0.0, 0.0);
      for (int j = 0; j < n; ++j) {
        if (j == i) {
          continue;
        }
        const std::complex<double> diff = z[i] - z[j];
        if (std::abs(diff) > 0.0) {
          S += 1.0 / diff;
        }
      }

      const std::complex<double> den = 1.0 - w * S;
      std::complex<double> corr = (std::abs(den) > 1e-300) ? w / den : w;
      if (!std::isfinite(corr.real()) || !std::isfinite(corr.imag())) {
        corr = std::complex<double>(0.0, 0.0);
      }
      z[i] -= corr;
      const double ac = std::abs(corr);
      corrMax = (ac > corrMax) ? ac : corrMax;
    }
  }
  if (finalCorr) {
    *finalCorr = corrMax;
  }
}

// --------------------------------------------------------------------------- //
//  The full 1D solve: modal coefficients in, real roots inside the element out.
//
//  Filtering, both parts required -- Aberth returns all n roots, most of them
//  complex or outside this element:
//    * |Im z| < imagTol                 -> real
//    * -1-margin <= Re z <= 1+margin    -> inside.  The margin keeps a root on
//      an element face; the neighbour finds it too and the caller flags it.
//
//  slope[i] = p'(root_i) is free (the iteration forms p/p' anyway): its sign is
//  the crossing direction, its magnitude flags a tangency.  Returns the count.
// --------------------------------------------------------------------------- //
struct AberthOpts {
  int niter = 12;
  double delta = 0.1;
  double aFloor = 1e-12;
  double imagTol = 1e-7;
  double margin = 1e-10;
};

inline int solveFibre(const double *a_in,
                      int N,
                      const AberthOpts &opt,
                      double *roots,
                      double *slopes,
                      double *finalCorr = nullptr)
{
  const int n = deflateDegree(a_in, N, opt.aFloor);
  if (finalCorr) {
    *finalCorr = 0.0;
  }
  if (n < 1) {
    return 0; // constant (or fully deflated): no isolated roots
  }

  std::complex<double> z[maxNq];
  aberthInit(n, opt.delta, z);
  aberthIterate(a_in, n, opt.niter, z, finalCorr);

  int m = 0;
  for (int i = 0; i < n; ++i) {
    if (std::fabs(z[i].imag()) >= opt.imagTol) {
      continue;
    }
    const double r = z[i].real();
    if (r < -1.0 - opt.margin || r > 1.0 + opt.margin) {
      continue;
    }
    double v, d;
    legendreValDer(a_in, n, r, v, d);
    roots[m] = r;
    slopes[m] = d;
    ++m;
  }

  // insertion sort: m <= N <= 20, and a sorted list lets the caller detect
  // multiple roots inside one sub-cell edge with a single pass
  for (int i = 1; i < m; ++i) {
    const double rk = roots[i], sk = slopes[i];
    int j = i - 1;
    while (j >= 0 && roots[j] > rk) {
      roots[j + 1] = roots[j];
      slopes[j + 1] = slopes[j];
      --j;
    }
    roots[j + 1] = rk;
    slopes[j + 1] = sk;
  }
  return m;
}

//  solveFibre on the REAL-PAIR kernel -- what v3_device.okl mirrors.
//
//  Identical deflation, filtering and sort to solveFibre above; only the
//  arithmetic differs (isoAberthReal.h, shared verbatim with the device).
//  Kept here, beside the complex version, so the two can be compared directly:
//  test_aberth_real.cpp requires 0 root-count mismatches over 20,400 fibres,
//  with roots agreeing to 2.8e-15.
inline int solveFibreReal(const double *a_in,
                          int N,
                          const AberthOpts &opt,
                          double *roots,
                          double *slopes,
                          double *finalCorr = nullptr)
{
  const int n = deflateDegree(a_in, N, opt.aFloor);
  if (finalCorr) {
    *finalCorr = 0.0;
  }
  if (n < 1) {
    return 0;
  }

  double zr[maxNq], zi[maxNq];
  isoAberthInitR(n, opt.delta, zr, zi);
  isoAberthIterateR(a_in, n, opt.niter, zr, zi, finalCorr);

  int m = 0;
  for (int i = 0; i < n; ++i) {
    if (std::fabs(zi[i]) >= opt.imagTol) {
      continue;
    }
    const double r = zr[i];
    if (r < -1.0 - opt.margin || r > 1.0 + opt.margin) {
      continue;
    }
    double v, d;
    legendreValDer(a_in, n, r, v, d);
    roots[m] = r;
    slopes[m] = d;
    ++m;
  }
  for (int i = 1; i < m; ++i) {
    const double rk = roots[i], sk = slopes[i];
    int j = i - 1;
    while (j >= 0 && roots[j] > rk) {
      roots[j + 1] = roots[j];
      slopes[j + 1] = slopes[j];
      --j;
    }
    roots[j + 1] = rk;
    slopes[j + 1] = sk;
  }
  return m;
}

//  Barycentric Lagrange evaluation of nodal data at an arbitrary point.
//  Used for (a) the transverse derivative, (b) the isoparametric map-out.
//  O(N).  bw are the barycentric weights of the GLL nodes; any common scale
//  cancels, so the closed form (-1)^j sqrt(w_j) may be used directly.
inline double baryEval(const double *fn, const double *xn, const double *bw, int Nq, double t)
{
  double num = 0.0, den = 0.0;
  for (int j = 0; j < Nq; ++j) {
    const double d = t - xn[j];
    if (std::fabs(d) < 1e-14) {
      return fn[j]; // exactly on a node
    }
    const double w = bw[j] / d;
    num += w * fn[j];
    den += w;
  }
  return num / den;
}

// Barycentric weights L_j(t) themselves (they sum to 1), for contracting a whole
// slice of nodal data in one pass.
inline void baryWeightsAt(const double *xn, const double *bw, int Nq, double t, double *L)
{
  for (int j = 0; j < Nq; ++j) {
    if (std::fabs(t - xn[j]) < 1e-14) {
      for (int i = 0; i < Nq; ++i) {
        L[i] = 0.0;
      }
      L[j] = 1.0;
      return;
    }
  }
  double den = 0.0;
  for (int j = 0; j < Nq; ++j) {
    L[j] = bw[j] / (t - xn[j]);
    den += L[j];
  }
  for (int j = 0; j < Nq; ++j) {
    L[j] /= den;
  }
}

//  Nodal -> modal (discrete Legendre projection), applied to one fibre.
//
//      a_k = (1/gamma_k) sum_j w_j f_j P_k(x_j)
//
//  NOTE gamma_N = 2/N, not 2/(2N+1): GLL quadrature is exact only to degree
//  2N-1, and P_N^2 has degree 2N.  Getting this wrong corrupts a_N -- precisely
//  the coefficient deflation inspects.  Report section 1.3.
//  `M` is the precomputed (N+1)x(N+1) matrix, row-major, built once at setup.
inline void nodalToModal(const double *M, const double *f, int Nq, double *a)
{
  for (int k = 0; k < Nq; ++k) {
    double s = 0.0;
    const double *Mk = M + static_cast<std::size_t>(k) * Nq;
    for (int j = 0; j < Nq; ++j) {
      s += Mk[j] * f[j];
    }
    a[k] = s;
  }
}

//  Derivative of the barycentric Lagrange basis, L_j'(t), EXACT.
//
//  With u_j = w_j/(t-x_j) and S = sum_k u_k we have L_j = u_j/S, so
//      u_j' = -u_j/(t-x_j),   L_j' = u_j'/S - L_j (S'/S).
//  O(N) and exact -- unlike the 1e-6 central difference used for the geometry
//  Jacobian, which is fine there (the rotation only needs geometry accuracy) but
//  would cap the Newton projection at sqrt(eps).  At a node the expression is
//  0/0; the standard node formula is used instead.
inline void baryDerivAt(const double *xn, const double *bw, int Nq, double t, double *dL)
{
  for (int m = 0; m < Nq; ++m) {
    if (std::fabs(t - xn[m]) < 1e-14) { // exactly on node m
      double acc = 0.0;
      for (int j = 0; j < Nq; ++j) {
        if (j == m) {
          continue;
        }
        dL[j] = (bw[j] / bw[m]) / (xn[m] - xn[j]);
        acc += dL[j];
      }
      dL[m] = -acc;
      return;
    }
  }
  double u[maxNq], S = 0.0, Sp = 0.0;
  for (int j = 0; j < Nq; ++j) {
    u[j] = bw[j] / (t - xn[j]);
    S += u[j];
    Sp += -u[j] / (t - xn[j]);
  }
  for (int j = 0; j < Nq; ++j) {
    const double up = -u[j] / (t - xn[j]);
    dL[j] = up / S - (u[j] / S) * (Sp / S);
  }
}

//  Value and reference gradient of a tensor-product nodal field at an arbitrary
//  (r,s,t).  One pass over the (N+1)^3 nodes, four accumulators.  `u` is the
//  element's nodal data in NekRS order: u[(ic*Nq + ib)*Nq + ia].
inline double tensorEvalGrad3(const double *u,
                              int Nq,
                              const double *Lr,
                              const double *Ls,
                              const double *Lt,
                              const double *Dr,
                              const double *Ds,
                              const double *Dt,
                              double *g)
{
  double f = 0.0, g0 = 0.0, g1 = 0.0, g2 = 0.0;
  for (int c = 0; c < Nq; ++c) {
    for (int b = 0; b < Nq; ++b) {
      const int row = (c * Nq + b) * Nq;
      const double LsLt = Ls[b] * Lt[c];
      const double DsLt = Ds[b] * Lt[c];
      const double LsDt = Ls[b] * Dt[c];
      for (int a = 0; a < Nq; ++a) {
        const double v = u[row + a];
        f += Lr[a] * LsLt * v;
        g0 += Dr[a] * LsLt * v;
        g1 += Lr[a] * DsLt * v;
        g2 += Lr[a] * LsDt * v;
      }
    }
  }
  if (g) {
    g[0] = g0;
    g[1] = g1;
    g[2] = g2;
  }
  return f;
}

//  Value of a tensor-product nodal field at one point -- no gradient.
//
//  Ordering as tensorEvalGrad3, u[(ic*Nq + ib)*Nq + ia], one accumulator not
//  four.  The weights Lr/Ls/Lt are field-independent: a caller sampling several
//  fields at one point computes them ONCE and calls this per component.
inline double tensorEval3(const double *u, int Nq, const double *Lr, const double *Ls, const double *Lt)
{
  double f = 0.0;
  for (int c = 0; c < Nq; ++c) {
    for (int b = 0; b < Nq; ++b) {
      const int row = (c * Nq + b) * Nq;
      const double w = Ls[b] * Lt[c];
      double inner = 0.0;
      for (int a = 0; a < Nq; ++a) {
        inner += Lr[a] * u[row + a];
      }
      f += w * inner;
    }
  }
  return f;
}

// --------------------------------------------------------------------------- //
//  Marching-cubes topology, GENERATED rather than tabulated by hand.
//
//  Corner c in [0,8) has (i,j,k) = (c&1, (c>>1)&1, (c>>2)&1); the 12 edges are
//  the corner pairs differing in one bit, sorted so the numbering is
//  reproducible.
//
//  Only the EDGE table is needed.  The 256x16 triangle table and MC33's ~730
//  subcases are not: cut points are connected by the ring walk, and face
//  ambiguity is settled by EVALUATING the interpolant (DEVELOPER.md).
// --------------------------------------------------------------------------- //
struct CubeEdge {
  int a, b, dir; // lo corner, hi corner, axis the edge runs along
};

inline const CubeEdge *cubeEdges()
{
  static CubeEdge E[12];
  static bool built = false;
  if (!built) {
    int n = 0;
    // the sorted order of (a,b) pairs; a < b always, since b = a | (1<<d)
    for (int a = 0; a < 8; ++a) {
      for (int d = 0; d < 3; ++d) {
        const int b = a | (1 << d);
        if (b != a) {
          E[n].a = a;
          E[n].b = b;
          E[n].dir = d;
          ++n;
        }
      }
    }
    built = true;
  }
  return E;
}

// bit e of mcEdgeTable()[idx] is set when edge e has endpoints of opposite sign,
// where bit c of idx is set when corner c has f < 0.
inline const int *mcEdgeTable()
{
  static int T[256];
  static bool built = false;
  if (!built) {
    const CubeEdge *E = cubeEdges();
    for (int idx = 0; idx < 256; ++idx) {
      int m = 0;
      for (int e = 0; e < 12; ++e) {
        if (((idx >> E[e].a) & 1) != ((idx >> E[e].b) & 1)) {
          m |= (1 << e);
        }
      }
      T[idx] = m;
    }
    built = true;
  }
  return T;
}

// Check the generated table against its DEFINING property rather than against a
// reference table: an all-same configuration cuts nothing, an isolated corner
// cuts exactly its three edges, and every bit agrees with the sign test.
inline bool verifyMcTables()
{
  const CubeEdge *E = cubeEdges();
  const int *T = mcEdgeTable();

  // the 12 edges are distinct, each joins corners differing in one bit
  int seen = 0;
  for (int e = 0; e < 12; ++e) {
    const int x = E[e].a ^ E[e].b;
    if (x != (1 << E[e].dir) || E[e].a >= E[e].b) {
      return false;
    }
    const int key = E[e].a * 8 + E[e].b;
    if (seen & (1 << (e))) {
      return false;
    }
    seen |= (1 << e);
    (void)key;
  }
  for (int e = 0; e < 12; ++e) {
    for (int g = e + 1; g < 12; ++g) {
      if (E[e].a == E[g].a && E[e].b == E[g].b) {
        return false;
      }
    }
  }

  if (T[0] != 0 || T[255] != 0) {
    return false;
  }
  for (int idx = 0; idx < 256; ++idx) {
    for (int e = 0; e < 12; ++e) {
      const bool cut = (((idx >> E[e].a) & 1) != ((idx >> E[e].b) & 1));
      if (((T[idx] >> e) & 1) != (cut ? 1 : 0)) {
        return false;
      }
    }
  }
  for (int c = 0; c < 8; ++c) { // an isolated corner cuts exactly 3 edges
    int nb = 0;
    for (int e = 0; e < 12; ++e) {
      if ((T[1 << c] >> e) & 1) {
        ++nb;
      }
    }
    if (nb != 3) {
      return false;
    }
  }
  for (int idx = 0; idx < 256; ++idx) { // complementary configs cut the same edges
    if (T[idx] != T[255 - idx]) {
      return false;
    }
  }
  return true;
}

// --------------------------------------------------------------------------- //
//  Newton projection onto the level set, FIXED iteration count.
//
//      (r,s,t) <- (r,s,t) - f * grad_r f / |grad_r f|^2
//
//  In the REFERENCE frame, so no inverse geometry map is needed.  No early exit:
//  identical work per point, so a warp never diverges -- a degenerate gradient
//  ZEROES the step rather than skipping it.  Returns true if all steps ran.
//
//  `lock` is a 3-bit mask of directions the step may not move.  A node whose two
//  parents sit on a shared element face must STAY on that face: constrained, the
//  step is driven by the field restricted to the face, which is identical from
//  both sides, so the neighbours agree to round-off and the edge closes.
//  Unconstrained they drift apart -- a visible crack at P2.  Clamping to
//  [-1,1]^3 is separate and does not fix it; the drifting node moves INWARD.
// --------------------------------------------------------------------------- //
inline bool newtonProject(const double *u,
                          int Nq,
                          const double *xn,
                          const double *bw,
                          int niter,
                          double gradRelTol,
                          double fScale,
                          double *rst,
                          int lock = 0)
{
  double Lr[maxNq], Ls[maxNq], Lt[maxNq], Dr[maxNq], Ds[maxNq], Dt[maxNq], g[3];
  bool full = true;
  const double gTol = gradRelTol * fScale;

  for (int it = 0; it < niter; ++it) {
    baryWeightsAt(xn, bw, Nq, rst[0], Lr);
    baryWeightsAt(xn, bw, Nq, rst[1], Ls);
    baryWeightsAt(xn, bw, Nq, rst[2], Lt);
    baryDerivAt(xn, bw, Nq, rst[0], Dr);
    baryDerivAt(xn, bw, Nq, rst[1], Ds);
    baryDerivAt(xn, bw, Nq, rst[2], Dt);

    const double f = tensorEvalGrad3(u, Nq, Lr, Ls, Lt, Dr, Ds, Dt, g);

    // Restrict the gradient to the FREE subspace.  Dropping the locked
    // components (rather than stepping and then re-imposing the constraint) is
    // what makes this a Newton step for the restricted problem: it is the exact
    // gradient of f viewed as a function on the face, so the quadratic rate is
    // kept.  Zeroing after the step would instead be a projected step, which
    // converges only linearly.
    for (int d = 0; d < 3; ++d) {
      if (lock & (1 << d)) {
        g[d] = 0.0;
      }
    }
    const double g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
    const double gn = std::sqrt(g2);

    // mask, not a branch out of the loop: a degenerate gradient gives a zero step
    const double scale = (gn > gTol) ? (f / g2) : 0.0;
    // A FULLY LOCKED point (lock == 0b111) has no free direction at all, so its
    // restricted gradient is zero BY CONSTRUCTION -- it sits at a cube corner,
    // where all three coordinates are pinned.  There is nothing to solve and
    // nothing wrong: reporting it as a degenerate gradient would conflate "the
    // field went flat here" with "this node was never free to move", and the
    // first is a warning while the second is the constraint doing its job.
    if (gn <= gTol && lock != 0x7) {
      full = false;
    }
    for (int d = 0; d < 3; ++d) {
      double v = rst[d] - scale * g[d]; // locked d: g[d] == 0, so v == rst[d]
      v = (v < -1.0) ? -1.0 : (v > 1.0 ? 1.0 : v);
      rst[d] = v;
    }
  }
  return full;
}

} // namespace isoDetail

// =========================================================================== //
//  INLINED from plugin/v1.hpp -- the data layer the writers and the caller see.
// =========================================================================== //
//  Parameters
struct IsoSurfaceParams {
  std::vector<dfloat> levels; // isolevels; typically <= 10

  // ---- element culling -----------------------------------------------------
  // Skip an element when the level lies outside
  //     [min_j f_j - theta*K_N*R,  max_j f_j + theta*K_N*R],  R = range of f_j.
  // theta = 1 is rigorous and the default; 0 is the naive nodal min/max, which
  // is NOT (an interpolant leaves its nodal range).  README.md has the trade.
  bool cullElements = true;
  dfloat cullTheta = 1.0;

  // ---- route C -------------------------------------------------------------
  dfloat aFloor = 1e-12; // degree deflation (report 3.3)
  int niter = 12;        // FIXED Aberth iterations (report 5.3)
  dfloat delta = 0.1;    // initial-guess offset off the real axis (notes/14)
  dfloat imagTol = 1e-7; // |Im z| < imagTol counts as real
  dfloat margin = 1e-10; // keep |r| <= 1+margin; roots beyond are other elements

  // ---- output --------------------------------------------------------------
  // Tangency detection, RELATIVE because a computed double root has |p'| at
  // ~sqrt(eps)*scale, not 0 -- an absolute floor misses it in a large field and
  // cries wolf in a small one:
  //     tangency  <=>  |p'(root)| < gradRelTol * max_j |f_j - c|
  dfloat gradRelTol = 1e-4;

  // ---- what gets written ---------------------------------------------------
  // THE CALLER DECIDES, by choosing a verb -- there are no dump flags:
  //     process()      -> <stem>_fld_<LL>_<SSSSS>.vtu    the deliverable
  //     writeSurface() -> <stem>_srf_<LL>_<SSSSS>.vtu    debug: geometry only
  //     writeRoots()   -> <stem>_roots_<LL>_<SSSSS>.vtu  debug: the raw crossings
  // Each verb writes exactly its own file and nothing else; a caller wanting two
  // of them calls two verbs.  The roots and the triangulation are intermediate
  // stages that exist so the fields file can be built and debugged.

  // Default for the 2-arg process(); process(t, step, solveGeom) overrides it.
  // FALSE reuses the surface built by the first process() -- right for a frozen
  // scalar, wrong for a changing one, and the plugin CANNOT tell the difference
  // (it holds the caller's array by handle).  A moving mesh is not detected.
  bool ifSolveGeom = false;

  // Multiplies the per-element device capacity in the v2/OKL path (plan E):
  // CAP_R = capScale * max(3(N+1)^2, (N+1)^3).  1.0 is measured-safe -- the
  // worst real case (poly5, every fibre carrying its full N roots) peaks at 83%
  // -- and the kernel reports an overflow rather than truncating.  It is a
  // RUNTIME kernel argument, so raising it needs no rebuild.  Ignored by v1.
  dfloat capScale = 1.0;

  // Precision of the interpolated field values in the .vtu.  Float32 halves the
  // file and is what a visualiser consumes; set true if the file feeds a solver.
  // (Geometry and normals are always Float32, as in the other two writers.)
  bool fieldsFloat64 = false;

  // ---- surface (milestones 2 and 3) ----------------------------------------
  // cellOrder = 1  flat triangles,      VTK_TRIANGLE (5),            3 nodes
  // cellOrder = 2  curved triangles,    VTK_QUADRATIC_TRIANGLE (22), 6 nodes
  //
  // Order 2 adds ONE node per triangulation edge, each Newton-projected onto the
  // level set.  Corners are shared, so the measured cost is 2.76 new points per
  // root (V:E:F = 1:2.76:1.77 on real extracted surfaces), and the in-triangle
  // error drops ~100x: 3.7e-2 flat -> 3.6e-4 curved on a 22.5-degree patch.
  // Order is RECORDED IN THE FILE (FieldData "cellOrder"), so a reader never has
  // to infer it from the cell types.  See plans/C_topology_and_curved_surface.md.
  int cellOrder = 1;

  // Newton projection of the midside nodes: FIXED sweeps, no convergence test,
  // one batched pass over every node this rank owns.  The start point is a chord
  // midpoint between two exact roots, i.e. already O(h^2) out, so convergence is
  // quadratic from the first step -- measured 4.4e-2 -> 1.7e-3 -> 2.5e-6 ->
  // 5.0e-12, with the quadratic constant 0.83 = 1/(2R) as theory says.  Three
  // sweeps reach the floor; a fourth is free if margin is wanted.
  int niterProj = 3;

  // Output stem.  Files are
  //     <fileName>_roots_<LL>_<SSSSS>.vtu    the crossings   (writeRoots)
  //     <fileName>_srf_<LL>_<SSSSS>.vtu      the surface     (writeSurface)
  //     <fileName>_fld_<LL>_<SSSSS>.vtu      surface+fields  (process)
  // LL is the isolevel index; SSSSS is a DENSE dump counter, not the timestep,
  // which is what keeps ParaView's file-series detection happy.  time and tstep
  // are inside as FieldData.
  std::string fileName = "iso";
};

//  One level's worth of roots.  Structure-of-arrays: this is what gets written
//  and what an OCCA port would fill directly.
struct RootSet {
  dfloat level = 0;
  dlong nRoots = 0;

  std::vector<dlong> elem;              // LOCAL element id (this rank)
  std::vector<hlong> eg;                // GLOBAL element id (unique across ranks)
  std::vector<uint8_t> dir;             // fibre direction: 0 = r, 1 = s, 2 = t
  std::vector<int16_t> j, k;            // the fibre's other two node indices
  std::vector<dfloat> r, s, t;          // REFERENCE coordinates (exact)
  std::vector<dfloat> x, y, z;          // PHYSICAL coordinates
  std::vector<dfloat> dfdr, dfds, dfdt; // gradient, reference frame
  std::vector<dfloat> dfdx, dfdy, dfdz; // gradient, physical frame
  std::vector<uint8_t> flag;            // see FlagBits

  // Per-element lookup, so a consumer can address one element's roots without
  // scanning the whole list:
  //     roots of LOCAL element e are the entries [eStart[e], eStart[e+1])
  // eCount[e] is the same information as a count.  Both are LOCAL and sized
  // Nelements (+1 for eStart), i.e. they cover culled elements too -- a culled
  // element simply has eCount == 0, which is exactly what a reader needs to
  // know.  Written to the .vtu as FieldData so the table survives the file.
  std::vector<dlong> eStart; // [Nelements+1], prefix sum of eCount
  std::vector<dlong> eCount; // [Nelements]

  enum FlagBits : uint8_t {
    Tangency = 1u << 0,  // |grad f| < gradFloor: sqrt(eps) accuracy, do not trust
    MultiEdge = 1u << 1, // >= 2 roots inside one sub-cell edge (report 6.3)
    OnFace = 1u << 2     // within margin of |r| = 1: the neighbour sees it too
  };

  void clear()
  {
    nRoots = 0;
    eStart.clear();
    eCount.clear();
    for (auto *v : {&elem}) {
      v->clear();
    }
    eg.clear();
    dir.clear();
    j.clear();
    k.clear();
    r.clear();
    s.clear();
    t.clear();
    x.clear();
    y.clear();
    z.clear();
    dfdr.clear();
    dfds.clear();
    dfdt.clear();
    dfdx.clear();
    dfdy.clear();
    dfdz.clear();
    flag.clear();
  }

  void reserve(dlong n)
  {
    elem.reserve(n);
    eg.reserve(n);
    dir.reserve(n);
    j.reserve(n);
    k.reserve(n);
    r.reserve(n);
    s.reserve(n);
    t.reserve(n);
    x.reserve(n);
    y.reserve(n);
    z.reserve(n);
    dfdr.reserve(n);
    dfds.reserve(n);
    dfdt.reserve(n);
    dfdx.reserve(n);
    dfdy.reserve(n);
    dfdz.reserve(n);
    flag.reserve(n);
  }
};

// --------------------------------------------------------------------------- //
//  One level's worth of triangulated surface -- the milestone 2/3 product.
//
//  Points come in two groups, in this order:
//      [0, nCorner)      CORNERS -- these ARE roots: corner i is root
//                        cornerRoot[i], and marching only decides which to
//                        connect.  A cut point on an inserted lattice plane
//                        follows no fibre and carries cornerRoot = -1.
//      [nCorner, nPoint) MIDSIDE nodes (cellOrder == 2), one per unique
//                        triangulation edge, Newton-projected onto the level set.
//
//  Structure-of-arrays, like RootSet: the writer packs SoA straight into the
//  appended VTK blocks.
// --------------------------------------------------------------------------- //
struct SurfaceSet {
  dfloat level = 0;
  uint8_t cellType = 5; // 5 = VTK_TRIANGLE, 22 = VTK_QUADRATIC_TRIANGLE
  int nodesPerCell = 3; // 3 or 6
  int cellOrder = 1;    // 1 or 2; written to the file as FieldData

  dlong nCorner = 0, nPoint = 0;
  std::vector<dfloat> px, py, pz; // physical
  std::vector<dfloat> pr, ps, pt; // REFERENCE, so a consumer can re-enter the element
  std::vector<dfloat> nx, ny, nz; // unit normal, grad f / |grad f|
  std::vector<dlong> pelem;       // owning LOCAL element
  std::vector<dlong> cornerRoot;  // -> index into the RootSet, or -1 for an
                                  //    inserted-plane cut point (see Inserted)
  std::vector<uint8_t> pflag;     // root flags; midside nodes inherit 0

  dlong nCell = 0;
  std::vector<dlong> conn;   // [nCell * nodesPerCell], cell-major
  std::vector<hlong> cellEg; // [nCell] global element id
  std::vector<uint8_t> cellFlag;

  // Cell flags and point flags share these bits where they overlap, so a reader
  // can test them the same way.
  enum CellFlagBits : uint8_t {
    Ambiguous = 1u << 0, // >4 cut points: fan resolved by EVALUATING the interpolant
    Inserted = 1u << 1,  // in a sub-cell created by lattice-plane insertion
    Tangency = 1u << 2   // a midside node hit a degenerate gradient, step zeroed
  };

  void clear()
  {
    nCorner = nPoint = nCell = 0;
    px.clear();
    py.clear();
    pz.clear();
    pr.clear();
    ps.clear();
    pt.clear();
    pelem.clear();
    nx.clear();
    ny.clear();
    nz.clear();
    cornerRoot.clear();
    pflag.clear();
    conn.clear();
    cellEg.clear();
    cellFlag.clear();
  }
};

// --------------------------------------------------------------------------- //
//  One level's worth of INTERPOLATED field values -- the milestone 5 product.
//
//  Values are sampled at every surface point, corners AND the cellOrder-2 midside
//  nodes.  A midside node is a real point on the surface carrying real physics;
//  filling it by averaging its two parents would put P1 data on a P2 cell, so the
//  geometry would be quadratic while the field stayed linear -- a curved triangle
//  rendered with a flat colour ramp.
//
//  val[f] is point-major, [nComp[f] * nPoint], so one component group is one
//  contiguous appended VTK block, matching how the other writers pack.
// --------------------------------------------------------------------------- //
struct FieldSet {
  std::vector<std::string> name; // [nField]
  std::vector<int> nComp;        // [nField]
  std::vector<std::vector<dfloat>> val;
  dlong nPoint = 0;

  int nField() const
  {
    return static_cast<int>(name.size());
  }

  void clear()
  {
    name.clear();
    nComp.clear();
    val.clear();
    nPoint = 0;
  }
};

struct IsoStats {
  dlong nElements = 0;       // local elements
  dlong nElementsSolved = 0; // survived culling (summed over levels)
  dlong nFibres = 0;         // fibre solves actually performed
  dlong nRoots = 0;          // roots kept, all levels
  dlong nTangency = 0;
  dlong nMultiEdge = 0;
  dlong nOnFace = 0;
  double maxFinalCorr = 0; // largest Aberth final correction: the safety net

  // milestone 2/3
  dlong nCell = 0;           // triangles, all levels
  dlong nCellAmbiguous = 0;  // fans with more than 4 cut points
  dlong nMidNode = 0;        // midside nodes created (cellOrder == 2)
  dlong nProjDegenerate = 0; // midside nodes whose Newton step was zeroed
  double maxProjResid = 0;   // largest |f| left at a midside node, relative

  // milestone 5
  dlong nInterpPoint = 0; // surface points sampled (summed over fields)
  dlong nSolveCalls = 0;  // times the root solve actually ran: the ifSolveGeom check

  // ---- per-element occupancy (notes/25) -----------------------------------
  // How close the worst element came to the capacity, which is the number that
  // says whether a run is comfortable or one fold away from losing data.  Two
  // aborts reached a user because the only thing being reported was a total:
  // a sum says nothing about the element that overflows.
  //
  // "Crossed" excludes elements with zero roots -- on a culled mesh most
  // elements are empty and they drag the average to meaninglessness.
  dlong nElemCrossed = 0; // elements with >= 1 root (summed over levels)
  dlong maxRootPerElem = 0;
  dlong minRootPerElem = 0; // over CROSSED elements only
  dlong maxCellPerElem = 0;
  dlong minCellPerElem = 0; // over elements with >= 1 triangle
  dlong nElemWithCell = 0;
};

// --------------------------------------------------------------------------- //
//  Per-element device capacity.  These are BOUNDS derived from the algorithm,
//  NOT measurements of what the test cases needed -- that distinction cost two
//  production aborts (notes/24, notes/25).  Do not lower them.
//
//      capRoot = capScale * 3N(N+1)^2   3(N+1)^2 fibres x at most N roots each
//      capCell = p_isoMaxTri * N^3      each sub-cell owns a FIXED slice, so a
//                                       smaller value writes past the buffer
//
//  p_isoMaxTri = p_isoMaxCut - 2 = 22, not the textbook marching-cubes 5: this
//  scheme has no triangle table, and a ring of np cut points emits np-2
//  triangles with np <= 24 (12 edges, each able to carry two crossings).
//  DEVELOPER.md derives all three.
// --------------------------------------------------------------------------- //
struct IsoCaps {
  dlong root = 0;
  dlong cell = 0;
  dlong point = 0;

  static IsoCaps forOrder(int N, double scale)
  {
    const dlong Nq = N + 1;
    const dlong fib = 3 * Nq * Nq;
    const dlong N3 = static_cast<dlong>(N) * N * N;
    IsoCaps c;
    // 3(N+1)^2 fibres x at most N roots each.  The max() with the fibre count
    // covers N = 0, where N*fib would be 0 but a constant field still needs a
    // slot per fibre to report "no root".
    const dlong bound = std::max(fib, static_cast<dlong>(N) * fib);
    c.root = static_cast<dlong>(scale * static_cast<double>(bound));
    // p_isoMaxTri * N^3, and p_isoMaxTri = p_isoMaxCut - 2 = 22: 24 cut points
    // (all 12 edges carrying a pair) in a SINGLE ring, since more rings means
    // FEWER triangles (np - 2R).
    //
    // It is pessimistic -- on the worst measured field only 83 of 4096 sub-cells
    // carry even one split edge, and the connectivity buffer goes 42 -> 93 MB at
    // N=7 on 512 elements.  It stays anyway.  The fixed-slice scheme gives each
    // sub-cell p_isoMaxTri slots unconditionally, so a smaller capCell is not a
    // tuning choice, it is a buffer overrun; and a capacity set from "what we
    // measured" rather than from the algorithm is exactly what produced the two
    // production aborts in notes/24 and notes/25.  Twice is enough.
    c.cell = 22 * N3;
    c.point = c.root + 3 * c.cell;
    return c;
  }
};

class isoSurface
{
public:
  // the fieldExtract convention, so a caller can reuse the same field list
  using field = std::tuple<std::string, std::vector<deviceMemory<dfloat>>>;

  isoSurface(mesh_t *mesh, const field &scalar, const IsoSurfaceParams &prm);
  isoSurface(mesh_t *mesh, const std::vector<dfloat> &scalar, const IsoSurfaceParams &prm);
  ~isoSurface() = default;

  void setField(const field &f);
  void setField(const std::vector<dfloat> &f);

  // Solve the fibre roots ON THE DEVICE: isoCullOKL -> isoSolveFibresOKL ->
  // isoCompactPointsOKL, then the host geometry the kernels do not produce.
  // Skipped when neither the field nor the mesh has changed.
  void solveRoots(bool force = false);

  bool rootsValid() const
  {
    return rootsValid_;
  }

  void markMeshChanged()
  {
    rootsValid_ = false;
    surfValid_ = false;
  }

  void setOutputFields(const std::vector<field> &flds);
  void setOutputFieldHost(const std::string &name, int nComp, const std::vector<dfloat> &vals);

  void buildTopology(bool force = false);
  void interpFields();

  void writeRoots(double time, int tstep);
  void writeSurface(double time, int tstep);
  void writeFields(double time, int tstep);
  // The deliverable: geometry (per prm.ifSolveGeom), fields sampled onto it,
  // one _fld_ file per level.  The scalar is re-read from the device each time
  // the geometry is solved, so a field changing in place needs no setField().
  void process(double time, int tstep);
  // Same, overriding prm.ifSolveGeom for this call: true re-solves the surface,
  // false reuses it.  This is the usual call in a time loop, where the scalar
  // changes every dump.
  void process(double time, int tstep, bool solveGeom);

  const RootSet &roots(int levelIdx) const
  {
    return rootSets_.at(levelIdx);
  }

  const SurfaceSet &surface(int levelIdx) const
  {
    return surfSets_.at(levelIdx);
  }

  const FieldSet &fields(int levelIdx) const
  {
    return fieldSets_.at(levelIdx);
  }

  const IsoStats &stats() const
  {
    return stats_;
  }

  const IsoSurfaceParams &params() const
  {
    return prm_;
  }

  const IsoCaps &caps() const
  {
    return caps_;
  }

private:
  void setupCommon(mesh_t *mesh);
  void allocDevice();
  void checkOverflow(const char *stage);

  enum Product { Roots = 0, Srf, Fld };

  std::string outName(int level, int counter, Product what) const;
  void writeVtm(Product what, int counter) const;
  void writeRootVtu(const RootSet &rs, int lev, const std::string &fname, double time, int tstep);
  void writeSurfVtu(const SurfaceSet &ss,
                    const FieldSet *fs,
                    int lev,
                    const std::string &fname,
                    double time,
                    int tstep);
  void ensureGeometry();
  void ensureGeometry(bool solveGeom);

  // Pull the scalar from the device into fld_.  Called at the top of every
  // actual solve, so the host geometry pass and the device upload both see the
  // caller's CURRENT values -- not whatever they held when setField() ran.
  void refreshFieldHost();

  // min/max/avg per element against a capacity, rank 0 prints.  Collective.
  void reportOccupancy(const char *what, dlong nCrossed, dlong locMin, dlong locMax, dlong locSum, dlong cap)
      const;

  // Residual distribution of the midside projection, so raising niterProj
  // is a decision made on data rather than a guess.  Collective; rank 0 prints.
  void reportProjection();

  // The COMPUTE the kernels do not cover -- see the file header.
  void gradAtRoot(dlong e, int dir, int j, int k, double rr, double slope, double *gref) const;
  void uploadRoots(const RootSet &rs);
  void uploadPoints(const SurfaceSet &ss);
  void pointNormal(dlong e, double r, double s, double t, const double *J, double *nrm) const;
  void mapToPhysical(dlong e, double rr, double ss, double tt, double *xyz, double *J) const;

  mesh_t *mesh_ = nullptr;
  IsoSurfaceParams prm_;
  IsoCaps caps_;
  std::string fieldName_;

  int N_ = 0, Nq_ = 0, Np_ = 0;
  std::vector<double> xn_, wn_, bw_, M_;
  double KN_ = 0;

  // The scalar's DEVICE array, kept as a handle rather than a snapshot: the
  // caller's o_field is re-read at the top of every solve, so a field that
  // changes in place needs no re-registration.  See refreshFieldHost().
  deviceMemory<dfloat> o_scalar_;
  bool haveScalarDev_ = false;

  // Host copies.  fld_ is still needed -- the host geometry pass reads nodal
  // values, and gx_/gy_/gz_ feed mapToPhysical.  It is a CACHE of o_scalar_ and
  // is refreshed from it by refreshFieldHost(), never trusted across a solve.
  std::vector<double> fld_;
  std::vector<double> gx_, gy_, gz_;

  std::vector<RootSet> rootSets_;
  std::vector<SurfaceSet> surfSets_;
  std::vector<FieldSet> fieldSets_;
  std::vector<field> outFields_;
  std::vector<std::string> hostName_;
  std::vector<int> hostComp_;
  std::vector<std::vector<double>> hostVal_;
  std::vector<std::vector<double>> outHost_;

  // ---- device state.  Allocated once in allocDevice(), reused every dump. ---
  deviceMemory<dfloat> o_xn_, o_bw_, o_M_, o_levels_, o_fld_;
  deviceMemory<int> o_elemCross_, o_overflow_;
  deviceMemory<int> o_rootCount_, o_rootDir_, o_rootJK_;
  deviceMemory<dfloat> o_rootR_, o_rootS_, o_rootT_, o_rootSlope_;
  deviceMemory<unsigned char> o_rootFlag_;

  // ---- step 2: topology ----------------------------------------------------
  deviceMemory<int> o_cellCount_, o_ptCount_, o_conn_, o_ptRoot_;
  deviceMemory<dfloat> o_ptR_, o_ptS_, o_ptT_;
  deviceMemory<unsigned char> o_ptFlag_, o_cellFlag_;
  deviceMemory<int> o_midCount_, o_midBase_, o_midLock_;
  // corners + midside per element.  isoInterpOKL loops [0, ptCount[e]) and
  // o_ptCount_ holds CORNERS ONLY, so passing it would silently skip every
  // midside node at cellOrder 2 and leave those field values at zero.
  deviceMemory<int> o_ptTotal_;
  deviceMemory<dfloat> o_src_, o_fldVal_;
  std::vector<int> ptTotal_;
  deviceMemory<int> o_nUnsep_, o_nDegen_;
  deviceMemory<int> o_projResidQ_, o_nProjShort_; // projection convergence report

  // compacted (dense) outputs of isoCompactPointsOKL
  deviceMemory<dlong> o_ptOffset_; // dlong, NOT int: the kernel says dlong*
  deviceMemory<int> o_outElem_, o_outRoot_;
  deviceMemory<dfloat> o_outR_, o_outS_, o_outT_;
  deviceMemory<unsigned char> o_outFlag_;

  bool surfValid_ = false;
  IsoStats stats_;
  bool printedSetup_ = false;
  bool rootsValid_ = false;
  int dumpCounter_ = 0;
  hlong egOffset_ = 0;
};

inline void isoSurface::setupCommon(mesh_t *mesh)
{
  const double tStart = MPI_Wtime();

  mesh_ = mesh;
  N_ = mesh->N;
  Nq_ = mesh->Nq;
  Np_ = mesh->Np;

  nekrsCheck(N_ < 1 || N_ > isoDetail::maxN,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "isoSurface: polynomial order %d outside the supported range [1,%d]\n",
             N_,
             isoDetail::maxN);

  // ---- 1D GLL nodes and weights (NekRS's own routine) ----------------------
  {
    std::vector<dfloat> xt(Nq_), wt(Nq_);
    JacobiGLL(N_, xt.data(), wt.data());
    xn_.resize(Nq_);
    wn_.resize(Nq_);
    for (int i = 0; i < Nq_; ++i) {
      xn_[i] = static_cast<double>(xt[i]);
      wn_[i] = static_cast<double>(wt[i]);
    }
  }

  // ---- barycentric weights.  The closed form b_j = (-1)^j sqrt(w_j) is only
  //      PROPORTIONAL to the true weights, but every formula here is homogeneous
  //      of degree 0 in b, so the common factor cancels (report section 1.3).
  bw_.resize(Nq_);
  for (int j = 0; j < Nq_; ++j) {
    bw_[j] = ((j % 2) ? -1.0 : 1.0) * std::sqrt(wn_[j]);
  }

  // ---- nodal -> modal.  gamma_N = 2/N, NOT 2/(2N+1): GLL quadrature is exact
  //      only to degree 2N-1 and P_N^2 has degree 2N.  Getting this wrong
  //      corrupts a_N, which is precisely what deflation inspects.
  M_.assign(static_cast<size_t>(Nq_) * Nq_, 0.0);
  for (int k = 0; k < Nq_; ++k) {
    const double gam = (k < N_) ? 2.0 / (2 * k + 1) : 2.0 / N_;
    for (int j = 0; j < Nq_; ++j) {
      double P0 = 1.0, P1 = xn_[j], Pk = (k == 0) ? 1.0 : xn_[j];
      for (int m = 1; m < k; ++m) {
        const double P2 = ((2 * m + 1) * xn_[j] * P1 - m * P0) / (m + 1);
        P0 = P1;
        P1 = P2;
        Pk = P2;
      }
      M_[static_cast<size_t>(k) * Nq_ + j] = wn_[j] * Pk / gam;
    }
  }

  KN_ = isoDetail::gllOvershoot(N_);

  // ---- geometry, once: the map-out needs it on the host --------------------
  {
    auto [hx, hy, hz] = mesh_->xyzHost();
    gx_.resize(mesh_->Nlocal);
    gy_.resize(mesh_->Nlocal);
    gz_.resize(mesh_->Nlocal);
    for (dlong i = 0; i < mesh_->Nlocal; ++i) {
      gx_[i] = static_cast<double>(hx[i]);
      gy_[i] = static_cast<double>(hy[i]);
      gz_[i] = static_cast<double>(hz[i]);
    }
  }

  // Global element numbering: NekRS carries node-level globalIds, not an
  // element-level one, so derive it the standard way.  eg = egOffset_ + e is
  // then unique across ranks and stable for a given decomposition.
  {
    hlong nel = static_cast<hlong>(mesh_->Nelements);
    egOffset_ = 0;
    MPI_Exscan(&nel, &egOffset_, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
    if (platform->comm.mpiRank() == 0) {
      egOffset_ = 0;
    }
  }

  fld_.assign(mesh_->Nlocal, 0.0);
  rootSets_.resize(prm_.levels.size());
  surfSets_.resize(prm_.levels.size());
  fieldSets_.resize(prm_.levels.size());
  for (size_t l = 0; l < prm_.levels.size(); ++l) {
    rootSets_[l].level = prm_.levels[l];
  }

  platform->timer.addUserStat("isoSurface::");
  platform->timer.set("isoSurface::setup", MPI_Wtime() - tStart);

  // ---- v3: size the per-element device capacity and allocate ---------------
  caps_ = IsoCaps::forOrder(N_, static_cast<double>(prm_.capScale));
  allocDevice();

  // ONE setup line.  The tunables that matter are theta (notes/15) and the
  // capacity; niter/aFloor/margin are fixed defaults nobody changes, and the
  // fibre count is implied by N.  Anything more belongs in --verbose.
  // printedSetup_ is per-object, so with several samplers in one .udf the banner
  // repeats verbatim.  Key it on the CONFIGURATION instead: identical settings
  // print once, a genuinely different one still announces itself.
  if (platform->comm.mpiRank() == 0 && !printedSetup_) {
    printedSetup_ = true;
    static std::string lastBanner;
    char cfg[256];
    snprintf(cfg,
             sizeof(cfg),
             "N=%d, %d level(s), theta=%.3g%s, cap %d/%d/%d (root/cell/point)",
             N_,
             static_cast<int>(prm_.levels.size()),
             static_cast<double>(prm_.cullTheta),
             prm_.cullElements ? "" : " (culling OFF)",
             static_cast<int>(caps_.root),
             static_cast<int>(caps_.cell),
             static_cast<int>(caps_.point));
    if (lastBanner != cfg) {
      lastBanner = cfg;
      printf("isoSurface[GPU]: %s\n", cfg);
    }
    if (platform->verbose()) {
      printf("isoSurface(%s): %d fibres/element, K_N=%.6f, niter=%d, aFloor=%.1e, margin=%.1e, "
             "niterProj=%d\n",
             fieldName_.c_str(),
             3 * Nq_ * Nq_,
             KN_,
             prm_.niter,
             static_cast<double>(prm_.aFloor),
             static_cast<double>(prm_.margin),
             prm_.niterProj);
    }
  }
}

inline isoSurface::isoSurface(mesh_t *mesh, const field &scalar, const IsoSurfaceParams &prm) : prm_(prm)
{
  setupCommon(mesh);
  setField(scalar);
}

inline isoSurface::isoSurface(mesh_t *mesh, const std::vector<dfloat> &scalar, const IsoSurfaceParams &prm)
    : prm_(prm)
{
  fieldName_ = "scalar";
  setupCommon(mesh);
  setField(scalar);
}

//  field ingestion
// =========================================================================== //
inline void isoSurface::setField(const field &f)
{
  const auto &comps = std::get<1>(f);
  nekrsCheck(comps.size() != 1,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "isoSurface: expected a SCALAR field, got %d components\n",
             static_cast<int>(comps.size()));

  fieldName_ = std::get<0>(f);
  rootsValid_ = false; // the scalar changed: the roots are stale
  surfValid_ = false;

  // Keep the HANDLE, not a snapshot.  Every solve re-reads through it, so a
  // caller whose field changes in place calls this once, in setup.
  o_scalar_ = comps[0];
  haveScalarDev_ = true;
  refreshFieldHost();
}

// Pull the device scalar into the host cache.  This is what makes a stored
// handle safe: fld_ (the host geometry pass) and o_fld_ (uploaded from fld_ in
// solveRoots) would otherwise both hold whatever the field contained when
// setField() ran, and a time loop would re-solve a stale surface every dump.
//
// A host-array caller (the other setField overload) has no device handle; its
// values are already in fld_ and this is a no-op.
inline void isoSurface::refreshFieldHost()
{
  if (!haveScalarDev_) {
    return;
  }
  // deviceMemory::copyTo takes a std::vector, not a raw pointer (the pointer
  // overloads are #if 0'd out in deviceMemory.hpp), and the size argument is a
  // COUNT of elements.
  std::vector<dfloat> tmp(mesh_->Nlocal);
  o_scalar_.copyTo(tmp, static_cast<size_t>(mesh_->Nlocal));
  for (dlong i = 0; i < mesh_->Nlocal; ++i) {
    fld_[i] = static_cast<double>(tmp[i]);
  }
}

inline void isoSurface::setField(const std::vector<dfloat> &f)
{
  nekrsCheck(static_cast<dlong>(f.size()) < mesh_->Nlocal,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "isoSurface: host field has %d entries, need %d\n",
             static_cast<int>(f.size()),
             static_cast<int>(mesh_->Nlocal));
  rootsValid_ = false; // the scalar changed: the roots are stale
  surfValid_ = false;
  haveScalarDev_ = false; // a host field supersedes any stored device handle
  for (dlong i = 0; i < mesh_->Nlocal; ++i) {
    fld_[i] = static_cast<double>(f[i]);
  }
}

// =========================================================================== //
//  geometry: (r,s,t) -> (x,y,z) and the Jacobian, by tensor barycentric
//  interpolation of the element's own geometry nodes.  The map is isoparametric,
//  so a root exact in (r,s,t) is exact in (x,y,z) to the accuracy of the
//  geometry representation -- no error is added (report section 6.1).
// =========================================================================== //
inline void isoSurface::mapToPhysical(dlong e, double rr, double ss, double tt, double *xyz, double *J) const
{
  double Lr[isoDetail::maxNq], Ls[isoDetail::maxNq], Lt[isoDetail::maxNq];
  double Dr[isoDetail::maxNq], Ds[isoDetail::maxNq], Dt[isoDetail::maxNq];

  isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, rr, Lr);
  isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, ss, Ls);
  isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, tt, Lt);

  if (J) {
    // dL_i/dt by a central difference on the barycentric weights: the map is a
    // smooth polynomial and this is only used for the gradient rotation, whose
    // accuracy needs to match the geometry, not machine epsilon.
    const double h = 1e-6;
    double Lp[isoDetail::maxNq], Lm[isoDetail::maxNq];
    auto deriv = [&](double c, double *out) {
      isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, std::min(c + h, 1.0), Lp);
      isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, std::max(c - h, -1.0), Lm);
      const double dd = std::min(c + h, 1.0) - std::max(c - h, -1.0);
      for (int i = 0; i < Nq_; ++i) {
        out[i] = (Lp[i] - Lm[i]) / dd;
      }
    };
    deriv(rr, Dr);
    deriv(ss, Ds);
    deriv(tt, Dt);
    for (int i = 0; i < 9; ++i) {
      J[i] = 0.0;
    }
  }

  xyz[0] = xyz[1] = xyz[2] = 0.0;
  const dlong off = e * Np_;
  for (int c = 0; c < Nq_; ++c) {
    for (int b = 0; b < Nq_; ++b) {
      const double wbc = Ls[b] * Lt[c];
      const dlong row = off + (c * Nq_ + b) * Nq_;
      for (int a = 0; a < Nq_; ++a) {
        const double w = Lr[a] * wbc;
        xyz[0] += w * gx_[row + a];
        xyz[1] += w * gy_[row + a];
        xyz[2] += w * gz_[row + a];
        if (J) {
          const double wr = Dr[a] * Ls[b] * Lt[c];
          const double ws = Lr[a] * Ds[b] * Lt[c];
          const double wt = Lr[a] * Ls[b] * Dt[c];
          J[0] += wr * gx_[row + a];
          J[1] += ws * gx_[row + a];
          J[2] += wt * gx_[row + a];
          J[3] += wr * gy_[row + a];
          J[4] += ws * gy_[row + a];
          J[5] += wt * gy_[row + a];
          J[6] += wr * gz_[row + a];
          J[7] += ws * gz_[row + a];
          J[8] += wt * gz_[row + a];
        }
      }
    }
  }
}

// =========================================================================== //
//  gradient at a root, in the REFERENCE frame.
//
//  Along the fibre the derivative is p'(root) -- ALREADY computed by the Aberth
//  iteration, which forms p/p' at every step; it is passed in as `slope`.
//  The two transverse derivatives come from one Lagrange contraction of nodal
//  values already held, followed by a Legendre derivative: a weighted sum of
//  data in hand, NOT a new field evaluation (report section 6.6).
// =========================================================================== //
inline void
isoSurface::gradAtRoot(dlong e, int dir, int j, int k, double rr, double slope, double *gref) const
{
  const dlong off = e * Np_;
  double L[isoDetail::maxNq];
  isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, rr, L);

  double line[isoDetail::maxNq], a[isoDetail::maxNq], v, d;

  // index of node (ia,ib,ic) inside element e:  off + (ic*Nq + ib)*Nq + ia
  auto at = [&](int ia, int ib, int ic) { return fld_[off + (ic * Nq_ + ib) * Nq_ + ia]; };

  if (dir == 0) { // fibre along r at (s_j, t_k)
    gref[0] = slope;
    for (int b = 0; b < Nq_; ++b) { // contract r at the root -> a line in s
      double sum = 0.0;
      for (int ia = 0; ia < Nq_; ++ia) {
        sum += L[ia] * at(ia, b, k);
      }
      line[b] = sum;
    }
    isoDetail::nodalToModal(M_.data(), line, Nq_, a);
    isoDetail::legendreValDer(a, N_, xn_[j], v, d);
    gref[1] = d;
    for (int c = 0; c < Nq_; ++c) { // ... and a line in t
      double sum = 0.0;
      for (int ia = 0; ia < Nq_; ++ia) {
        sum += L[ia] * at(ia, j, c);
      }
      line[c] = sum;
    }
    isoDetail::nodalToModal(M_.data(), line, Nq_, a);
    isoDetail::legendreValDer(a, N_, xn_[k], v, d);
    gref[2] = d;
  } else if (dir == 1) { // fibre along s at (r_j, t_k)
    gref[1] = slope;
    for (int ia = 0; ia < Nq_; ++ia) {
      double sum = 0.0;
      for (int b = 0; b < Nq_; ++b) {
        sum += L[b] * at(ia, b, k);
      }
      line[ia] = sum;
    }
    isoDetail::nodalToModal(M_.data(), line, Nq_, a);
    isoDetail::legendreValDer(a, N_, xn_[j], v, d);
    gref[0] = d;
    for (int c = 0; c < Nq_; ++c) {
      double sum = 0.0;
      for (int b = 0; b < Nq_; ++b) {
        sum += L[b] * at(j, b, c);
      }
      line[c] = sum;
    }
    isoDetail::nodalToModal(M_.data(), line, Nq_, a);
    isoDetail::legendreValDer(a, N_, xn_[k], v, d);
    gref[2] = d;
  } else { // fibre along t at (r_j, s_k)
    gref[2] = slope;
    for (int ia = 0; ia < Nq_; ++ia) {
      double sum = 0.0;
      for (int c = 0; c < Nq_; ++c) {
        sum += L[c] * at(ia, k, c);
      }
      line[ia] = sum;
    }
    isoDetail::nodalToModal(M_.data(), line, Nq_, a);
    isoDetail::legendreValDer(a, N_, xn_[j], v, d);
    gref[0] = d;
    for (int b = 0; b < Nq_; ++b) {
      double sum = 0.0;
      for (int c = 0; c < Nq_; ++c) {
        sum += L[c] * at(j, b, c);
      }
      line[b] = sum;
    }
    isoDetail::nodalToModal(M_.data(), line, Nq_, a);
    isoDetail::legendreValDer(a, N_, xn_[k], v, d);
    gref[1] = d;
  }
}

// =========================================================================== //
//  allocDevice -- one allocation pass, reused for every dump.
//
//  Sizes come from IsoCaps, so they depend only on N and prm.capScale, not on
//  the field.  The basis tables are uploaded here too: they are constant for
//  the run, and re-uploading them per dump would be pure waste.
// =========================================================================== //
inline void isoSurface::allocDevice()
{
  const dlong Nel = mesh_->Nelements;
  const int nLev = static_cast<int>(prm_.levels.size());

  // basis tables.  dfloat may be float in this build, so narrow explicitly
  // rather than relying on an implicit conversion inside copyFrom.
  auto upload = [](deviceMemory<dfloat> &o, const std::vector<double> &h) {
    std::vector<dfloat> tmp(h.size());
    for (size_t i = 0; i < h.size(); ++i) {
      tmp[i] = static_cast<dfloat>(h[i]);
    }
    o.resize(tmp.size());
    o.copyFrom(tmp.data(), tmp.size());
  };
  upload(o_xn_, xn_);
  upload(o_bw_, bw_);
  upload(o_M_, M_);

  std::vector<dfloat> lv(prm_.levels.begin(), prm_.levels.end());
  o_levels_.resize(std::max<size_t>(lv.size(), 1));
  if (!lv.empty()) {
    o_levels_.copyFrom(lv.data(), lv.size());
  }

  o_fld_.resize(mesh_->Nlocal);
  o_elemCross_.resize(static_cast<size_t>(Nel) * std::max(nLev, 1));
  o_overflow_.resize(1);

  const size_t nR = static_cast<size_t>(Nel) * caps_.root;
  o_rootCount_.resize(Nel);
  o_rootDir_.resize(nR);
  o_rootJK_.resize(nR);
  o_rootR_.resize(nR);
  o_rootS_.resize(nR);
  o_rootT_.resize(nR);
  o_rootSlope_.resize(nR);
  o_rootFlag_.resize(nR);

  const size_t nP = static_cast<size_t>(Nel) * caps_.point;
  const size_t nC = static_cast<size_t>(Nel) * caps_.cell;
  o_cellCount_.resize(Nel);
  o_ptCount_.resize(Nel);
  o_conn_.resize(nC * 6); // room for P2; P1 uses the first 3 of each
  o_ptRoot_.resize(nP);
  o_ptR_.resize(nP);
  o_ptS_.resize(nP);
  o_ptT_.resize(nP);
  o_ptFlag_.resize(nP);
  o_cellFlag_.resize(nC);
  o_midCount_.resize(Nel);
  o_midBase_.resize(Nel);
  o_midLock_.resize(nP);
  o_ptTotal_.resize(Nel);
  o_nUnsep_.resize(1);
  o_nDegen_.resize(1);
  o_projResidQ_.resize(7); // one counter per decade of residual
  o_nProjShort_.resize(1);

  o_ptOffset_.resize(Nel);
  o_outElem_.resize(nR);
  o_outRoot_.resize(nR);
  o_outR_.resize(nR);
  o_outS_.resize(nR);
  o_outT_.resize(nR);
  o_outFlag_.resize(nR);
}

//  checkOverflow -- a capacity overflow is a WRONG ANSWER, not a warning.
//
//  The kernels increment o_overflow_ rather than setting a bit (CUDA can only
//  transform @atomic on a compound assignment -- notes/22), so the value is the
//  number of slots lost.  Aborting names the stage and what to do, because a
//  truncated surface looks plausible in ParaView and would be believed.
inline void isoSurface::checkOverflow(const char *stage)
{
  std::vector<int> ovf(1, 0);
  o_overflow_.copyTo(ovf, 1);
  long long g = ovf[0];
  MPI_Allreduce(MPI_IN_PLACE, &g, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm());
  nekrsCheck(g != 0,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "isoSurface: %s overflowed the per-element capacity, dropping %lld\n"
             "  root(s)/triangle(s).  The surface has holes.\n"
             "\n"
             "  Both capacities are BOUNDS derived from the algorithm, so neither\n"
             "  should ever overflow at capScale = 1 (currently %.2f):\n"
             "    capRoot = capScale * 3N(N+1)^2  -- %d fibres x at most N roots\n"
             "    capCell = p_isoMaxTri * N^3     -- a fan of <= 12 cuts gives <= 10 tris\n"
             "  Reaching this means a bound is wrong, which is a BUG -- please report\n"
             "  it with the order (N = %d) and the field.  Raising capScale will get\n"
             "  you moving again, and it is a runtime argument so no rebuild is\n"
             "  needed, but it is a workaround and the holes may recur.\n",
             stage,
             g,
             static_cast<double>(prm_.capScale),
             3 * (N_ + 1) * (N_ + 1),
             N_);
}

// =========================================================================== //
//  solveRoots -- STEP 1 of plan F: the fibre solve runs on the DEVICE.
//
//      isoCullOKL          the theta*K_N band, one pass over the nodes
//      isoSolveFibresOKL   one fibre per thread, fixed-iteration Aberth
//      isoCompactPointsOKL gather the sparse per-element slabs into dense
//
//  Then the host does what no kernel produces: the physical coordinates, both
//  gradients, and the MultiEdge flag.  That round trip is deliberate -- see the
//  file header; the kernels do the O(fibres) work and the host does O(roots).
// =========================================================================== //
inline void isoSurface::solveRoots(bool force)
{
  if (rootsValid_ && !force) {
    if (platform->comm.mpiRank() == 0 && platform->verbose()) {
      printf("isoSurface(%s): roots still valid, skipping the solve\n", fieldName_.c_str());
    }
    return;
  }

  platform->timer.tic("isoSurface::solve");

  // The caller's device array may have been overwritten in place since the last
  // solve, so re-read it.  Everything below -- the host geometry pass and the
  // o_fld_ upload -- reads fld_.
  refreshFieldHost();

  const dlong solveCalls = stats_.nSolveCalls;
  stats_ = IsoStats();
  stats_.nSolveCalls = solveCalls + 1;
  stats_.nElements = mesh_->Nelements;

  const dlong Nel = mesh_->Nelements;
  const int nLev = static_cast<int>(prm_.levels.size());

  // the scalar, host -> device (fld_ is kept on the host for the geometry pass)
  {
    std::vector<dfloat> tmp(mesh_->Nlocal);
    for (dlong i = 0; i < mesh_->Nlocal; ++i) {
      tmp[i] = static_cast<dfloat>(fld_[i]);
    }
    o_fld_.copyFrom(tmp.data(), tmp.size());
  }
  {
    std::vector<int> zero(1, 0);
    o_overflow_.copyFrom(zero.data(), 1);
  }

  // ---- cull, all levels at once -------------------------------------------
  isoCullOKL(Nel,
             static_cast<int>(nLev),
             static_cast<int>(N_),
             static_cast<dfloat>(prm_.cullElements ? prm_.cullTheta : -1.0),
             o_levels_,
             o_fld_,
             o_elemCross_);

  std::vector<int> cross(static_cast<size_t>(Nel) * std::max(nLev, 1));
  o_elemCross_.copyTo(cross, cross.size());

  std::vector<int> rc(Nel), rdir, rjk;
  std::vector<dfloat> rr, rs, rt, rsl;
  std::vector<unsigned char> rfl;

  for (int l = 0; l < nLev; ++l) {
    auto &rset = rootSets_[l];
    rset.clear();
    rset.level = prm_.levels[l];
    rset.eCount.assign(Nel, 0);
    rset.eStart.assign(Nel + 1, 0);

    // ---- the fibre solve ---------------------------------------------------
    isoSolveFibresOKL(Nel,
                      static_cast<int>(l),
                      static_cast<int>(nLev),
                      static_cast<int>(N_),
                      static_cast<dfloat>(prm_.levels[l]),
                      caps_.root,
                      static_cast<int>(prm_.niter),
                      static_cast<dfloat>(prm_.delta),
                      static_cast<dfloat>(prm_.aFloor),
                      static_cast<dfloat>(prm_.imagTol),
                      static_cast<dfloat>(prm_.margin),
                      static_cast<dfloat>(prm_.gradRelTol),
                      o_elemCross_,
                      o_xn_,
                      o_bw_,
                      o_M_,
                      o_fld_,
                      o_rootCount_,
                      o_rootR_,
                      o_rootS_,
                      o_rootT_,
                      o_rootSlope_,
                      o_rootDir_,
                      o_rootJK_,
                      o_rootFlag_,
                      o_overflow_);
    checkOverflow("the fibre solve");

    // ---- counts back, host prefix sum, then compact ------------------------
    o_rootCount_.copyTo(rc, Nel);
    std::vector<dlong> off(Nel, 0); // dlong: matches the kernel's ptOffset type
    int total = 0;
    for (dlong e = 0; e < Nel; ++e) {
      off[e] = total;
      total += rc[e];
      rset.eCount[e] = rc[e];
      if (cross[e * nLev + l]) {
        ++stats_.nElementsSolved;
      }
    }
    o_ptOffset_.copyFrom(off.data(), off.size());

    isoCompactPointsOKL(Nel,
                        caps_.root,
                        o_rootCount_,
                        o_ptOffset_,
                        o_rootR_,
                        o_rootS_,
                        o_rootT_,
                        o_rootDir_,
                        o_rootFlag_,
                        o_outR_,
                        o_outS_,
                        o_outT_,
                        o_outElem_,
                        o_outRoot_,
                        o_outFlag_);

    if (total == 0) {
      continue;
    }

    rr.resize(total);
    rs.resize(total);
    rt.resize(total);
    rsl.resize(total);
    rdir.assign(total, 0);
    rjk.assign(total, 0);
    rfl.resize(total);
    o_outR_.copyTo(rr, total);
    o_outS_.copyTo(rs, total);
    o_outT_.copyTo(rt, total);
    // outElem/outRoot are the ELEMENT ID and the slot within it -- not dir/jk.
    // dir and (j,k) stay in the sparse slabs and are gathered below.
    std::vector<int> oelem(total), oroot(total);
    o_outElem_.copyTo(oelem, total);
    o_outRoot_.copyTo(oroot, total);
    o_outFlag_.copyTo(rfl, total);
    // the slope and (dir, j, k) still live in the sparse slabs; pull the ones
    // we need through the same offsets
    {
      std::vector<dfloat> sslope(static_cast<size_t>(Nel) * caps_.root);
      std::vector<int> sdir(sslope.size()), sjk(sslope.size());
      o_rootSlope_.copyTo(sslope, sslope.size());
      o_rootDir_.copyTo(sdir, sdir.size());
      o_rootJK_.copyTo(sjk, sjk.size());
      for (dlong e = 0; e < Nel; ++e) {
        for (int q = 0; q < rc[e]; ++q) {
          const size_t src = static_cast<size_t>(e) * caps_.root + q;
          const int dst = off[e] + q;
          rsl[dst] = sslope[src];
          rdir[dst] = sdir[src];
          rjk[dst] = sjk[src];
        }
      }
    }

    // ---- THE HOST GAP: physical coordinates, gradients, MultiEdge ----------
    rset.reserve(total);
    for (dlong e = 0; e < Nel; ++e) {
      const int n = rc[e];
      for (int q = 0; q < n; ++q) {
        const int i = off[e] + q;
        const double rc_ = rr[i], sc_ = rs[i], tc_ = rt[i];
        const int dir = rdir[i];
        const int jj = rjk[i] / Nq_, kk = rjk[i] % Nq_;
        const double along = (dir == 0) ? rc_ : ((dir == 1) ? sc_ : tc_);

        double gref[3] = {0, 0, 0};
        gradAtRoot(e, dir, jj, kk, along, static_cast<double>(rsl[i]), gref);
        double xyz[3], Jm[9];
        mapToPhysical(e, rc_, sc_, tc_, xyz, Jm);

        // grad_x u = J^-T grad_r u.  Only writeRoots() consumes this; the
        // SURFACE normal is computed independently in pointNormal().
        double gphys[3] = {0, 0, 0};
        {
          const double det = Jm[0] * (Jm[4] * Jm[8] - Jm[5] * Jm[7]) -
                             Jm[1] * (Jm[3] * Jm[8] - Jm[5] * Jm[6]) +
                             Jm[2] * (Jm[3] * Jm[7] - Jm[4] * Jm[6]);
          if (std::fabs(det) > 1e-300) {
            const double iJ[9] = {(Jm[4] * Jm[8] - Jm[5] * Jm[7]) / det,
                                  (Jm[2] * Jm[7] - Jm[1] * Jm[8]) / det,
                                  (Jm[1] * Jm[5] - Jm[2] * Jm[4]) / det,
                                  (Jm[5] * Jm[6] - Jm[3] * Jm[8]) / det,
                                  (Jm[0] * Jm[8] - Jm[2] * Jm[6]) / det,
                                  (Jm[2] * Jm[3] - Jm[0] * Jm[5]) / det,
                                  (Jm[3] * Jm[7] - Jm[4] * Jm[6]) / det,
                                  (Jm[1] * Jm[6] - Jm[0] * Jm[7]) / det,
                                  (Jm[0] * Jm[4] - Jm[1] * Jm[3]) / det};
            gphys[0] = iJ[0] * gref[0] + iJ[3] * gref[1] + iJ[6] * gref[2];
            gphys[1] = iJ[1] * gref[0] + iJ[4] * gref[1] + iJ[7] * gref[2];
            gphys[2] = iJ[2] * gref[0] + iJ[5] * gref[1] + iJ[8] * gref[2];
          }
        }

        uint8_t fl = rfl[i]; // Tangency | OnFace, set by the kernel
        // MultiEdge: the kernel does not set it.  The roots of one fibre are
        // adjacent in the dense list and sorted, so a neighbour in the same
        // sub-interval is the previous entry with the same (dir, j, k).
        if (q > 0 && rdir[i - 1] == dir && rjk[i - 1] == rjk[i]) {
          const double prev = (dir == 0) ? rr[i - 1] : ((dir == 1) ? rs[i - 1] : rt[i - 1]);
          int ea = 0, eb = 0;
          while (ea < N_ && prev > xn_[ea + 1]) {
            ++ea;
          }
          while (eb < N_ && along > xn_[eb + 1]) {
            ++eb;
          }
          if (ea == eb) {
            fl |= RootSet::MultiEdge;
            ++stats_.nMultiEdge;
          }
        }
        if (fl & RootSet::Tangency) {
          ++stats_.nTangency;
        }
        if (fl & RootSet::OnFace) {
          ++stats_.nOnFace;
        }

        rset.elem.push_back(e);
        rset.eg.push_back(egOffset_ + static_cast<hlong>(e));
        rset.dir.push_back(static_cast<uint8_t>(dir));
        rset.j.push_back(static_cast<int16_t>(jj));
        rset.k.push_back(static_cast<int16_t>(kk));
        rset.r.push_back(static_cast<dfloat>(rc_));
        rset.s.push_back(static_cast<dfloat>(sc_));
        rset.t.push_back(static_cast<dfloat>(tc_));
        rset.x.push_back(static_cast<dfloat>(xyz[0]));
        rset.y.push_back(static_cast<dfloat>(xyz[1]));
        rset.z.push_back(static_cast<dfloat>(xyz[2]));
        rset.dfdr.push_back(static_cast<dfloat>(gref[0]));
        rset.dfds.push_back(static_cast<dfloat>(gref[1]));
        rset.dfdt.push_back(static_cast<dfloat>(gref[2]));
        rset.dfdx.push_back(static_cast<dfloat>(gphys[0]));
        rset.dfdy.push_back(static_cast<dfloat>(gphys[1]));
        rset.dfdz.push_back(static_cast<dfloat>(gphys[2]));
        rset.flag.push_back(fl);
        ++rset.nRoots;
      }
    }
    for (dlong e = 0; e < Nel; ++e) {
      rset.eStart[e + 1] = rset.eStart[e] + rset.eCount[e];
    }
    stats_.nRoots += rset.nRoots;

    // per-element occupancy, over the elements that actually cross
    for (dlong e = 0; e < Nel; ++e) {
      const dlong n = rset.eCount[e];
      if (n <= 0) {
        continue;
      }
      ++stats_.nElemCrossed;
      stats_.maxRootPerElem = std::max(stats_.maxRootPerElem, n);
      stats_.minRootPerElem = (stats_.minRootPerElem == 0) ? n : std::min(stats_.minRootPerElem, n);
    }
  }

  rootsValid_ = true;
  platform->timer.toc("isoSurface::solve");

  hlong gRoots = stats_.nRoots, gSolved = stats_.nElementsSolved, gEl = stats_.nElements;
  hlong gTan = stats_.nTangency, gMul = stats_.nMultiEdge;
  MPI_Allreduce(MPI_IN_PLACE, &gRoots, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gSolved, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gEl, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gTan, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gMul, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());

  if (platform->comm.mpiRank() == 0) {
    const double cull = (gEl * nLev > 0) ? 100.0 * (1.0 - double(gSolved) / (gEl * nLev)) : 0.0;
    printf("isoSurface(%s): [GPU] %lld roots, %lld/%lld element-levels solved (%.1f%% culled)",
           fieldName_.c_str(),
           (long long)gRoots,
           (long long)gSolved,
           (long long)(gEl * nLev),
           cull);
    if (gTan || gMul) {
      printf(", %lld tangency, %lld multi-edge", (long long)gTan, (long long)gMul);
    }
    printf("\n");
  }
  reportOccupancy("roots",
                  stats_.nElemCrossed,
                  stats_.minRootPerElem,
                  stats_.maxRootPerElem,
                  stats_.nRoots,
                  caps_.root);
}

// --------------------------------------------------------------------------- //
//  reportOccupancy -- min/max/avg per element against the capacity.
//
//  A total says nothing about the element that overflows, and two capacity
//  aborts reached a user before this existed (notes/24, notes/25).  What matters
//  is the WORST element as a fraction of the cap: at 30% a fold in the field is
//  harmless, at 90% the next timestep may abort.
//
//  The average is over elements that actually CROSS the level.  Including the
//  empty ones makes it a statement about how much of the mesh the surface
//  touches, which is a different question and hides the crowding.
//
//  Collective: every rank must call it.
// --------------------------------------------------------------------------- //
inline void isoSurface::reportOccupancy(const char *what,
                                        dlong nCrossed,
                                        dlong locMin,
                                        dlong locMax,
                                        dlong locSum,
                                        dlong cap) const
{
  MPI_Comm comm = platform->comm.mpiComm();
  hlong gCrossed = nCrossed, gSum = locSum;
  hlong gMax = locMax;
  // a rank with no crossed element must not win the min
  hlong gMin = (nCrossed > 0) ? locMin : std::numeric_limits<hlong>::max();
  MPI_Allreduce(MPI_IN_PLACE, &gCrossed, 1, MPI_HLONG, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &gSum, 1, MPI_HLONG, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &gMax, 1, MPI_HLONG, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &gMin, 1, MPI_HLONG, MPI_MIN, comm);

  if (platform->comm.mpiRank() != 0 || gCrossed == 0) {
    return;
  }
  const double avg = double(gSum) / double(gCrossed);
  const double pk = (cap > 0) ? 100.0 * double(gMax) / double(cap) : 0.0;
  printf("isoSurface(%s): %s/element over %lld crossed: min %lld, avg %.1f, max %lld"
         "  (cap %lld, peak %.0f%%)%s\n",
         fieldName_.c_str(),
         what,
         (long long)gCrossed,
         (long long)gMin,
         avg,
         (long long)gMax,
         (long long)cap,
         pk,
         (pk >= 80.0) ? "  [!] near capacity" : "");
}

//  uploadRoots -- put ONE level's roots back into the per-element slabs.
//
//  solveRoots leaves only the LAST level's roots on the device, so a multi-level
//  run would otherwise march the wrong crossings.  The RootSet is dense and
//  element-major; the kernels want the sparse e*capRoot layout.
inline void isoSurface::uploadRoots(const RootSet &rs)
{
  const dlong Nel = mesh_->Nelements;
  const size_t slab = static_cast<size_t>(Nel) * caps_.root;
  std::vector<dfloat> rr(slab, 0), rsv(slab, 0), rt(slab, 0);
  std::vector<int> rdir(slab, 0), rjk(slab, 0), rcount(Nel, 0);
  std::vector<unsigned char> rfl(slab, 0);

  for (dlong e = 0; e < Nel; ++e) {
    const dlong a = rs.eStart[e], b = rs.eStart[e + 1];
    rcount[e] = static_cast<int>(b - a);
    for (dlong q = a; q < b; ++q) {
      const size_t o = static_cast<size_t>(e) * caps_.root + (q - a);
      rr[o] = rs.r[q];
      rsv[o] = rs.s[q];
      rt[o] = rs.t[q];
      rdir[o] = rs.dir[q];
      rjk[o] = rs.j[q] * Nq_ + rs.k[q];
      rfl[o] = rs.flag[q];
    }
  }
  o_rootCount_.copyFrom(rcount.data(), rcount.size());
  o_rootR_.copyFrom(rr.data(), rr.size());
  o_rootS_.copyFrom(rsv.data(), rsv.size());
  o_rootT_.copyFrom(rt.data(), rt.size());
  o_rootDir_.copyFrom(rdir.data(), rdir.size());
  o_rootJK_.copyFrom(rjk.data(), rjk.size());
  o_rootFlag_.copyFrom(rfl.data(), rfl.size());
}

//  pointNormal -- the UNIT physical normal at a surface point.
//
//  grad_x f = J^{-T} grad_r f, with grad_r from one tensor contraction of the
//  element's nodal values.  The kernels do not produce this: they emit reference
//  coordinates only, and every writer reads ss.nx/ny/nz.
inline void isoSurface::pointNormal(dlong e, double r, double s, double t, const double *J, double *nrm) const
{
  double Lr[isoDetail::maxNq], Ls[isoDetail::maxNq], Lt[isoDetail::maxNq];
  double Dr[isoDetail::maxNq], Ds[isoDetail::maxNq], Dt[isoDetail::maxNq], g[3];
  isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, r, Lr);
  isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, s, Ls);
  isoDetail::baryWeightsAt(xn_.data(), bw_.data(), Nq_, t, Lt);
  isoDetail::baryDerivAt(xn_.data(), bw_.data(), Nq_, r, Dr);
  isoDetail::baryDerivAt(xn_.data(), bw_.data(), Nq_, s, Ds);
  isoDetail::baryDerivAt(xn_.data(), bw_.data(), Nq_, t, Dt);
  isoDetail::tensorEvalGrad3(&fld_[e * Np_], Nq_, Lr, Ls, Lt, Dr, Ds, Dt, g);

  double gp[3] = {g[0], g[1], g[2]};
  const double det = J[0] * (J[4] * J[8] - J[5] * J[7]) - J[1] * (J[3] * J[8] - J[5] * J[6]) +
                     J[2] * (J[3] * J[7] - J[4] * J[6]);
  if (std::fabs(det) > 1e-300) {
    const double iJ[9] = {(J[4] * J[8] - J[5] * J[7]) / det,
                          (J[2] * J[7] - J[1] * J[8]) / det,
                          (J[1] * J[5] - J[2] * J[4]) / det,
                          (J[5] * J[6] - J[3] * J[8]) / det,
                          (J[0] * J[8] - J[2] * J[6]) / det,
                          (J[2] * J[3] - J[0] * J[5]) / det,
                          (J[3] * J[7] - J[4] * J[6]) / det,
                          (J[1] * J[6] - J[0] * J[7]) / det,
                          (J[0] * J[4] - J[1] * J[3]) / det};
    gp[0] = iJ[0] * g[0] + iJ[3] * g[1] + iJ[6] * g[2];
    gp[1] = iJ[1] * g[0] + iJ[4] * g[1] + iJ[7] * g[2];
    gp[2] = iJ[2] * g[0] + iJ[5] * g[1] + iJ[8] * g[2];
  }
  const double gn = std::sqrt(gp[0] * gp[0] + gp[1] * gp[1] + gp[2] * gp[2]);
  const double inv = (gn > 1e-300) ? 1.0 / gn : 0.0;
  nrm[0] = gp[0] * inv;
  nrm[1] = gp[1] * inv;
  nrm[2] = gp[2] * inv;
}

// =========================================================================== //
//  buildTopology -- STEP 2 of plan F: marching cubes runs on the DEVICE.
//
//      isoTopologyOKL      the N^3 sub-cell march, corners are the roots
//      (host)              the P2 midside nodes: the unique-edge dedup, the
//                          conn 3->6 growth and the face-lock mask.  No kernel
//                          does this -- isoProjectOKL takes midCount/midLock as
//                          GIVEN inputs.
//      isoProjectOKL       the fixed-sweep Newton projection, face-locked
//      (host)              mapToPhysical for every point, and the normals.
//                          The kernels emit only reference (r,s,t).
//
//  The device keeps a FIXED lattice where v1 inserts a plane for a fibre with
//  two roots in one interval.  Both crossings are still kept: the ring walk
//  connects them (DEVELOPER.md).  nUnseparated counts such edges -- a measure of
//  how folded the field is, not of lost data.
// =========================================================================== //
inline void isoSurface::buildTopology(bool force)
{
  if (!rootsValid_) {
    solveRoots();
  }
  if (surfValid_ && !force) {
    if (platform->comm.mpiRank() == 0 && platform->verbose()) {
      printf("isoSurface(%s): surface still valid, skipping the rebuild\n", fieldName_.c_str());
    }
    return;
  }

  platform->timer.tic("isoSurface::topology");

  nekrsCheck(prm_.cellOrder < 1 || prm_.cellOrder > 2,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "isoSurface: cellOrder must be 1 or 2, got %d\n",
             prm_.cellOrder);

  stats_.nCell = 0;
  stats_.nCellAmbiguous = 0;
  stats_.nMidNode = 0;
  stats_.nProjDegenerate = 0;
  stats_.maxProjResid = 0.0;

  const dlong Nel = mesh_->Nelements;
  const int nLev = static_cast<int>(prm_.levels.size());
  const int npc = (prm_.cellOrder == 2) ? 6 : 3;

  {
    std::vector<int> zero(1, 0);
    o_overflow_.copyFrom(zero.data(), 1);
    o_nUnsep_.copyFrom(zero.data(), 1);
    o_nDegen_.copyFrom(zero.data(), 1);
    o_nProjShort_.copyFrom(zero.data(), 1);
    std::vector<int> zero7(7, 0);
    o_projResidQ_.copyFrom(zero7.data(), 7);
  }

  for (int l = 0; l < nLev; ++l) {
    SurfaceSet &ss = surfSets_[l];
    const RootSet &rs = rootSets_[l];
    ss.clear();
    ss.level = prm_.levels[l];
    ss.cellOrder = prm_.cellOrder;
    ss.cellType = (prm_.cellOrder == 2) ? 22 : 5;
    ss.nodesPerCell = npc;

    // the roots for THIS level have to be back on the device: solveRoots left
    // the last level's there, and a multi-level run would otherwise march the
    // wrong crossings.
    uploadRoots(rs);

    isoTopologyOKL(Nel,
                   static_cast<int>(N_),
                   static_cast<dfloat>(prm_.levels[l]),
                   caps_.root,
                   caps_.cell,
                   caps_.point,
                   static_cast<int>(prm_.cellOrder),
                   o_xn_,
                   o_bw_,
                   o_fld_,
                   o_rootCount_,
                   o_rootR_,
                   o_rootS_,
                   o_rootT_,
                   o_rootDir_,
                   o_rootJK_,
                   o_rootFlag_,
                   o_cellCount_,
                   o_ptCount_,
                   o_conn_,
                   o_ptR_,
                   o_ptS_,
                   o_ptT_,
                   o_ptRoot_,
                   o_ptFlag_,
                   o_cellFlag_,
                   o_nUnsep_,
                   o_overflow_);
    checkOverflow("marching cubes");

    // ---- pull the per-element counts and the corner points -----------------
    std::vector<int> pc(Nel), cc(Nel);
    o_ptCount_.copyTo(pc, Nel);
    o_cellCount_.copyTo(cc, Nel);

    const size_t slabP = static_cast<size_t>(Nel) * caps_.point;
    const size_t slabC = static_cast<size_t>(Nel) * caps_.cell;
    std::vector<dfloat> pr(slabP), ps(slabP), pt(slabP);
    std::vector<int> proot(slabP), conn(slabC * 6);
    std::vector<unsigned char> pflag(slabP), cflag(slabC);
    o_ptR_.copyTo(pr, slabP);
    o_ptS_.copyTo(ps, slabP);
    o_ptT_.copyTo(pt, slabP);
    o_ptRoot_.copyTo(proot, slabP);
    o_ptFlag_.copyTo(pflag, slabP);
    o_conn_.copyTo(conn, slabC * 6);
    o_cellFlag_.copyTo(cflag, slabC);

    // ---- HOST: the P2 midside nodes ----------------------------------------
    // Hash on the ordered corner pair so an edge shared by two triangles gets
    // ONE node: that sharing is what makes the cost 2.76 new points per root
    // rather than 3 per triangle.  No kernel does this.
    std::vector<int> midCount(Nel, 0), midBase(Nel, 0);
    std::vector<dfloat> midR(slabP, 0), midS(slabP, 0), midT(slabP, 0);
    std::vector<int> midLock(slabP, 0);
    std::vector<int> conn6(slabC * 6, -1);

    if (prm_.cellOrder == 2) {
      for (dlong e = 0; e < Nel; ++e) {
        if (cc[e] == 0) {
          continue;
        }
        midBase[e] = pc[e]; // midside nodes start after this element's corners
        std::unordered_map<unsigned long long, int> edgeNode;
        edgeNode.reserve(static_cast<size_t>(3 * cc[e]));
        int nMid = 0;
        for (int c = 0; c < cc[e]; ++c) {
          const size_t src = (static_cast<size_t>(e) * caps_.cell + c) * 6;
          const int v[3] = {conn[src + 0], conn[src + 1], conn[src + 2]};
          conn6[src + 0] = v[0];
          conn6[src + 1] = v[1];
          conn6[src + 2] = v[2];
          // VTK_QUADRATIC_TRIANGLE: node 3 = mid(0,1), 4 = mid(1,2), 5 = mid(2,0)
          for (int en = 0; en < 3; ++en) {
            const int a = v[en], b = v[(en + 1) % 3];
            const int lo = std::min(a, b), hi = std::max(a, b);
            const unsigned long long key =
                (static_cast<unsigned long long>(lo) << 32) | static_cast<unsigned long long>(hi);
            auto it = edgeNode.find(key);
            int node;
            if (it == edgeNode.end()) {
              node = midBase[e] + nMid;
              edgeNode.emplace(key, node);
              const size_t oa = static_cast<size_t>(e) * caps_.point + a;
              const size_t ob = static_cast<size_t>(e) * caps_.point + b;
              const size_t om = static_cast<size_t>(e) * caps_.point + nMid;
              midR[om] = 0.5 * (pr[oa] + pr[ob]);
              midS[om] = 0.5 * (ps[oa] + ps[ob]);
              midT[om] = 0.5 * (pt[oa] + pt[ob]);
              // WATERTIGHTNESS: if both parents sit on the same reference face,
              // the node must stay ON it, or the two elements sharing the face
              // walk their copies apart -- a visible crack (notes/17).
              int lk = 0;
              const dfloat ca[3] = {pr[oa], ps[oa], pt[oa]};
              const dfloat cb[3] = {pr[ob], ps[ob], pt[ob]};
              for (int d = 0; d < 3; ++d) {
                if (std::fabs(std::fabs(ca[d]) - 1.0) <= static_cast<double>(prm_.margin) &&
                    std::fabs(ca[d] - cb[d]) <= static_cast<double>(prm_.margin)) {
                  lk |= (1 << d);
                }
              }
              midLock[om] = lk;
              ++nMid;
            } else {
              node = it->second;
            }
            conn6[src + 3 + en] = node;
          }
        }
        midCount[e] = nMid;
      }

      // the midside starts go to the device, then the projection
      o_midCount_.copyFrom(midCount.data(), midCount.size());
      o_midBase_.copyFrom(midBase.data(), midBase.size());
      o_midLock_.copyFrom(midLock.data(), midLock.size());
      // seed the slab slots the kernel will refine
      for (dlong e = 0; e < Nel; ++e) {
        for (int m = 0; m < midCount[e]; ++m) {
          const size_t src = static_cast<size_t>(e) * caps_.point + m;
          const size_t dst = static_cast<size_t>(e) * caps_.point + midBase[e] + m;
          pr[dst] = midR[src];
          ps[dst] = midS[src];
          pt[dst] = midT[src];
          pflag[dst] = 0;
        }
      }
      o_ptR_.copyFrom(pr.data(), pr.size());
      o_ptS_.copyFrom(ps.data(), ps.size());
      o_ptT_.copyFrom(pt.data(), pt.size());
      o_ptFlag_.copyFrom(pflag.data(), pflag.size());

      isoProjectOKL(Nel,
                    caps_.point,
                    static_cast<int>(prm_.niterProj),
                    o_midBase_,
                    static_cast<dfloat>(prm_.gradRelTol),
                    static_cast<dfloat>(prm_.levels[l]),
                    o_xn_,
                    o_bw_,
                    o_fld_,
                    o_midCount_,
                    o_midLock_,
                    o_ptR_,
                    o_ptS_,
                    o_ptT_,
                    o_ptFlag_,
                    o_nDegen_,
                    o_projResidQ_,
                    o_nProjShort_);

      o_ptR_.copyTo(pr, slabP);
      o_ptS_.copyTo(ps, slabP);
      o_ptT_.copyTo(pt, slabP);
      o_ptFlag_.copyTo(pflag, slabP);
    }

    // ---- HOST: physical coordinates and normals, then assemble the SoA -----
    for (dlong e = 0; e < Nel; ++e) {
      const int nCorner = pc[e];
      const int nTotal = nCorner + (prm_.cellOrder == 2 ? midCount[e] : 0);
      for (int m = 0; m < nTotal; ++m) {
        const size_t o = static_cast<size_t>(e) * caps_.point + m;
        double xyz[3], Jm[9], g[3] = {0, 0, 0};
        mapToPhysical(e, pr[o], ps[o], pt[o], xyz, Jm);
        pointNormal(e, pr[o], ps[o], pt[o], Jm, g);
        ss.px.push_back(static_cast<dfloat>(xyz[0]));
        ss.py.push_back(static_cast<dfloat>(xyz[1]));
        ss.pz.push_back(static_cast<dfloat>(xyz[2]));
        ss.pr.push_back(pr[o]);
        ss.ps.push_back(ps[o]);
        ss.pt.push_back(pt[o]);
        ss.pelem.push_back(e);
        ss.nx.push_back(static_cast<dfloat>(g[0]));
        ss.ny.push_back(static_cast<dfloat>(g[1]));
        ss.nz.push_back(static_cast<dfloat>(g[2]));
        ss.cornerRoot.push_back(m < nCorner ? proot[o] : -1);
        ss.pflag.push_back(pflag[o]);
        ++ss.nPoint;
        if (m < nCorner) {
          ++ss.nCorner;
        }
      }
      stats_.nMidNode += (prm_.cellOrder == 2) ? midCount[e] : 0;
    }

    // the per-element TOTAL, for step 3's interpolation
    ptTotal_.assign(Nel, 0);
    for (dlong e = 0; e < Nel; ++e) {
      ptTotal_[e] = pc[e] + (prm_.cellOrder == 2 ? midCount[e] : 0);
    }
    o_ptTotal_.copyFrom(ptTotal_.data(), ptTotal_.size());

    // ---- assemble the cells, element-local ids -> dense -------------------
    std::vector<dlong> base(Nel, 0);
    {
      dlong acc = 0;
      for (dlong e = 0; e < Nel; ++e) {
        base[e] = acc;
        acc += pc[e] + (prm_.cellOrder == 2 ? midCount[e] : 0);
      }
    }
    for (dlong e = 0; e < Nel; ++e) {
      for (int c = 0; c < cc[e]; ++c) {
        const size_t src = (static_cast<size_t>(e) * caps_.cell + c) * 6;
        for (int m = 0; m < npc; ++m) {
          const int local = (prm_.cellOrder == 2) ? conn6[src + m] : conn[src + m];
          ss.conn.push_back(base[e] + local);
        }
        ss.cellEg.push_back(egOffset_ + static_cast<hlong>(e));
        unsigned char cf = cflag[static_cast<size_t>(e) * caps_.cell + c];
        if (cf & SurfaceSet::Ambiguous) {
          ++stats_.nCellAmbiguous;
        }
        // a midside node whose Newton step was zeroed marks its cell
        if (prm_.cellOrder == 2) {
          for (int en = 3; en < 6; ++en) {
            const int p = base[e] + conn6[src + en];
            if (p >= 0 && p < static_cast<int>(ss.pflag.size()) && (ss.pflag[p] & RootSet::Tangency)) {
              cf |= SurfaceSet::Tangency;
              break;
            }
          }
        }
        ss.cellFlag.push_back(cf);
        ++ss.nCell;
      }
    }
    stats_.nCell += ss.nCell;

    // per-element triangle occupancy, over the elements that carry a surface
    for (dlong e = 0; e < Nel; ++e) {
      const dlong n = cc[e];
      if (n <= 0) {
        continue;
      }
      ++stats_.nElemWithCell;
      stats_.maxCellPerElem = std::max(stats_.maxCellPerElem, n);
      stats_.minCellPerElem = (stats_.minCellPerElem == 0) ? n : std::min(stats_.minCellPerElem, n);
    }
  }

  surfValid_ = true;
  platform->timer.toc("isoSurface::topology");

  reportOccupancy("tris ",
                  stats_.nElemWithCell,
                  stats_.minCellPerElem,
                  stats_.maxCellPerElem,
                  stats_.nCell,
                  caps_.cell);

  std::vector<int> uns(1, 0), deg(1, 0);
  o_nUnsep_.copyTo(uns, 1);
  o_nDegen_.copyTo(deg, 1);

  hlong gCell = stats_.nCell, gAmb = stats_.nCellAmbiguous, gMid = stats_.nMidNode;
  hlong gUns = uns[0], gDeg = deg[0];
  MPI_Allreduce(MPI_IN_PLACE, &gCell, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gAmb, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gMid, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gUns, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  MPI_Allreduce(MPI_IN_PLACE, &gDeg, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());

  if (platform->comm.mpiRank() == 0) {
    printf("isoSurface(%s): [GPU] %lld triangles (order %d)",
           fieldName_.c_str(),
           (long long)gCell,
           prm_.cellOrder);
    if (gAmb) {
      printf(", %lld ambiguous", (long long)gAmb);
    }
    if (prm_.cellOrder == 2) {
      printf(", %lld midside nodes", (long long)gMid);
      if (gDeg) {
        printf(", %lld degenerate", (long long)gDeg);
      }
    }
    if (gUns) {
      // NOT a loss.  Before plan H this counted roots the marching mask could
      // not see and the kernel DROPPED -- holes in the surface.  The ring walk
      // keeps both crossings and connects them, so this is now a measure of how
      // folded the field is against the lattice, not of damage.  The oracle
      // asserts every root becomes a point.
      printf(", %lld edge(s) carrying two crossings", (long long)gUns);
    }
    printf("\n");
  }
  reportProjection();
}

// --------------------------------------------------------------------------- //
//  reportProjection -- did niterProj actually converge?
//
//  niterProj is a FIXED step count, not an iteration to a tolerance, so where
//  the surface curves sharply across a sub-cell a node can stop short -- in the
//  wrong place, with nothing in the output to say so.  Printing the residual
//  distribution puts that decision on a measurement.  Collective.
// --------------------------------------------------------------------------- //
inline void isoSurface::reportProjection()
{
  if (prm_.cellOrder != 2) {
    return; // no midside nodes, nothing projected
  }
  MPI_Comm comm = platform->comm.mpiComm();
  std::vector<int> q(7, 0), sh(1, 0);
  o_projResidQ_.copyTo(q, 7);
  o_nProjShort_.copyTo(sh, 1);

  hlong g[7];
  for (int i = 0; i < 7; ++i) {
    g[i] = q[i];
  }
  hlong gShort = sh[0];
  MPI_Allreduce(MPI_IN_PLACE, g, 7, MPI_HLONG, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &gShort, 1, MPI_HLONG, MPI_SUM, comm);
  if (platform->comm.mpiRank() != 0) {
    return;
  }
  hlong tot = 0;
  for (int i = 0; i < 7; ++i) {
    tot += g[i];
  }
  if (tot == 0) {
    return;
  }
  // the worst decade that actually has nodes in it
  int worst = 0;
  for (int i = 6; i >= 0; --i) {
    if (g[i]) {
      worst = i;
      break;
    }
  }
  static const char *edge[7] = {"<=1e-9", "1e-8", "1e-7", "1e-6", "1e-5", "1e-4", ">1e-3"};
  printf("isoSurface(%s): projection niterProj=%d, residual |f|/max|f| over %lld node(s):"
         " worst %s",
         fieldName_.c_str(),
         prm_.niterProj,
         (long long)tot,
         edge[worst]);
  if (gShort) {
    printf(", %lld above 1e-6 (%.2f%%)", (long long)gShort, 100.0 * double(gShort) / double(tot));
  }
  printf("\n");
  if (worst >= 5) {
    printf("isoSurface(%s):   [!] midside nodes are not converging -- raise "
           "IsoSurfaceParams::niterProj (currently %d; 16 clears every demo_3d case)\n",
           fieldName_.c_str(),
           prm_.niterProj);
  }
}

//  uploadPoints -- put ONE level's surface points back into the slabs.
//
//  The mirror of uploadRoots: buildTopology leaves only the LAST level's points
//  on the device, so a multi-level run would interpolate onto the wrong surface.
//  The SurfaceSet is dense; the kernels want the sparse e*capPoint layout.
inline void isoSurface::uploadPoints(const SurfaceSet &ss)
{
  const dlong Nel = mesh_->Nelements;
  const size_t slab = static_cast<size_t>(Nel) * caps_.point;
  std::vector<dfloat> pr(slab, 0), ps(slab, 0), pt(slab, 0);
  std::vector<int> cnt(Nel, 0);

  for (dlong i = 0; i < ss.nPoint; ++i) {
    const dlong e = ss.pelem[i];
    const int m = cnt[e]++;
    const size_t o = static_cast<size_t>(e) * caps_.point + m;
    pr[o] = ss.pr[i];
    ps[o] = ss.ps[i];
    pt[o] = ss.pt[i];
  }
  o_ptR_.copyFrom(pr.data(), pr.size());
  o_ptS_.copyFrom(ps.data(), ps.size());
  o_ptT_.copyFrom(pt.data(), pt.size());
  // points are stored element-major in the SurfaceSet, so cnt[] is exactly the
  // per-element total the interp kernel must loop over
  ptTotal_.assign(cnt.begin(), cnt.end());
  o_ptTotal_.copyFrom(ptTotal_.data(), ptTotal_.size());
}

// =========================================================================== //
//  interpFields -- STEP 3 of plan F: sampling runs on the DEVICE.
//
//  The host prologue stays: pulling the registered device fields into one
//  contiguous buffer, and concatenating any host-registered ones.  That is
//  bookkeeping over a handful of arrays, not per-point work, and there is no
//  kernel for it.
//
//  The per-point work -- one set of basis weights, reused across every field
//  and component -- is isoInterpOKL.  The weights are the expensive part and do
//  not depend on the field, which is why the kernel loops components inside.
//
//  NOTE it is passed o_ptTotal_, NOT o_ptCount_.  The topology kernel leaves
//  CORNERS ONLY in ptCount; at cellOrder 2 the midside nodes sit after them, and
//  passing ptCount would silently interpolate nothing there -- a .vtu with zeros
//  in half its field array, which looks plausible in ParaView.
// =========================================================================== //
inline void isoSurface::interpFields()
{
  platform->timer.tic("isoSurface::interp");

  const int nDev = static_cast<int>(outFields_.size());
  const int nHost = static_cast<int>(hostName_.size());
  const int nF = nDev + nHost;

  const dlong Nel = mesh_->Nelements;
  const int nLev = static_cast<int>(prm_.levels.size());

  // ---- host prologue: gather every registered field, device and host -------
  std::vector<std::string> fName(nF);
  std::vector<int> fComp(nF);
  outHost_.clear();
  outHost_.resize(nF);
  for (int f = 0; f < nDev; ++f) {
    const auto &comps = std::get<1>(outFields_[f]);
    const int nc = static_cast<int>(comps.size());
    fName[f] = std::get<0>(outFields_[f]);
    fComp[f] = nc;
    outHost_[f].assign(static_cast<size_t>(nc) * mesh_->Nlocal, 0.0);
    std::vector<dfloat> tmp(mesh_->Nlocal);
    for (int c = 0; c < nc; ++c) {
      comps[c].copyTo(tmp, static_cast<size_t>(mesh_->Nlocal));
      double *dst = outHost_[f].data() + static_cast<size_t>(c) * mesh_->Nlocal;
      for (dlong i = 0; i < mesh_->Nlocal; ++i) {
        dst[i] = static_cast<double>(tmp[i]);
      }
    }
  }
  for (int h = 0; h < nHost; ++h) {
    const int f = nDev + h;
    fName[f] = hostName_[h];
    fComp[f] = hostComp_[h];
    outHost_[f] = hostVal_[h];
  }

  const size_t slabP = static_cast<size_t>(Nel) * caps_.point;

  for (int l = 0; l < nLev; ++l) {
    const SurfaceSet &ss = surfSets_[l];
    FieldSet &fs = fieldSets_[l];
    fs.clear();
    fs.nPoint = ss.nPoint;
    for (int f = 0; f < nF; ++f) {
      fs.name.push_back(fName[f]);
      fs.nComp.push_back(fComp[f]);
      fs.val.emplace_back(static_cast<size_t>(fComp[f]) * ss.nPoint, static_cast<dfloat>(0));
    }
    if (nF == 0 || ss.nPoint == 0) {
      continue;
    }

    // the surface points for THIS level have to be on the device: buildTopology
    // leaves only the last level's there, exactly as solveRoots does for roots.
    uploadPoints(ss);

    for (int f = 0; f < nF; ++f) {
      const int nc = fComp[f];

      // the source field, host -> device, component-major
      std::vector<dfloat> src(static_cast<size_t>(nc) * mesh_->Nlocal);
      for (size_t q = 0; q < src.size(); ++q) {
        src[q] = static_cast<dfloat>(outHost_[f][q]);
      }
      o_src_.resize(src.size());
      o_src_.copyFrom(src.data(), src.size());
      o_fldVal_.resize(slabP * nc);

      isoInterpOKL(Nel,
                   caps_.point,
                   static_cast<dlong>(mesh_->Nlocal),
                   static_cast<int>(nc),
                   o_xn_,
                   o_bw_,
                   o_src_,
                   o_ptTotal_, // corners AND midside -- see the note above
                   o_ptR_,
                   o_ptS_,
                   o_ptT_,
                   o_fldVal_);

      // ---- back to the dense FieldSet the writers consume ------------------
      std::vector<dfloat> val(slabP * nc);
      o_fldVal_.copyTo(val, val.size());
      dlong dense = 0;
      for (dlong e = 0; e < Nel; ++e) {
        for (int m = 0; m < ptTotal_[e]; ++m) {
          const size_t src_o = (static_cast<size_t>(e) * caps_.point + m) * nc;
          for (int c = 0; c < nc; ++c) {
            fs.val[f][static_cast<size_t>(dense) * nc + c] = val[src_o + c];
          }
          ++dense;
        }
      }
      stats_.nInterpPoint += dense;
    }
  }

  platform->timer.toc("isoSurface::interp");

  hlong gPts = 0;
  for (int l = 0; l < nLev; ++l) {
    gPts += fieldSets_[l].nPoint;
  }
  MPI_Allreduce(MPI_IN_PLACE, &gPts, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm());
  if (platform->comm.mpiRank() == 0 && nF > 0) {
    printf("isoSurface(%s): [GPU] interpolated %d field(s) onto %lld surface points\n",
           fieldName_.c_str(),
           nF,
           (long long)gPts);
  }
}

inline void isoSurface::ensureGeometry()
{
  ensureGeometry(prm_.ifSolveGeom);
}

// solveGeom == true forces a re-solve; false reuses a surface already built.
// Note that `false` is NOT "use the old field": the first call still solves,
// because rootsValid_ is false until one has.
inline void isoSurface::ensureGeometry(bool solveGeom)
{
  if (solveGeom) {
    solveRoots(true);
    buildTopology(true);
  } else {
    solveRoots();
    buildTopology();
  }
}

//  process -- the deliverable, and only the deliverable.
//
//  Geometry per prm.ifSolveGeom, fields sampled onto it, one _fld_ file per
//  isolevel.  A caller wanting the debug products calls writeRoots() or
//  writeSurface() as well; there are no dump flags.  dumpCounter_ advances HERE
//  ONLY, so a debug write shares the dump index of the fields file it belongs to.
inline void isoSurface::process(double time, int tstep)
{
  process(time, tstep, prm_.ifSolveGeom);
}

inline void isoSurface::process(double time, int tstep, bool solveGeom)
{
  ensureGeometry(solveGeom);
  interpFields();
  writeFields(time, tstep);
  ++dumpCounter_;
}

inline std::string isoSurface::outName(int level, int counter, Product what) const
{
  const char *tag = (what == Srf) ? "_srf_" : (what == Fld) ? "_fld_" : "_roots_";
  std::ostringstream o;
  o << prm_.fileName << tag << std::setw(2) << std::setfill('0') << level << "_" << std::setw(5)
    << std::setfill('0') << counter << ".vtu";
  return o.str();
}

//  writeVtm -- one .vtm per dump and product, grouping that dump's levels.
//
//  Rank 0 only, after every level's collective MPI_File_close, so no block names a
//  missing .vtu.  Block index == level index; the name carries the level value.
//  Written only for more than one level.  Members are basenames (same directory).
//  The XML itself is a pure function so test_vtm.sh can extract and test it.

// BEGIN isoVtmXml
inline std::string isoVtmXml(const std::vector<std::string> &names, const std::vector<std::string> &files)
{
  std::ostringstream x;
  x << "<?xml version=\"1.0\"?>\n"
    << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
    << "  <vtkMultiBlockDataSet>\n";
  for (size_t l = 0; l < files.size(); ++l) {
    // names hold fieldName_ and a %g number: no '&', '<' or '"' can occur.
    x << "    <DataSet index=\"" << l << "\" name=\"" << names[l] << "\" file=\"" << files[l] << "\"/>\n";
  }
  x << "  </vtkMultiBlockDataSet>\n</VTKFile>\n";
  return x.str();
}
// END isoVtmXml

inline void isoSurface::writeVtm(Product what, int counter) const
{
  const int nLev = static_cast<int>(prm_.levels.size());
  if (nLev < 2 || platform->comm.mpiRank() != 0) {
    return;
  }
  const char *tag = (what == Srf) ? "_srf_" : (what == Fld) ? "_fld_" : "_roots_";
  std::ostringstream fn;
  fn << prm_.fileName << tag << std::setw(5) << std::setfill('0') << counter << ".vtm";

  std::vector<std::string> names, files;
  for (int l = 0; l < nLev; ++l) {
    std::string member = outName(l, counter, what);
    const size_t slash = member.find_last_of('/');
    files.push_back(slash == std::string::npos ? member : member.substr(slash + 1));
    char nm[160];
    snprintf(nm, sizeof nm, "level_%02d %s=%.9g", l, fieldName_.c_str(), static_cast<double>(prm_.levels[l]));
    names.push_back(nm);
  }
  std::ofstream f(fn.str());
  f << isoVtmXml(names, files);
  printf("isoSurface(%s): wrote %s  (%d levels)\n", fieldName_.c_str(), fn.str().c_str(), nLev);
}

// =========================================================================== //
//  writeRoots -- one .vtu per isolevel, collective MPI-IO.
//
//  The file is an UnstructuredGrid of VTK_VERTEX cells: a point cloud that opens
//  directly in ParaView and carries, per point, everything milestone 2+ needs to
//  rebuild a curved surface WITHOUT re-reading the element data:
//
//      Points     physical x,y,z
//      rst        reference coordinates -- lets a later stage re-enter the element
//      gradPhys   grad_x f ; the unit normal is gradPhys/|gradPhys|
//      gradRef    grad_r f ; the Hermite tangent along a fibre
//      eg         GLOBAL element id.  The rank-local `elem` is deliberately NOT
//                 written: it is meaningless once the per-rank output is merged
//                 into one file, and `eg` is what element_table() groups by.
//      dirjk      fibre direction and its two node indices
//      flag       tangency / multi-edge / on-face bits
//
//  The isolevel, time and timestep are FieldData, so a reader never has to parse
//  the filename (the fieldExtract convention).
//
//  Layout follows fieldExtract::writeVts: an ASCII XML header written by rank 0,
//  then raw appended binary blocks written collectively with MPI_File_write_at_all,
//  each prefixed by a UInt64 byte count (VTK header_type="UInt64").
// =========================================================================== //
inline void
isoSurface::writeRootVtu(const RootSet &rs, int lev, const std::string &fname, double time, int tstep)
{
  MPI_Comm comm = platform->comm.mpiComm();
  const int rank = platform->comm.mpiRank();

  // ---- global point count and this rank's offset --------------------------
  const long long nLocal = static_cast<long long>(rs.nRoots);
  long long nElemWithRoots = 0;
  for (dlong e = 0; e < static_cast<dlong>(rs.eCount.size()); ++e) {
    if (rs.eCount[e] > 0) {
      ++nElemWithRoots;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &nElemWithRoots, 1, MPI_LONG_LONG, MPI_SUM, comm);
  long long nGlobal = 0, myOffset = 0;
  MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, comm);
  MPI_Exscan(&nLocal, &myOffset, 1, MPI_LONG_LONG, MPI_SUM, comm);
  if (rank == 0) {
    myOffset = 0;
  }

  // ---- pack the payload.  float32 for geometry and gradients (plenty: these
  //      are consumed by a visualiser), int32/uint8 for the indices and flags.
  std::vector<float> pos(3 * nLocal), rst(3 * nLocal), gP(3 * nLocal), gR(3 * nLocal);
  std::vector<int> idirjk(3 * nLocal);
  std::vector<long long> ieg(nLocal);
  std::vector<int> conn(nLocal), offs(nLocal);
  std::vector<unsigned char> types(nLocal), flg(nLocal);
  for (long long i = 0; i < nLocal; ++i) {
    pos[3 * i + 0] = static_cast<float>(rs.x[i]);
    pos[3 * i + 1] = static_cast<float>(rs.y[i]);
    pos[3 * i + 2] = static_cast<float>(rs.z[i]);
    rst[3 * i + 0] = static_cast<float>(rs.r[i]);
    rst[3 * i + 1] = static_cast<float>(rs.s[i]);
    rst[3 * i + 2] = static_cast<float>(rs.t[i]);
    gP[3 * i + 0] = static_cast<float>(rs.dfdx[i]);
    gP[3 * i + 1] = static_cast<float>(rs.dfdy[i]);
    gP[3 * i + 2] = static_cast<float>(rs.dfdz[i]);
    gR[3 * i + 0] = static_cast<float>(rs.dfdr[i]);
    gR[3 * i + 1] = static_cast<float>(rs.dfds[i]);
    gR[3 * i + 2] = static_cast<float>(rs.dfdt[i]);
    ieg[i] = static_cast<long long>(rs.eg[i]);
    idirjk[3 * i + 0] = static_cast<int>(rs.dir[i]);
    idirjk[3 * i + 1] = static_cast<int>(rs.j[i]);
    idirjk[3 * i + 2] = static_cast<int>(rs.k[i]);
    conn[i] = static_cast<int>(myOffset + i); // one VTK_VERTEX per point
    offs[i] = static_cast<int>(myOffset + i + 1);
    types[i] = 1; // VTK_VERTEX
    flg[i] = rs.flag[i];
  }

  // ---- appended-block bookkeeping -----------------------------------------
  // Each block is [UInt64 nbytes][payload]; offsets are relative to the byte
  // after the "_" that opens <AppendedData>.
  struct Blk {
    const char *name;
    const char *type;
    int nComp;
    int itemBytes;
    const void *data;
    MPI_Datatype mpiType;
    bool isCell;
  };

  const std::vector<Blk> blk = {{"Position", "Float32", 3, 4, pos.data(), MPI_FLOAT, false},
                                {"rst", "Float32", 3, 4, rst.data(), MPI_FLOAT, false},
                                {"gradPhys", "Float32", 3, 4, gP.data(), MPI_FLOAT, false},
                                {"gradRef", "Float32", 3, 4, gR.data(), MPI_FLOAT, false},
                                {"eg", "Int64", 1, 8, ieg.data(), MPI_LONG_LONG, false},
                                {"dirjk", "Int32", 3, 4, idirjk.data(), MPI_INT, false},
                                {"flag", "UInt8", 1, 1, flg.data(), MPI_UNSIGNED_CHAR, false},
                                {"connectivity", "Int32", 1, 4, conn.data(), MPI_INT, true},
                                {"offsets", "Int32", 1, 4, offs.data(), MPI_INT, true},
                                {"types", "UInt8", 1, 1, types.data(), MPI_UNSIGNED_CHAR, true}};

  std::vector<long long> blkOffset(blk.size(), 0);
  long long running = 0;
  for (size_t b = 0; b < blk.size(); ++b) {
    blkOffset[b] = running;
    running += 8 + static_cast<long long>(blk[b].nComp) * blk[b].itemBytes * nGlobal;
  }
  const long long totalAppended = running;

  // ---- XML header ----------------------------------------------------------
  auto ind = [](int n) { return std::string(2 * n, ' '); };
  auto darray = [&](const Blk &b, long long off) {
    return "<DataArray type=\"" + std::string(b.type) + "\" Name=\"" + b.name + "\" NumberOfComponents=\"" +
           std::to_string(b.nComp) + "\" format=\"appended\" offset=\"" + std::to_string(off) + "\"/>\n";
  };

  std::string msg;
  msg += "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
         "header_type=\"UInt64\">\n";
  msg += ind(1) + "<UnstructuredGrid>\n";

  // metadata: the isolevel lives here, not in the filename
  msg += ind(2) + "<FieldData>\n";
  msg += ind(3) +
         "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(time) + " </DataArray>\n";
  msg += ind(3) +
         "<DataArray type=\"Int32\" Name=\"CYCLE\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(tstep) + " </DataArray>\n";
  msg += ind(3) +
         "<DataArray type=\"Float64\" Name=\"isoLevel\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(static_cast<double>(rs.level)) + " </DataArray>\n";
  msg += ind(3) +
         "<DataArray type=\"Int32\" Name=\"isoLevelIndex\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(lev) + " </DataArray>\n";
  msg += ind(3) +
         "<DataArray type=\"Int32\" Name=\"polynomialOrder\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(N_) + " </DataArray>\n";
  msg += ind(3) +
         "<DataArray type=\"Int64\" Name=\"numPoints\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(nGlobal) + " </DataArray>\n";
  msg += ind(3) +
         "<DataArray type=\"Float64\" Name=\"cullTheta\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(static_cast<double>(prm_.cullTheta)) + " </DataArray>\n";
  // how many elements survived culling and actually carry a crossing, against the
  // whole mesh: the pair is the cull rate, recoverable from the file alone.
  msg += ind(3) +
         "<DataArray type=\"Int64\" Name=\"numElementsWithRoots\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(nElemWithRoots) + " </DataArray>\n";
  msg += ind(3) +
         "<DataArray type=\"Int64\" Name=\"numElementsGlobal\" NumberOfTuples=\"1\" "
         "format=\"ascii\"> " +
         std::to_string(static_cast<long long>(mesh_->NelementsGlobal)) + " </DataArray>\n";
  msg += ind(2) + "</FieldData>\n";

  msg += ind(2) + "<Piece NumberOfPoints=\"" + std::to_string(nGlobal) + "\" NumberOfCells=\"" +
         std::to_string(nGlobal) + "\">\n";
  msg += ind(3) + "<Points>\n";
  msg += ind(4) + darray(blk[0], blkOffset[0]);
  msg += ind(3) + "</Points>\n";
  msg += ind(3) + "<PointData>\n";
  for (size_t b = 1; b <= 6; ++b) {
    msg += ind(4) + darray(blk[b], blkOffset[b]);
  }
  msg += ind(3) + "</PointData>\n";
  msg += ind(3) + "<Cells>\n";
  for (size_t b = 7; b < blk.size(); ++b) {
    msg += ind(4) + darray(blk[b], blkOffset[b]);
  }
  msg += ind(3) + "</Cells>\n";
  msg += ind(2) + "</Piece>\n";
  msg += ind(1) + "</UnstructuredGrid>\n";
  msg += ind(1) + "<AppendedData encoding=\"raw\">\n";
  msg += "_";

  const MPI_Offset appendedStart = static_cast<MPI_Offset>(msg.size());

  MPI_File fh;
  const int err = MPI_File_open(comm, fname.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  nekrsCheck(err != MPI_SUCCESS,
             comm,
             EXIT_FAILURE,
             "isoSurface: cannot open %s for writing\n",
             fname.c_str());
  MPI_File_set_size(fh, 0);

  if (rank == 0) {
    MPI_File_write_at(fh, 0, msg.c_str(), msg.size(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  for (size_t b = 0; b < blk.size(); ++b) {
    const long long itemsLocal = static_cast<long long>(blk[b].nComp) * nLocal;
    const MPI_Offset blockStart = appendedStart + blkOffset[b];
    if (rank == 0) {
      const unsigned long long nbytes =
          static_cast<unsigned long long>(blk[b].nComp) * blk[b].itemBytes * nGlobal;
      MPI_File_write_at(fh, blockStart, &nbytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_STATUS_IGNORE);
    }
    const MPI_Offset dataStart =
        blockStart + 8 + static_cast<MPI_Offset>(blk[b].nComp) * blk[b].itemBytes * myOffset;
    MPI_File_write_at_all(fh,
                          dataStart,
                          const_cast<void *>(blk[b].data),
                          static_cast<int>(itemsLocal),
                          blk[b].mpiType,
                          MPI_STATUS_IGNORE);
  }

  const std::string footer = std::string("\n") + ind(1) + "</AppendedData>\n</VTKFile>\n";
  if (rank == 0) {
    MPI_File_write_at(fh,
                      appendedStart + static_cast<MPI_Offset>(totalAppended),
                      footer.c_str(),
                      footer.size(),
                      MPI_CHAR,
                      MPI_STATUS_IGNORE);
  }
  MPI_File_close(&fh);
}

inline void isoSurface::writeRoots(double time, int tstep)
{
  platform->timer.tic("isoSurface::writeRoots");
  for (size_t l = 0; l < rootSets_.size(); ++l) {
    const std::string fn = outName(static_cast<int>(l), dumpCounter_, Roots);
    writeRootVtu(rootSets_[l], static_cast<int>(l), fn, time, tstep);
    if (platform->comm.mpiRank() == 0) {
      long long n = 0;
      const long long loc = rootSets_[l].nRoots;
      MPI_Reduce(const_cast<long long *>(&loc), &n, 1, MPI_LONG_LONG, MPI_SUM, 0, platform->comm.mpiComm());
      printf("isoSurface(%s): wrote %s  level=%g  points=%lld\n",
             fieldName_.c_str(),
             fn.c_str(),
             static_cast<double>(rootSets_[l].level),
             n);
    } else {
      const long long loc = rootSets_[l].nRoots;
      MPI_Reduce(const_cast<long long *>(&loc),
                 nullptr,
                 1,
                 MPI_LONG_LONG,
                 MPI_SUM,
                 0,
                 platform->comm.mpiComm());
    }
  }
  writeVtm(Roots, dumpCounter_);
  platform->timer.toc("isoSurface::writeRoots");
}

inline void isoSurface::setOutputFields(const std::vector<field> &flds)
{
  outFields_ = flds;

  // NOTE the asymmetry with setField(): that one hard-errors on more than one
  // component, because the EXTRACTION scalar must be a scalar.  Here a vector is
  // the common case -- {"velocity", {o_ux, o_uy, o_uz}} -- so only the empty
  // case is rejected.
  for (const auto &f : outFields_) {
    nekrsCheck(std::get<1>(f).empty(),
               platform->comm.mpiComm(),
               EXIT_FAILURE,
               "isoSurface: output field '%s' has no components\n",
               std::get<0>(f).c_str());
  }

  if (platform->comm.mpiRank() == 0 && !outFields_.empty()) {
    printf("isoSurface(%s): %d output field(s):", fieldName_.c_str(), (int)outFields_.size());
    for (const auto &f : outFields_) {
      printf(" %s(%d)", std::get<0>(f).c_str(), (int)std::get<1>(f).size());
    }
    printf("\n");
  }
}

inline void
isoSurface::setOutputFieldHost(const std::string &name, int nComp, const std::vector<dfloat> &vals)
{
  nekrsCheck(nComp < 1,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "isoSurface: output field '%s' needs >= 1 component\n",
             name.c_str());
  nekrsCheck(static_cast<dlong>(vals.size()) < nComp * mesh_->Nlocal,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "isoSurface: output field '%s' has %d values, expected %d\n",
             name.c_str(),
             (int)vals.size(),
             (int)(nComp * mesh_->Nlocal));
  hostName_.push_back(name);
  hostComp_.push_back(nComp);
  hostVal_.emplace_back(static_cast<size_t>(nComp) * mesh_->Nlocal);
  for (size_t q = 0; q < hostVal_.back().size(); ++q) {
    hostVal_.back()[q] = static_cast<double>(vals[q]);
  }
}

// (Re)build roots + topology according to prm.ifSolveGeom.  Shared by process()

// =========================================================================== //
//  writeSurface -- the triangulated surface, one .vtu per isolevel.
//
//  Same MPI-IO scaffolding as writeRootVtu (ASCII header from rank 0, appended
//  raw blocks written collectively, UInt64 block prefixes); what differs is the
//  cell block, which now carries real triangles instead of one VTK_VERTEX per
//  point.
//
//  Cells are VTK_TRIANGLE (5) at cellOrder 1 and VTK_QUADRATIC_TRIANGLE (22) at
//  cellOrder 2, with node order 0,1,2 = corners and 3,4,5 = midside of (0,1),
//  (1,2), (2,0) -- read out of VTK itself, not from documentation.  `offsets` is
//  the CUMULATIVE END index into connectivity, one entry per cell.
//
//  cellOrder is recorded as FieldData so a reader takes it from the file rather
//  than inferring it from the types array.
// =========================================================================== //
inline void isoSurface::writeSurfVtu(const SurfaceSet &ss,
                                     const FieldSet *fs,
                                     int lev,
                                     const std::string &fname,
                                     double time,
                                     int tstep)
{
  MPI_Comm comm = platform->comm.mpiComm();
  const int rank = platform->comm.mpiRank();
  const int npc = ss.nodesPerCell;

  const long long nPtLocal = static_cast<long long>(ss.nPoint);
  const long long nClLocal = static_cast<long long>(ss.nCell);
  long long nPtGlobal = 0, ptOffset = 0, nClGlobal = 0, clOffset = 0;
  MPI_Allreduce(&nPtLocal, &nPtGlobal, 1, MPI_LONG_LONG, MPI_SUM, comm);
  MPI_Allreduce(&nClLocal, &nClGlobal, 1, MPI_LONG_LONG, MPI_SUM, comm);
  MPI_Exscan(&nPtLocal, &ptOffset, 1, MPI_LONG_LONG, MPI_SUM, comm);
  MPI_Exscan(&nClLocal, &clOffset, 1, MPI_LONG_LONG, MPI_SUM, comm);
  if (rank == 0) {
    ptOffset = 0;
    clOffset = 0;
  }

  // ---- pack ----------------------------------------------------------------
  std::vector<float> pos(3 * nPtLocal), nor(3 * nPtLocal);
  std::vector<unsigned char> pflg(nPtLocal);
  for (long long i = 0; i < nPtLocal; ++i) {
    pos[3 * i + 0] = static_cast<float>(ss.px[i]);
    pos[3 * i + 1] = static_cast<float>(ss.py[i]);
    pos[3 * i + 2] = static_cast<float>(ss.pz[i]);
    nor[3 * i + 0] = static_cast<float>(ss.nx[i]);
    nor[3 * i + 1] = static_cast<float>(ss.ny[i]);
    nor[3 * i + 2] = static_cast<float>(ss.nz[i]);
    pflg[i] = ss.pflag[i];
  }
  // isCorner: 1 for a point that is a root, 0 for a midside node.  Cheap, and it
  // is what lets a reader recover the root correspondence after the merge.
  std::vector<unsigned char> isCorner(nPtLocal, 0);
  for (long long i = 0; i < static_cast<long long>(ss.nCorner); ++i) {
    isCorner[i] = 1;
  }

  // connectivity is GLOBAL: local point ids shift by this rank's point offset
  std::vector<int> conn(static_cast<size_t>(npc) * nClLocal), offs(nClLocal);
  std::vector<unsigned char> types(nClLocal), cflg(nClLocal);
  std::vector<long long> ceg(nClLocal);
  for (long long c = 0; c < nClLocal; ++c) {
    for (int m = 0; m < npc; ++m) {
      conn[npc * c + m] = static_cast<int>(ptOffset + ss.conn[npc * c + m]);
    }
    offs[c] = static_cast<int>(npc * (clOffset + c + 1)); // cumulative END index
    types[c] = ss.cellType;
    cflg[c] = ss.cellFlag[c];
    ceg[c] = static_cast<long long>(ss.cellEg[c]);
  }

  // A block's XML SECTION is a tag on the block, not its position in the list.
  // It used to be a hardcoded index range (0 / 1..3 / 4..5 / 6..end), which meant
  // appending a per-point field array silently emitted it inside <Cells> -- a
  // corrupt file.  With the tag, fields are appended freely.
  enum Sect { SPoints = 0, SPointData, SCellData, SCells };

  struct Blk {
    const char *name;
    const char *type;
    int nComp;
    int itemBytes;
    const void *data;
    MPI_Datatype mpiType;
    int kind; // 0 = per point, 1 = per cell, 2 = per cell x nodesPerCell
    int sect; // which XML section it belongs in
  };

  // The DIAGNOSTIC blocks -- normal, isCorner, flag, cellFlag -- describe how the
  // surface was built, not what the field does on it.  They belong to
  // writeSurface(), the debug product, and are omitted from the _fld_ file,
  // where they are dead weight: 11 bytes per point and 1 per cell, about a third
  // of a P2 file that carries one vector field.  `fs` is non-null exactly for
  // the field file, so it is the discriminator.
  const bool diag = (fs == nullptr);

  std::vector<Blk> blk = {{"Position", "Float32", 3, 4, pos.data(), MPI_FLOAT, 0, SPoints}};
  if (diag) {
    blk.push_back({"normal", "Float32", 3, 4, nor.data(), MPI_FLOAT, 0, SPointData});
    blk.push_back({"isCorner", "UInt8", 1, 1, isCorner.data(), MPI_UNSIGNED_CHAR, 0, SPointData});
    blk.push_back({"flag", "UInt8", 1, 1, pflg.data(), MPI_UNSIGNED_CHAR, 0, SPointData});
  }
  blk.push_back({"eg", "Int64", 1, 8, ceg.data(), MPI_LONG_LONG, 1, SCellData});
  if (diag) {
    blk.push_back({"cellFlag", "UInt8", 1, 1, cflg.data(), MPI_UNSIGNED_CHAR, 1, SCellData});
  }
  blk.push_back({"connectivity", "Int32", 1, 4, conn.data(), MPI_INT, 2, SCells});
  blk.push_back({"offsets", "Int32", 1, 4, offs.data(), MPI_INT, 1, SCells});
  blk.push_back({"types", "UInt8", 1, 1, types.data(), MPI_UNSIGNED_CHAR, 1, SCells});

  // ---- the interpolated fields, as PointData ------------------------------
  // Packed here (not in the FieldSet) so the narrowing to Float32 and the
  // lifetime of the buffers stay local to the write, as for every other block.
  std::vector<std::vector<float>> fpack32;
  std::vector<std::vector<double>> fpack64;
  if (fs && fs->nField() > 0) {
    nekrsCheck(fs->nPoint != ss.nPoint,
               comm,
               EXIT_FAILURE,
               "isoSurface: field set has %lld points but the surface has %lld\n",
               (long long)fs->nPoint,
               (long long)ss.nPoint);
    fpack32.resize(fs->nField());
    fpack64.resize(fs->nField());
    for (int f = 0; f < fs->nField(); ++f) {
      const int nc = fs->nComp[f];
      const size_t n = static_cast<size_t>(nc) * nPtLocal;
      if (prm_.fieldsFloat64) {
        fpack64[f].resize(n);
        for (size_t q = 0; q < n; ++q) {
          fpack64[f][q] = static_cast<double>(fs->val[f][q]);
        }
        blk.push_back({fs->name[f].c_str(), "Float64", nc, 8, fpack64[f].data(), MPI_DOUBLE, 0, SPointData});
      } else {
        fpack32[f].resize(n);
        for (size_t q = 0; q < n; ++q) {
          fpack32[f][q] = static_cast<float>(fs->val[f][q]);
        }
        blk.push_back({fs->name[f].c_str(), "Float32", nc, 4, fpack32[f].data(), MPI_FLOAT, 0, SPointData});
      }
    }
  }

  auto globalCount = [&](int kind) -> long long {
    return (kind == 0) ? nPtGlobal : (kind == 1) ? nClGlobal : nClGlobal * npc;
  };
  auto localCount = [&](int kind) -> long long {
    return (kind == 0) ? nPtLocal : (kind == 1) ? nClLocal : nClLocal * npc;
  };
  auto localOffset = [&](int kind) -> long long {
    return (kind == 0) ? ptOffset : (kind == 1) ? clOffset : clOffset * npc;
  };

  std::vector<long long> blkOffset(blk.size(), 0);
  long long running = 0;
  for (size_t b = 0; b < blk.size(); ++b) {
    blkOffset[b] = running;
    running += 8 + static_cast<long long>(blk[b].nComp) * blk[b].itemBytes * globalCount(blk[b].kind);
  }
  const long long totalAppended = running;

  auto ind = [](int n) { return std::string(2 * n, ' '); };
  auto darray = [&](const Blk &b, long long off) {
    return "<DataArray type=\"" + std::string(b.type) + "\" Name=\"" + b.name + "\" NumberOfComponents=\"" +
           std::to_string(b.nComp) + "\" format=\"appended\" offset=\"" + std::to_string(off) + "\"/>\n";
  };
  auto fdScalar = [&](const char *type, const char *name, const std::string &val) {
    return ind(3) + "<DataArray type=\"" + type + "\" Name=\"" + name +
           "\" NumberOfTuples=\"1\" format=\"ascii\"> " + val + " </DataArray>\n";
  };

  std::string msg;
  msg += "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
         "header_type=\"UInt64\">\n";
  msg += ind(1) + "<UnstructuredGrid>\n";
  msg += ind(2) + "<FieldData>\n";
  msg += fdScalar("Float64", "TimeValue", std::to_string(time));
  msg += fdScalar("Int32", "CYCLE", std::to_string(tstep));
  msg += fdScalar("Float64", "isoLevel", std::to_string(static_cast<double>(ss.level)));
  msg += fdScalar("Int32", "isoLevelIndex", std::to_string(lev));
  msg += fdScalar("Int32", "polynomialOrder", std::to_string(N_));
  msg += fdScalar("Int32", "cellOrder", std::to_string(ss.cellOrder));
  msg += fdScalar("Int32", "cellType", std::to_string(static_cast<int>(ss.cellType)));
  msg += fdScalar("Int32", "nodesPerCell", std::to_string(npc));
  msg += fdScalar("Int64", "numPoints", std::to_string(nPtGlobal));
  msg += fdScalar("Int64", "numCells", std::to_string(nClGlobal));
  msg += ind(2) + "</FieldData>\n";

  msg += ind(2) + "<Piece NumberOfPoints=\"" + std::to_string(nPtGlobal) + "\" NumberOfCells=\"" +
         std::to_string(nClGlobal) + "\">\n";
  // emit by SECTION TAG, so a block's position in the list no longer decides
  // which XML element it lands in
  auto section = [&](int want, const char *tag) {
    std::string out = ind(3) + "<" + tag + ">\n";
    for (size_t b = 0; b < blk.size(); ++b) {
      if (blk[b].sect == want) {
        out += ind(4) + darray(blk[b], blkOffset[b]);
      }
    }
    return out + ind(3) + "</" + tag + ">\n";
  };
  msg += section(SPoints, "Points");
  msg += section(SPointData, "PointData");
  msg += section(SCellData, "CellData");
  msg += section(SCells, "Cells");
  msg += ind(2) + "</Piece>\n";
  msg += ind(1) + "</UnstructuredGrid>\n";
  msg += ind(1) + "<AppendedData encoding=\"raw\">\n";
  msg += "_";

  const MPI_Offset appendedStart = static_cast<MPI_Offset>(msg.size());

  MPI_File fh;
  const int err = MPI_File_open(comm, fname.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  nekrsCheck(err != MPI_SUCCESS,
             comm,
             EXIT_FAILURE,
             "isoSurface: cannot open %s for writing\n",
             fname.c_str());
  MPI_File_set_size(fh, 0);

  if (rank == 0) {
    MPI_File_write_at(fh, 0, msg.c_str(), msg.size(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  for (size_t b = 0; b < blk.size(); ++b) {
    const MPI_Offset blockStart = appendedStart + blkOffset[b];
    if (rank == 0) {
      const unsigned long long nbytes =
          static_cast<unsigned long long>(blk[b].nComp) * blk[b].itemBytes * globalCount(blk[b].kind);
      MPI_File_write_at(fh, blockStart, &nbytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_STATUS_IGNORE);
    }
    const MPI_Offset dataStart =
        blockStart + 8 + static_cast<MPI_Offset>(blk[b].nComp) * blk[b].itemBytes * localOffset(blk[b].kind);
    MPI_File_write_at_all(fh,
                          dataStart,
                          const_cast<void *>(blk[b].data),
                          static_cast<int>(static_cast<long long>(blk[b].nComp) * localCount(blk[b].kind)),
                          blk[b].mpiType,
                          MPI_STATUS_IGNORE);
  }

  const std::string footer = std::string("\n") + ind(1) + "</AppendedData>\n</VTKFile>\n";
  if (rank == 0) {
    MPI_File_write_at(fh,
                      appendedStart + static_cast<MPI_Offset>(totalAppended),
                      footer.c_str(),
                      footer.size(),
                      MPI_CHAR,
                      MPI_STATUS_IGNORE);
  }
  MPI_File_close(&fh);
}

inline void isoSurface::writeSurface(double time, int tstep)
{
  if (!surfValid_) {
    buildTopology();
  }
  platform->timer.tic("isoSurface::writeSurface");
  for (size_t l = 0; l < surfSets_.size(); ++l) {
    const std::string fn = outName(static_cast<int>(l), dumpCounter_, Srf);
    writeSurfVtu(surfSets_[l], nullptr, static_cast<int>(l), fn, time, tstep);

    const long long loc[2] = {static_cast<long long>(surfSets_[l].nPoint),
                              static_cast<long long>(surfSets_[l].nCell)};
    long long tot[2] = {0, 0};
    MPI_Reduce(const_cast<long long *>(loc), tot, 2, MPI_LONG_LONG, MPI_SUM, 0, platform->comm.mpiComm());
    if (platform->comm.mpiRank() == 0) {
      printf("isoSurface(%s): wrote %s  level=%g  points=%lld  cells=%lld (order %d)\n",
             fieldName_.c_str(),
             fn.c_str(),
             static_cast<double>(surfSets_[l].level),
             tot[0],
             tot[1],
             surfSets_[l].cellOrder);
    }
  }
  writeVtm(Srf, dumpCounter_);
  platform->timer.toc("isoSurface::writeSurface");
}

// --------------------------------------------------------------------------- //
//  writeFields -- THE DELIVERABLE.  The surface, carrying the interpolated
//  fields as PointData, one file per isolevel.  Same writer as writeSurface,
//  with N more per-point blocks.

inline void isoSurface::writeFields(double time, int tstep)
{
  // LAZY, never forced: solveRoots()/buildTopology() are no-ops when valid, so a
  // standalone writeFields() still works, while a call from process() -- which
  // has already applied the ifSolveGeom policy -- does not solve a second time.
  solveRoots();
  buildTopology();
  if (fieldSets_.empty() || fieldSets_[0].nPoint != surfSets_[0].nPoint) {
    interpFields(); // not sampled yet, or the surface moved under it
  }

  platform->timer.tic("isoSurface::writeFields");
  for (size_t l = 0; l < surfSets_.size(); ++l) {
    const std::string fn = outName(static_cast<int>(l), dumpCounter_, Fld);
    writeSurfVtu(surfSets_[l], &fieldSets_[l], static_cast<int>(l), fn, time, tstep);

    const long long loc[2] = {static_cast<long long>(surfSets_[l].nPoint),
                              static_cast<long long>(surfSets_[l].nCell)};
    long long tot[2] = {0, 0};
    MPI_Reduce(const_cast<long long *>(loc), tot, 2, MPI_LONG_LONG, MPI_SUM, 0, platform->comm.mpiComm());
    if (platform->comm.mpiRank() == 0) {
      printf("isoSurface(%s): wrote %s  level=%g  points=%lld  cells=%lld  fields=%d (order %d)\n",
             fieldName_.c_str(),
             fn.c_str(),
             static_cast<double>(surfSets_[l].level),
             tot[0],
             tot[1],
             fieldSets_[l].nField(),
             surfSets_[l].cellOrder);
    }
  }
  writeVtm(Fld, dumpCounter_);
  platform->timer.toc("isoSurface::writeFields");
}

// Milestone 4.  Named now so the call site is stable; interpolation onto the
// roots is independent of the solve, which is why it is a separate verb.
// =========================================================================== //
//  MILESTONE 5 -- field interpolation onto the surface.

#endif // nekrs_isosurface_hpp_
