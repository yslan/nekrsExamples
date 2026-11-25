"""
route_c.py -- ROUTE C: Ehrlich-Aberth iteration on the Legendre series (MATRIX-FREE).

    GLL nodal values -> Legendre coefficients -> deflate -> simultaneous root iteration

No matrix is ever formed.  Designed for the 3D/GPU target: ONE FIBER PER THREAD,
fixed iteration count, no branches in the hot loop, ~3N doubles of state.

------------------------------------------------------------------------------
MATH
------------------------------------------------------------------------------
Steps 1-2 are shared with route A (see demo_A/route_a.py): project the GLL nodal
values onto the Legendre basis,

    p(x) = sum_{k=0}^{n} a_k P_k(x),     a = M f,   gamma_N = 2/N (not 2/(2N+1))

and deflate trailing negligible coefficients so n is the polynomial's TRUE degree.
Route A then builds an n x n comrade matrix and calls an eigensolver.  Route C does
not build a matrix at all.

EHRLICH-ABERTH.  Given current estimates z_1..z_n of ALL n roots simultaneously,
each is updated by a Newton step corrected for the pull of the other roots:

                       w_i
    z_i  <-  z_i  -  ---------- ,        w_i = p(z_i) / p'(z_i)
                     1 - w_i S_i

                     S_i = sum_{j != i}  1 / (z_i - z_j)

Interpretation: Newton's method applied to p(x) / prod_{j != i}(x - z_j), i.e. to the
polynomial with the OTHER roots divided out.  The S_i term is what stops two estimates
collapsing onto the same root -- the property a sequence of independent Newton solves
lacks, and the reason this finds ALL roots rather than the same one n times.
Convergence is CUBIC for simple roots.

EVALUATING p AND p' WITHOUT A MATRIX.  The same Sturm-Liouville three-term recurrence
that builds the comrade matrix, run forwards as an evaluation:

    P_0 = 1,  P_1 = x,    P_{k+1} = ((2k+1) x P_k - k P_{k-1}) / (k+1)
    P'_0 = 0, P'_1 = 1,   P'_{k+1} = ((2k+1)(P_k + x P'_k) - k P'_{k-1}) / (k+1)

accumulating p = sum a_k P_k and p' = sum a_k P'_k as it goes.  O(n) work, O(1)
storage, no allocation.  Verified against numpy's legval/legder: value error 8.9e-16,
derivative error 3.6e-15.

INITIALIZATION (the one real trap).  The textbook monomial initialization -- n points
on a circle of radius ~ the Cauchy bound -- DIVERGES here: measured corrections grew
from 1e+00 to 8e+01 over 12 iterations.  The roots of interest lie in [-1,1], so start
on Chebyshev points nudged off the real axis:

    z_i = cos(pi (i+1/2)/n)  +  0.1i * sin(pi (i+1/2)/n),    i = 0..n-1

The imaginary part matters: with all estimates exactly real and collinear the repulsion
term S_i stays real, so the iteration is confined to the real axis and cannot reach a
complex root (these exist routinely -- ATAP p.139-140, Bernstein-ellipse clustering).

MEASURED (scratch/30, 4000 random fibers, root-count mismatches vs route B at a fixed
12 iterations):
    circle R=1.0/1.1/1.5/2.0 (the textbook monomial init)   294 / 294 / 305 / 298
    Chebyshev, delta = 0                                    1161   <- worst
    Chebyshev, delta = 0.01 / 0.05 / 0.1 / 0.2                 0 / 0 / 0 / 0
    Chebyshev, delta = 0.5 / 1.0                               9 / 11
    GLL nodes or equispaced, delta = 0.1                       0 / 0
So the requirement is only (a) real parts spread over [-1,1] and (b) delta nonzero and
not too large.  delta in [0.05,0.2] is a broad plateau; 0.1 sits mid-plateau.  The point
DISTRIBUTION barely matters -- Chebyshev is kept because it is what the roots of a
near-minimax polynomial follow, and it costs one cos/sin pair.

ISOLEVEL.  Subtract c from the nodal values before the transform, exactly as route A.

------------------------------------------------------------------------------
"EIGENSOLVER" CHARACTERISATION  (required by the project brief)
------------------------------------------------------------------------------
    Problem type      NONE -- there is no eigenproblem and no matrix.  This is a
                      simultaneous fixed-point iteration on n complex numbers.
                      (For the record: Ehrlich-Aberth IS mathematically a Newton
                      iteration on the comrade matrix's characteristic polynomial,
                      so it solves the same problem route A's eigensolver does.)
    State per fiber   a[n+1] real + z[n] complex = ~3n doubles.
                      N= 5: 18   N= 8: 27   N=11: 36   N=15: 48
                      Compare route A's n^2 matrix (25/64/121/225) and route B's
                      two (N+2)^2 matrices (98/200/338/578).  A CUDA thread holds
                      ~255 32-bit registers ~ 128 doubles before spilling: route B
                      spills at every N, route A around N=11, route C never.
    Nonzero pattern   n/a (no matrix).  The only all-pairs term is S_i, O(n^2) per
                      iteration, n <= 15.
    Work              niter * [ n evaluations x O(n) recurrence + n repulsion sums
                      x O(n) ] = O(niter * n^2).
                      At niter=12:  N=5: 600 flops   N=8: 1536   N=15: 5400.
                      Dense QR by comparison is O(n^3) ~ 10n^3: 1250 / 5120 / 33750.
                      Ratio QR/Aberth: 2.1x (N=5) to 6.3x (N=15).
    Branches          NONE in the hot loop when niter is fixed.  No pivoting, no
                      shift strategy, no deflation-during-iteration, no convergence
                      test.  This is the property that matters for SIMT.
    Iteration count   FIXED.  Measured (scratch/23, 20,000 random fibers, compared
                      against route B):
                          niter    root-count mismatches    max |.-B|
                             6              57               5.7e-08
                             8               2               4.0e-08
                            10               0               2.0e-13
                            12               0               2.0e-13
                            24               0               2.0e-13
                      DEFAULT niter = 12 (10 suffices; 12 is margin).  Convergence
                      is cubic and then STALLS at machine precision, so extra
                      iterations are harmless -- exactly what a fixed-count kernel
                      needs.
    Accuracy          Matches route B to 2.0e-13 or better on every non-tangency
                      case tested.  On the 8 analytic fields x orders 3..15 the only
                      disagreements are 11/104, ALL on 'tangent' -- a double root,
                      which is sqrt(eps)-conditioned for every method (notes/05).

------------------------------------------------------------------------------
WHY NOT JUST WRITE A QR?  (the question this route had to answer)
------------------------------------------------------------------------------
It is not that a QR is unwritable -- it is ~39 lines for an already-Hessenberg matrix
(route A's last-column form is Hessenberg for free, ATAP Ex.18.2), and a hand-written
complex single-shift Givens QR does work.  It is that at EQUAL fixed cost it is less
accurate and less reliable.  Head-to-head on 4000 random fibers vs route B:

    method                 root-count mismatches     max |.-B|
    QR, 4 sweeps/root              169                1.4e+00
    QR, 6 sweeps/root               25                2.3e-07
    QR, 8 sweeps/root                1                2.2e-08
    Aberth, 10 iterations            0                7.1e-14
    Aberth, 12 iterations            0                7.1e-14

QR needs ~8 sweeps per root to approach reliability and still does not reach zero,
while costing 2.5-7.5x more flops.  It also carries a Wilkinson shift formula, a
deflation criterion (a data-dependent branch), and per-fiber-varying work.  Aberth
reaches machine precision at a fixed 10-12 iterations with none of that.

So the recommendation is not "a QR is too hard to write" -- it is "a QR is the wrong
tool for n <= 15 with a fixed budget".  demo_A keeps the matrix route for reference.

------------------------------------------------------------------------------
LIMITATIONS
------------------------------------------------------------------------------
  * No theorem guarantees all roots are found at a FIXED iteration count.  Empirically
    0 failures in 20,000 fibers at niter >= 10, but a production kernel should carry a
    cheap safety net -- see roots_route_c(..., return_info=True), which reports the
    final correction magnitude; a large value means that fiber did not settle.
  * Double roots are sqrt(eps)-accurate, as for every method (notes/05).  On an
    isosurface these are tangencies: detect them (a returned pair separated by ~1e-8),
    do not trust either root to full precision.
  * Complex arithmetic doubles the register count relative to a real-only scheme.
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "demo_A"))
from gllbase import nodal_to_modal                       # noqa: E402
from route_a import deflate_degree, A_FLOOR_DEFAULT      # noqa: E402

__all__ = ["legendre_val_der", "aberth_init", "aberth_iterate", "roots_route_c",
           "NITER_DEFAULT", "state_doubles"]

# Fixed iteration count.  Justified by a 20,000-fiber sweep -- see the module
# docstring.  10 is enough; 12 is margin.  Cost is exactly linear in this number.
NITER_DEFAULT = 12


def legendre_val_der(a, x):
    """Evaluate p = sum_k a_k P_k and p' at (complex) x by forward recurrence.

    O(n) work, O(1) storage, no matrix.  This is the GPU-shaped kernel primitive:
    the same three-term recurrence that builds route A's comrade matrix, run as an
    evaluation instead.  `x` may be an array (vectorized over roots).
    """
    a = np.asarray(a, dtype=float)
    n = len(a) - 1
    x = np.asarray(x)
    if n == 0:
        return a[0] * np.ones_like(x), np.zeros_like(x)

    P0 = np.ones_like(x)
    P1 = x
    dP0 = np.zeros_like(x)
    dP1 = np.ones_like(x)
    v = a[0] * P0 + a[1] * P1
    d = a[1] * dP1
    for k in range(1, n):
        P2 = ((2 * k + 1) * x * P1 - k * P0) / (k + 1)
        dP2 = ((2 * k + 1) * (P1 + x * dP1) - k * dP0) / (k + 1)
        v = v + a[k + 1] * P2
        d = d + a[k + 1] * dP2
        P0, P1 = P1, P2
        dP0, dP1 = dP1, dP2
    return v, d


def aberth_init(n, offset=0.1):
    """Initial root estimates: Chebyshev points nudged off the real axis.

    The textbook circle-of-radius-R initialization DIVERGES for these polynomials
    (see the module docstring).  The roots of interest are in [-1,1]; `offset`
    lifts the estimates off the real line so the repulsion term is well defined and
    genuinely complex roots remain reachable.
    """
    t = np.pi * (np.arange(n) + 0.5) / n
    return (np.cos(t) + 1j * offset * np.sin(t)).astype(complex)


def aberth_iterate(a, z, niter=NITER_DEFAULT, track=False):
    """Run `niter` Ehrlich-Aberth updates on the estimates z.  No convergence test.

    Returns z (and, if track=True, the per-iteration max |correction| history).
    """
    a = np.asarray(a, dtype=float)
    n = len(z)
    z = np.asarray(z, dtype=complex).copy()
    hist = []
    for _ in range(niter):
        v, d = legendre_val_der(a, z)
        w = np.where(np.abs(d) > 0.0, v / np.where(np.abs(d) > 0.0, d, 1.0), 0.0)
        if n > 1:
            diff = z[:, None] - z[None, :]
            np.fill_diagonal(diff, np.inf)               # excludes j == i
            S = np.sum(1.0 / diff, axis=1)
        else:
            S = np.zeros(1, dtype=complex)
        den = 1.0 - w * S
        corr = np.where(np.abs(den) > 1e-300, w / np.where(np.abs(den) > 1e-300, den, 1.0), w)
        corr = np.where(np.isfinite(corr), corr, 0.0)
        z = z - corr
        if track:
            hist.append(float(np.max(np.abs(corr))) if n else 0.0)
    return (z, hist) if track else z


def state_doubles(N):
    """Per-fiber state in doubles: a[n+1] real + z[n] complex.  For the cost table."""
    return (N + 1) + 2 * N


def roots_route_c(N, f, niter=NITER_DEFAULT, a_floor=A_FLOOR_DEFAULT,
                  imag_tol=1e-7, margin=1e-10, return_info=False, track=False):
    """All real roots in [-1,1] of the GLL interpolant of f (isolevel folded in).

    Matrix-free Ehrlich-Aberth.  Shares route A's transform and degree deflation.

    Parameters
        N         polynomial order (len(f) == N+1)
        f         nodal values of (field - isolevel)
        niter     FIXED iteration count (default 12; >=10 measured safe)
        a_floor   degree-deflation threshold, relative to max|a| (default 1e-12)
        imag_tol  |Im z| < imag_tol counts as a real root
        margin    keep roots on [-1-margin, 1+margin]
        return_info  also return dict(degree=, niter=, final_corr=, spectrum=)
                     -- `final_corr` is the safety net: a large value means this
                     fiber did not settle within `niter`
        track     also return the per-iteration correction history

    Returns sorted real roots, plus whatever the flags request.
    """
    f = np.asarray(f, dtype=float)
    a_full = nodal_to_modal(N) @ f
    a = deflate_degree(a_full, a_floor)
    n = len(a) - 1

    if n < 1:
        empty = np.array([])
        info = dict(degree=0, niter=0, final_corr=0.0, spectrum=np.array([]))
        if return_info and track:
            return empty, info, []
        if return_info:
            return empty, info
        if track:
            return empty, []
        return empty

    z0 = aberth_init(n)
    out = aberth_iterate(a, z0, niter, track=True)
    z, hist = out

    r = np.sort(z[np.abs(z.imag) < imag_tol].real)
    r = r[(r >= -1.0 - margin) & (r <= 1.0 + margin)]

    res = (r,)
    if return_info:
        res += (dict(degree=n, niter=niter,
                     final_corr=hist[-1] if hist else 0.0,
                     spectrum=z),)
    if track:
        res += (hist,)
    return res[0] if len(res) == 1 else res
