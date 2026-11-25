"""
route_a.py -- ROUTE A: modal Legendre comrade matrix, WITH DEGREE DEFLATION.

    GLL nodal values  ->  discrete Legendre coefficients  ->  DEFLATE  ->  comrade  ->  eig

DEFAULT IS "A + DEFLATION" (a_floor = 1e-12).  Plain route A (a_floor = 0) divides by
the top Legendre coefficient a_N and is silently wrong when a_N sits at the noise
floor; deflation drops trailing negligible coefficients first, which removes that
failure entirely (0 disagreements with route B over 20,000 randomized draws) while
keeping a STANDARD eigenproblem.  See "DEGREE DEFLATION" below and notes/07, notes/08.

This is Trefethen ATAP Thm 18.1 (colleague, Chebyshev) transplanted to the Legendre
basis, i.e. ATAP Exercise 18.4 ("comrade matrix", Barnett 1975).

------------------------------------------------------------------------------
MATH
------------------------------------------------------------------------------
Step 1.  Nodal -> modal.
    The GLL interpolant of order N is written in the Legendre basis
        p(x) = sum_{k=0}^{N} a_k P_k(x)
    with the discrete projection
        a_k = (1/gamma_k) sum_j w_j f_j P_k(x_j),     gamma_k = <P_k,P_k>_GLL.
    GLL quadrature is exact to degree 2N-1, so gamma_k = 2/(2k+1) for k < N, but
        gamma_N = 2/N        (NOT 2/(2N+1))
    because P_N^2 has degree 2N.  This is the classic GLL top-mode inflation.  It
    matters here more than usual: a_N is the number step 3 divides by.

Step 2.  Multiplication by x.
    Legendre obeys the Sturm-Liouville three-term recurrence
        (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x),
    equivalently
        x P_k(x) = (k/(2k+1)) P_{k-1}(x) + ((k+1)/(2k+1)) P_{k+1}(x),
        x P_0(x) = P_1(x).
    So on span{P_0..P_{N-1}} the operator "multiply by x" is TRIDIAGONAL:
        C[k, k-1] = k/(2k+1),      C[k, k+1] = (k+1)/(2k+1),     C[0,1] = 1.

Step 3.  Close the system (the spill-over row).
    Row k = N-1 needs P_N, which is outside the span.  At a root of p,
        P_N = -(1/a_N) sum_{k<N} a_k P_k,
    so substituting into  x P_{N-1} = ((N-1)/(2N-1)) P_{N-2} + (N/(2N-1)) P_N  gives
        C[N-1, :] -= (N/(2N-1)) * a[0:N] / a_N.
    The eigenvalues of the resulting NxN matrix C are exactly the roots of p.

The isolevel c enters once, by subtracting it from the nodal values before step 1
(equivalently a_0 -= c, since P_0 = 1).

------------------------------------------------------------------------------
EIGENPROBLEM CHARACTERISATION  (required by the project brief)
------------------------------------------------------------------------------
    Problem type      standard, non-symmetric:  C v = lambda v
    Size              N x N              (degree N interpolant -> N roots)
    Structure         tridiagonal EXCEPT the last row, which is full.
                      -> the matrix is already upper Hessenberg apart from the
                         subdiagonal, and is in fact a rank-1 perturbation of a
                         symmetrizable tridiagonal (Jacobi-like) matrix:
                             C = J_hat - (N/(2N-1)) e_{N-1} (a_0..a_{N-1})^T / a_N
    Nonzero pattern   nnz = 3N - 3 exactly (measured, generic coefficients).
                      Derivation: the tridiagonal part has a ZERO DIAGONAL (the
                      Legendre recurrence has no P_k term in x*P_k), so it carries
                      2(N-1) entries, not 3N-2; the full last row adds N, of which
                      one was already counted.  2(N-1) + (N-1) = 3N-3.
                          N= 5: 12/25  = 48.0%      N=10: 27/100 = 27.0%
                          N= 8: 21/64  = 32.8%      N=15: 42/225 = 18.7%
                      So the matrix is sparse in principle, but at these sizes a
                      sparse solver would lose to dense LAPACK on overhead alone.
    Symmetry          NONE.  The tridiagonal part is symmetrizable (it is the
                      Jacobi matrix of the Legendre measure up to a diagonal
                      similarity) but the last-row update destroys that, so no
                      symmetric solver applies.
    Eigenvalues       real ones are the roots in the real line; complex ones come
                      in conjugate pairs and are discarded (they are genuine roots
                      of the interpolant, just not on the real axis -- ATAP p.139-140,
                      the Walsh/Bernstein-ellipse clustering).
    LAPACK path       scipy.linalg.eig -> dgeev = dgebal (balance) + dgehrd
                      (reduce to Hessenberg) + dhseqr (Francis QR).
                      Balancing is ESSENTIAL here (ATAP p.139, Parlett-Reinsch 1969):
                      the last row can be huge when a_N is small, and balancing is
                      what keeps the computed roots accurate.  Never pass
                      'nobalance'.
    Cost              dgehrd  (10/3) N^3 flops   -- wasted: C is nearly Hessenberg
                                                    already, but LAPACK has no
                                                    "skip it" entry point exposed
                                                    by scipy (no dhseqr binding).
                      dhseqr  ~10 N^3 flops (eigenvalues only, ~2-3 Francis sweeps
                                             per eigenvalue)
                      TOTAL   O(N^3), ~13 N^3 flops.
                      Step 1 (nodal->modal) adds 2(N+1)^2 -- negligible.
                      For N = 5..15 this is 1.6k..44k flops, i.e. entirely
                      latency/overhead-bound in practice.  The measured time is
                      NOT the flop count; see bench_A.py.
    Exploitable?      Yes in principle: a comrade matrix is diagonal-plus-rank-one
                      in disguise, and O(N^2) structured QR variants exist
                      (Aurentz-Mach-Vandebril-Watkins et al.).  Not pursued here --
                      at N<=15 the dense solve is already overhead-dominated.

------------------------------------------------------------------------------
THE FAILURE MODE, AND THE FIX (DEGREE DEFLATION)
------------------------------------------------------------------------------
Step 3 divides by a_N.  The dangerous case is NOT simply "a_N is small" -- a
well-resolved element has a small a_N and route A handles it fine (0 discrepancies
in 4000 randomized trials with a realistic spectral decay, notes/07).  The failure
is a_N small WHILE THE LOW MODES STAY O(1): a polynomial whose true degree is < N
being forced through a degree-N comrade matrix.  That happens for real -- e.g. a
degree-5 field sampled at N=6 -- and it fails SILENTLY: at N=6 on 'poly5' route A
returns five plausible real roots in [-1,1] with a max error of 4.2e-01, while
cond(C) = 2.4e17 and the last row carries entries ~1e15.

THE FIX.  Before building C, drop trailing coefficients that are negligible:

    keep a[0:d+1]  where  d = max{k : |a_k| > a_floor * max_j |a_j|}

then build the d x d comrade matrix.  The polynomial's actual degree is used, so
the division is by a genuine leading coefficient rather than by roundoff.

Threshold.  a_floor = 1e-12 by default.  Measured (scratch/16), 20,000 randomized
draws with 25% forced degree-deficient, "mismatch" = disagreeing with route B on the
NUMBER of roots:

    a_floor    mismatches   max|A-B|    deflation fired
    0 (off)        24        2.8e-01        0.0%      <- plain route A, broken
    1e-16          17        2.8e-01        8.3%
    1e-14           0        3.1e-04       54.1%
    1e-12           0        2.8e-13       61.1%      <- DEFAULT
    1e-10           0        2.6e-11       68.0%
    1e-08           9        7.5e-09       75.0%
    1e-06         792        1.1e-06       82.1%      <- discards real leading terms
    1e-02        6810        1.8e-02       96.1%

1e-12 and 1e-10 are both in the safe plateau; 1e-12 is chosen as the tighter of the
two.  Above 1e-8 the threshold starts deleting genuine leading coefficients.  On the
8 analytic test fields x orders 3..15, deflation at 1e-12 disagrees with route B on
ZERO non-tangency cases (plain route A disagrees on 2).

Cost of deflation: one O(N) reduction plus a smaller eigenproblem.  19.3 -> 23.7 us
at N=5; still cheaper than route B (28.0 us).

WHAT DEFLATION DOES NOT FIX.  A double root is sqrt(eps)-conditioned for every method
(notes/05) -- deflation is irrelevant there.  And on GPU the deflated size is
data-dependent, so batches must be bucketed by degree; measured degree distributions
say this is rare (at N=8 with a realistic decay the degree is full ~100% of the
time).  See notes/08.
"""
import numpy as np
import scipy.linalg as sla
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from gllbase import gll, nodal_to_modal    # noqa: E402

__all__ = ["modal_coeffs", "comrade_matrix", "comrade_nnz", "roots_route_a",
           "deflate_degree", "A_FLOOR_DEFAULT"]

# Default degree-deflation threshold.  Justified by a 20,000-draw sweep -- see the
# module docstring and scratch/16.  Set a_floor=0.0 to recover plain route A.
A_FLOOR_DEFAULT = 1e-12


def modal_coeffs(N, f):
    """Discrete Legendre coefficients a (length N+1) from GLL nodal values f."""
    return nodal_to_modal(N) @ np.asarray(f, dtype=float)


def comrade_matrix(a):
    """Legendre comrade matrix C (n x n, n = deg) for p = sum_k a_k P_k.

    Raises ZeroDivisionError-like ValueError if the leading coefficient vanishes:
    the construction is undefined there.  That is route A's structural limitation,
    not an implementation bug -- see the module docstring.
    """
    a = np.asarray(a, dtype=float)
    n = len(a) - 1
    if n < 1:
        raise ValueError("need degree >= 1")
    if a[n] == 0.0:
        raise ValueError("leading coefficient a_N is exactly zero: "
                         "the comrade matrix is undefined (use route B)")
    C = np.zeros((n, n))
    for k in range(n):                      # x P_k = k/(2k+1) P_{k-1} + (k+1)/(2k+1) P_{k+1}
        if k > 0:
            C[k, k - 1] = k / (2 * k + 1)
        if k + 1 <= n - 1:
            C[k, k + 1] = (k + 1) / (2 * k + 1)
    C[0, 1] = 1.0 if n > 1 else 0.0         # x P_0 = P_1
    C[n - 1, :] -= (n / (2 * n - 1)) * a[:n] / a[n]     # spill-over row
    return C


def deflate_degree(a, a_floor=A_FLOOR_DEFAULT):
    """Trim trailing negligible Legendre coefficients; return the truncated vector.

    d = max{k : |a_k| > a_floor * max_j |a_j|}, and a[0:d+1] is returned.  With
    a_floor = 0 the input is returned unchanged (plain route A).  A constant
    polynomial deflates to length 1, which the caller treats as "no roots".
    Cost: O(N), one pass.
    """
    a = np.asarray(a, dtype=float)
    if a_floor <= 0.0:
        return a
    scale = np.max(np.abs(a))
    if scale == 0.0:
        return a[:1]
    nz = np.nonzero(np.abs(a) > a_floor * scale)[0]
    return a[:(nz[-1] + 1)] if len(nz) else a[:1]


def comrade_nnz(N):
    """(nnz, dense_entries) of the NxN comrade matrix -- for the cost table."""
    C = comrade_matrix(np.r_[np.ones(N), 1.0])
    return int(np.count_nonzero(C)), N * N


def roots_route_a(N, f, imag_tol=1e-8, margin=1e-12, return_all=False,
                  a_floor=A_FLOOR_DEFAULT, return_info=False):
    """All real roots in [-1,1] of the GLL interpolant of f (isolevel folded into f).

    Route A + degree deflation.  Deflation is ON by default; it is what makes this
    route safe to use (see the module docstring).

    Parameters
        N         polynomial order (len(f) == N+1)
        f         nodal values of (field - isolevel) at the GLL nodes
        imag_tol  an eigenvalue counts as real if |Im| <= imag_tol * max(1,|Re|)
        margin    roots are kept on [-1-margin, 1+margin] so a root sitting exactly
                  on a face is not lost to rounding
        a_floor   degree-deflation threshold, relative to max|a| (default 1e-12).
                  Pass 0.0 for PLAIN route A -- fragile, for comparison only.
        return_all   also return the raw complex spectrum
        return_info  also return dict(degree=, deflated=, a_N_ratio=)

    Returns sorted real roots, plus whatever the return_* flags ask for.
    """
    f = np.asarray(f, dtype=float)
    a_full = modal_coeffs(N, f)
    a = deflate_degree(a_full, a_floor)

    scale = np.max(np.abs(a_full))
    info = dict(degree=len(a) - 1,
                deflated=(len(a) - 1) < N,
                a_N_ratio=(abs(a_full[N]) / scale if scale > 0 else 0.0))

    def _ret(r, lam=None):
        out = (r,)
        if return_all:
            out += (lam if lam is not None else np.array([]),)
        if return_info:
            out += (info,)
        return out[0] if len(out) == 1 else out

    if len(a) < 2 or a[-1] == 0.0:
        # constant (or fully deflated): no isolated roots
        return _ret(np.array([]), np.array([]))

    C = comrade_matrix(a)
    lam = sla.eig(C, right=False)           # dgeev: balance + Hessenberg + Francis QR
    keep = np.abs(lam.imag) <= imag_tol * np.maximum(1.0, np.abs(lam.real))
    r = np.sort(lam[keep].real)
    r = r[(r >= -1.0 - margin) & (r <= 1.0 + margin)]
    return _ret(r, lam)
