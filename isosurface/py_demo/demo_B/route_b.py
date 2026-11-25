"""
route_b.py -- ROUTE B: nodal barycentric (Lagrange-basis) arrowhead pencil.

    GLL nodal values  ->  arrowhead matrix pencil (A,B)  ->  generalized eig

No basis transform.  This is the "eigenvalue problems constructed from function
values themselves" idea Trefethen flags as under-exploited at ATAP p.139
(Amiraslani, Corless, Lancaster 2004/2006; Corless-Watt; Lawrence-Corless).

------------------------------------------------------------------------------
MATH
------------------------------------------------------------------------------
Expand in the LAGRANGE cardinal basis on the GLL nodes, l_j(x_i) = delta_ij:

        p(x) = sum_j f_j l_j(x)

so the expansion coefficients ARE the nodal values -- no projection, no transform,
no gamma_N subtlety, nothing to get wrong.  That is the whole appeal.

The obstruction: {l_j} has NO three-term recurrence, because x*l_j(x) has degree
N+1 and is not a combination of the l_i (degree N).  So the "multiply by x"
operator is not a matrix on span{l_j} -- Step 2 of the colleague recipe fails.

The fix: allow a GENERALIZED eigenproblem A v = lambda B v with B singular.

Barycentric form.  With
        b_j = 1 / prod_{i != j} (x_j - x_i)          (barycentric weights)
        l(x) = prod_i (x - x_i)                      (node polynomial)
the interpolant is
        p(x) = l(x) * sum_j b_j f_j / (x - x_j) .
Away from the nodes, roots of p are zeros of the SECULAR function
        s(x) = sum_j b_j f_j / (x - x_j) .
Clearing denominators of a secular equation is the classical ARROWHEAD
eigenproblem:

            [ 0    f_0  f_1  ...  f_N ]           [ 0                 ]
            [ b_0  x_0                ]           [    1              ]
      A  =  [ b_1       x_1           ] ,   B  =  [       1           ]
            [  :             ...      ]           [           ...     ]
            [ b_N                 x_N ]           [               1   ]

both (N+2) x (N+2).  B is the identity with a single zero in the (0,0) entry.

Determinant identity (expand along the first row/column of the arrowhead):
        det(A - lambda B) = - sum_j b_j f_j prod_{i != j} (x_i - lambda)
                          = -/+ p(lambda)
i.e. the pencil's characteristic polynomial IS the interpolant.  This is an
IDENTITY, the same status as ATAP Thm 18.1 -- verified to 12 digits in
notes/02 and re-checked by run_demo.py --study identity.

Therefore: the FINITE eigenvalues of (A,B) are exactly the roots of p.  The
singular (0,0) entry of B sends one eigenvalue to infinity, so of the N+2
eigenvalues, N+1 - (number of degree deficiencies) are finite.

Root exactly AT a node.  The secular function seems to have a pole at x_j, but
when f_j = 0 the residue b_j f_j vanishes with it; and the pencil is assembled
from f_j and b_j directly, never from a quotient.  Verified to 1e-15 (notes/02).

Isolevel c:  put f_j - c in the TOP ROW.  Do not shift the diagonal (the diagonal
holds the nodes, not values).

------------------------------------------------------------------------------
EIGENPROBLEM CHARACTERISATION  (required by the project brief)
------------------------------------------------------------------------------
    Problem type      GENERALIZED, non-symmetric, SINGULAR B:  A v = lambda B v
    Size              (N+2) x (N+2)      -- note: larger than route A's N x N
    Structure         A: ARROWHEAD (bordered diagonal) -- diagonal x_j, first row
                         f_j, first column b_j, zero elsewhere.
                      B: diagonal, = diag(0, 1, 1, ..., 1).  Exactly singular.
    Nonzero pattern   nnz(A) = 3N + 3  (measured: N+1 diagonal + N+1 row + N+1 col)
                          N= 5: 18/49   = 36.7%     N=10: 33/144 = 22.9%
                          N= 8: 27/100  = 27.0%     N=15: 48/289 = 16.6%
                      nnz(B) = N+1.
    Symmetry          A is NOT symmetric (row holds f, column holds b).  It becomes
                      symmetric iff f_j = b_j.  A diagonal similarity D A D^-1 with
                      D = diag(1, sqrt(|f_j/b_j|)) symmetrizes it when the signs of
                      f_j b_j are all positive -- NOT generally true, so no symmetric
                      solver applies.  (Do not "fix" this by scaling: the sign
                      pattern of f is the data.)
    Eigenvalues       N+2 total: >= 1 at infinity (always), plus one extra infinite
                      eigenvalue per degree deficiency.  Finite ones = roots.
                      DEGREE DEFICIENCY IS REPORTED, NOT DIVIDED BY -- the central
                      advantage over route A.
    LAPACK path       scipy.linalg.eig(A, B) -> dggev = dggbal + dgghrd (reduce to
                      generalized Hessenberg-triangular) + dhgeqz (QZ iteration).
                      Use homogeneous_eigvals=True to get (alpha, beta) and filter
                      infinities by the RELATIVE test |beta| > tol*|alpha|.  Never
                      compute alpha/beta first and threshold on magnitude.
    Cost              dgghrd ~ (20/3) n^3,  dhgeqz ~ 30 n^3 (eigenvalues only)
                      with n = N+2.  TOTAL O(N^3) but with a constant roughly
                      3x route A's, on a matrix 2-3 ranks larger:
                          flops(B)/flops(A) ~ (37/13) * ((N+2)/N)^3
                          N= 5: ~7.8x       N=10: ~4.9x      N=15: ~4.2x
                      MEASURED wall time is ~4x flat across N=5..15 (bench_B.py),
                      because at this size both are overhead-bound.
    Exploitable?      Yes, more so than route A: arrowhead eigenproblems admit
                      O(n^2) (even O(n log n) with fast multipole) solvers via
                      secular-equation deflation -- the same machinery as the
                      divide-and-conquer symmetric tridiagonal eigensolver.
                      Not pursued here; noted as the obvious optimisation if the
                      isosurface sweep ever becomes eigensolve-bound.

------------------------------------------------------------------------------
WHY THIS ROUTE ON GLL SPECIFICALLY
------------------------------------------------------------------------------
The barycentric form is only as good as the conditioning of the weights b_j.  The
spread max|b|/min|b| grows like O(N) on GLL nodes but like O(4^N/N) on equispaced
nodes (measured, notes/02):  N=64 gives 10.1 on GLL vs 1.8e18 equispaced.  So this
construction is well conditioned precisely on the node set SEM already uses.
"""
import numpy as np
import scipy.linalg as sla
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from gllbase import gll, gll_bary_weights      # noqa: E402

__all__ = ["arrowhead_pencil", "pencil_nnz", "roots_route_b"]


def arrowhead_pencil(N, f):
    """Assemble the (N+2)x(N+2) arrowhead pencil (A, B) for GLL nodal values f.

    Cost: O(N) to fill (the matrices are mostly zero); allocation dominates.
    """
    x, _ = gll(N)
    b = gll_bary_weights(N)
    f = np.asarray(f, dtype=float)
    n = N + 2
    A = np.zeros((n, n))
    A[0, 1:] = f
    A[1:, 0] = b
    A[1:, 1:] = np.diag(x)
    Bm = np.eye(n)
    Bm[0, 0] = 0.0
    return A, Bm


def pencil_nnz(N):
    """(nnz(A), nnz(B), dense_entries) -- for the cost table."""
    A, Bm = arrowhead_pencil(N, np.ones(N + 1))
    return int(np.count_nonzero(A)), int(np.count_nonzero(Bm)), (N + 2) ** 2


def roots_route_b(N, f, imag_tol=1e-8, margin=1e-12, inf_tol=1e-8,
                  return_info=False):
    """All real roots in [-1,1] of the GLL interpolant of f (isolevel folded in).

    Parameters
        N         polynomial order (len(f) == N+1)
        f         nodal values of (field - isolevel)
        imag_tol  |Im| <= imag_tol * max(1,|Re|) counts as real
        margin    keep roots on [-1-margin, 1+margin]
        inf_tol   an eigenvalue is INFINITE if |beta| <= inf_tol * |alpha|
                  (relative test on the homogeneous pair -- see module docstring)
        return_info  also return dict(n_finite=..., n_infinite=..., spectrum=...)

    Never divides by a leading coefficient: degree deficiency simply raises
    n_infinite above its floor of 1.
    """
    A, Bm = arrowhead_pencil(N, f)
    alpha, beta = sla.eig(A, Bm, right=False, homogeneous_eigvals=True)

    finite = np.abs(beta) > inf_tol * np.abs(alpha)
    lam = alpha[finite] / beta[finite]

    keep = np.abs(lam.imag) <= imag_tol * np.maximum(1.0, np.abs(lam.real))
    r = np.sort(lam[keep].real)
    r = r[(r >= -1.0 - margin) & (r <= 1.0 + margin)]

    if return_info:
        return r, dict(n_finite=int(finite.sum()),
                       n_infinite=int((~finite).sum()),
                       spectrum=lam)
    return r
