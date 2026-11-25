"""
gllbase.py -- shared primitives for the GLL global-rootfinding demos.

Single reference element  E = [-1, 1]  with GLL nodes of order N (N+1 nodes).
Nothing here is route-specific: routes live in demo_A/route_a.py, demo_B/route_b.py.

Contents
    gll(N)                  GLL nodes + quadrature weights
    bary_weights(x)         barycentric weights  b_j = 1/prod_{i!=j}(x_j-x_i)
    gll_bary_weights(N)     the same, via the closed form  (-1)^j sqrt(w_j)
    legendre_vandermonde    V[j,k] = P_k(x_j)
    nodal_to_modal(N)       M with  a = M f   (discrete Legendre coefficients)
    modal_to_nodal(N)       V with  f = V a
    TESTS                   catalogue of test fields with ANALYTIC exact roots
    match_roots             error metric that tolerates missed/spurious roots
    ref_roots_companion     reference: monomial companion  (numpy.roots)
    ref_roots_brent         reference: dense bracketing + Brent
    timeit_us               median wall time per call, microseconds
"""
import numpy as np
from numpy.polynomial import legendre as L
from functools import lru_cache
import time

__all__ = ["gll", "bary_weights", "gll_bary_weights", "legendre_vandermonde",
           "nodal_to_modal", "modal_to_nodal", "TESTS", "TestField",
           "match_roots", "ref_roots_companion", "ref_roots_brent", "timeit_us",
           "interp_error", "root_error_floor"]


# --------------------------------------------------------------------------
# GLL basis
# --------------------------------------------------------------------------
@lru_cache(maxsize=64)
def gll(N):
    """N+1 Gauss-Lobatto-Legendre nodes and weights on [-1,1] (polynomial order N).

    Interior nodes are the roots of P'_N; weights are w_j = 2/(N(N+1) P_N(x_j)^2).
    Returns (x, w) as read-only arrays, cached per N.
    """
    if N < 1:
        raise ValueError("N >= 1 required")
    if N == 1:
        x = np.array([-1.0, 1.0])
    else:
        dP = L.legder([0] * N + [1])
        x = np.concatenate(([-1.0], np.sort(L.legroots(dP)), [1.0]))
    PN = L.legval(x, [0] * N + [1])
    w = 2.0 / (N * (N + 1) * PN ** 2)
    x.flags.writeable = False
    w.flags.writeable = False
    return x, w


def bary_weights(x):
    """Barycentric weights  b_j = 1 / prod_{i != j} (x_j - x_i).   O(n^2)."""
    n = len(x)
    b = np.empty(n)
    for j in range(n):
        b[j] = 1.0 / np.prod(x[j] - np.delete(x, j))
    return b


@lru_cache(maxsize=64)
def gll_bary_weights(N):
    """Barycentric weights of the GLL nodes in closed form:  b_j = (-1)^j sqrt(w_j).

    PROPORTIONAL to the true barycentric weights, not equal: b_true = kappa(N) * b_here
    with kappa = 2.77 (N=4), 33.5 (N=8), -3241 (N=15).  The ratio is constant across j
    to 1e-14, and route B's pencil is homogeneous of degree 1 in b, so kappa cancels
    and every root is unaffected.  It does show up in det(A - lam B) -- see
    demo_B/run_demo.py --study identity, and notes/02.  Use bary_weights(x) if an
    unscaled weight is ever required.
    """
    _, w = gll(N)
    b = np.array([(-1) ** j * np.sqrt(w[j]) for j in range(N + 1)])
    b.flags.writeable = False
    return b


@lru_cache(maxsize=64)
def legendre_vandermonde(N):
    """V[j,k] = P_k(x_j), k = 0..N, at the GLL nodes.  Modal -> nodal:  f = V a."""
    x, _ = gll(N)
    V = np.empty((N + 1, N + 1))
    for k in range(N + 1):
        V[:, k] = L.legval(x, [0] * k + [1])
    V.flags.writeable = False
    return V


@lru_cache(maxsize=64)
def _discrete_norms(N):
    """gamma_k = <P_k,P_k>_GLL.  NOTE the top mode: 2/N, not 2/(2N+1).

    GLL quadrature is exact to degree 2N-1; P_N^2 has degree 2N, so the top mode
    alone is integrated inexactly.  Getting this wrong corrupts a_N -- i.e. exactly
    the coefficient route A divides by.  See notes/01.
    """
    g = np.array([2.0 / (2 * k + 1) for k in range(N + 1)])
    g[N] = 2.0 / N
    g.flags.writeable = False
    return g


@lru_cache(maxsize=64)
def nodal_to_modal(N):
    """M with  a = M f : discrete Legendre coefficients from GLL nodal values.

    a_k = (1/gamma_k) sum_j w_j f_j P_k(x_j)   ->   M = diag(1/gamma) V^T diag(w).
    Dense (N+1)x(N+1); applying it costs 2(N+1)^2 flops.  Exact (not a least-squares
    fit): the GLL interpolant and this discrete projection are the same polynomial.
    """
    _, w = gll(N)
    V = legendre_vandermonde(N)
    M = (V * w[:, None]).T / _discrete_norms(N)[:, None]
    M.flags.writeable = False
    return M


def modal_to_nodal(N):
    return legendre_vandermonde(N)


# --------------------------------------------------------------------------
# Test fields with ANALYTIC exact roots  (no numerical oracle needed)
# --------------------------------------------------------------------------
class TestField:
    """A scalar field f on [-1,1] and an isolevel c, with the exact roots of f-c."""

    def __init__(self, key, desc, f, exact_roots, level=0.0, poly_degree=None):
        self.key = key
        self.desc = desc
        self._f = f
        self._exact = exact_roots
        self.level = level
        self.poly_degree = poly_degree          # None if not a polynomial

    def values(self, x):
        """Nodal values of  f - c  (the isolevel is folded in once, here)."""
        return self._f(np.asarray(x)) - self.level

    def raw(self, x):
        return self._f(np.asarray(x))

    def exact(self):
        r = np.sort(np.asarray(self._exact, dtype=float))
        return r[(r >= -1.0) & (r <= 1.0)]

    def __repr__(self):
        return f"<TestField {self.key}: {self.desc}, {len(self.exact())} roots>"


def _shifted_sine_roots(k, x0, c):
    """All x in [-1,1] with sin(k*pi*(x-x0)) = c."""
    out = []
    a = np.arcsin(c)
    for m in range(-2 * k - 4, 2 * k + 5):
        for t in (a + 2 * np.pi * m, np.pi - a + 2 * np.pi * m):
            x = t / (k * np.pi) + x0
            if -1.0 <= x <= 1.0:
                out.append(x)
    return np.unique(np.round(out, 14))


def _build_tests():
    """Test fields with ANALYTIC roots.

    DESIGN NOTE (learned the hard way -- see notes/03).  Every field here is
    deliberately NOT even and NOT odd about x=0.  A field with definite parity has
    a_N = 0 identically for every N of the wrong parity, which (a) makes route A
    hit its degenerate branch for a reason that has nothing to do with resolution,
    and (b) makes a convergence table zig-zag by orders of magnitude between
    consecutive N.  The centres below (0.15, 0.2, ...) exist purely to break parity.
    """
    T = {}

    # -- exactly representable polynomial: the machine-precision check ---------
    r5 = np.array([-0.9, -0.42, 0.07, 0.37, 0.88])
    T["poly5"] = TestField(
        "poly5", "degree-5 polynomial, 5 prescribed roots (exact for N>=5)",
        lambda x, r=r5: np.prod([x - ri for ri in r], axis=0), r5, 0.0, 5)

    # -- smooth, entire, spectrally easy (off-centre: breaks parity) ----------
    a_g, x0g, c_g = 4.0, 0.15, 0.5
    dg = np.sqrt(-np.log(c_g) / a_g)
    T["gauss"] = TestField(
        "gauss", "exp(-4(x-0.15)^2) = 0.5   (2 roots, entire, fast spectral)",
        lambda x: np.exp(-a_g * (x - x0g) ** 2),
        np.array([x0g - dg, x0g + dg]), c_g)

    # -- steep internal front (boundary-layer-like) ---------------------------
    a_t, x0t, c_t = 5.0, 0.2, 0.3
    T["tanh"] = TestField(
        "tanh", "tanh(5(x-0.2)) = 0.3   (1 root, steep front)",
        lambda x: np.tanh(a_t * (x - x0t)),
        np.array([x0t + np.arctanh(c_t) / a_t]), c_t)

    # -- Runge: analytic but poles at +-i/5 -> slow, algebraic-looking rate ----
    #    shifted off-centre so a_N does not vanish by parity
    c_r, x0r = 0.5, 0.1
    dr = np.sqrt((1.0 / c_r - 1.0) / 25.0)
    T["runge"] = TestField(
        "runge", "1/(1+25(x-0.1)^2) = 0.5   (2 roots, poles near [-1,1])",
        lambda x: 1.0 / (1.0 + 25.0 * (x - x0r) ** 2),
        np.array([x0r - dr, x0r + dr]), c_r)

    # -- oscillatory: several roots inside ONE element ------------------------
    T["sine2"] = TestField(
        "sine2", "sin(2 pi (x-0.05)) = 0.4   (4 roots, resolved by N ~ 11)",
        lambda x: np.sin(2 * np.pi * (x - 0.05)),
        _shifted_sine_roots(2, 0.05, 0.4), 0.4)

    T["sine3"] = TestField(
        "sine3", "sin(3 pi (x-0.05)) = 0.5   (6 roots, resolved by N ~ 15)",
        lambda x: np.sin(3 * np.pi * (x - 0.05)),
        _shifted_sine_roots(3, 0.05, 0.5), 0.5)

    # -- deliberately UNDER-resolved at every N tested ------------------------
    T["sine5"] = TestField(
        "sine5", "sin(5 pi (x-0.03)) = 0.3   (10 roots, under-resolved for N<=15)",
        lambda x: np.sin(5 * np.pi * (x - 0.03)),
        _shifted_sine_roots(5, 0.03, 0.3), 0.3)

    # -- degenerate: double root (tangency to the isolevel) -------------------
    T["tangent"] = TestField(
        "tangent", "(x-0.25)^2 (x+0.6) = 0   (DOUBLE root at 0.25 + simple at -0.6)",
        lambda x: (x - 0.25) ** 2 * (x + 0.6),
        np.array([-0.6, 0.25, 0.25]), 0.0, 3)

    return T


TESTS = _build_tests()


# --------------------------------------------------------------------------
# Error metric
# --------------------------------------------------------------------------
def match_roots(found, exact, dedup_tol=1e-8):
    """Compare a computed root set against the exact one.

    Returns dict:
        max_err   largest |found - exact| over the matched pairs (inf if none matched)
        n_found   number of (deduplicated) computed roots
        n_exact   number of exact roots (multiplicity collapsed the same way)
        n_missed  exact roots with no computed root within dedup_tol... reported as
                  a count of unmatched exact roots after greedy nearest matching
        n_spurious computed roots not matched to any exact root
    Multiple roots are collapsed on BOTH sides before matching, so a double root
    counts once -- an eigensolver returns a cluster, not a labelled multiplicity.
    """
    f = np.sort(np.asarray(found, dtype=float).ravel())
    e = np.sort(np.asarray(exact, dtype=float).ravel())
    f = _collapse(f, dedup_tol)
    e = _collapse(e, dedup_tol)
    if len(f) == 0 or len(e) == 0:
        return dict(max_err=np.inf if len(e) else 0.0, n_found=len(f),
                    n_exact=len(e), n_missed=len(e), n_spurious=len(f))
    # greedy nearest matching (sets are tiny and sorted)
    used = np.zeros(len(f), dtype=bool)
    errs, missed = [], 0
    for ei in e:
        d = np.abs(f - ei)
        d[used] = np.inf
        k = int(np.argmin(d))
        if np.isinf(d[k]):
            missed += 1
        else:
            used[k] = True
            errs.append(d[k])
    return dict(max_err=max(errs) if errs else np.inf,
                n_found=len(f), n_exact=len(e), n_missed=missed,
                n_spurious=int(np.sum(~used)))


def _collapse(r, tol):
    if len(r) == 0:
        return r
    out = [r[0]]
    for v in r[1:]:
        if abs(v - out[-1]) > tol:
            out.append(v)
    return np.array(out)


def interp_error(N, field, nsample=6001):
    """max |GLL interpolant - true field| on [-1,1]: the error floor for any root.

    A rootfinder cannot beat the interpolant it is given.  This is the quantity to
    compare a root error against -- an error at this level means the ROOTFINDER is
    exact and the INTERPOLATION is the limit.
    """
    from numpy.polynomial import legendre as L
    x, _ = gll(N)
    a = nodal_to_modal(N) @ field.values(x)
    xs = np.linspace(-1.0, 1.0, nsample)
    return float(np.max(np.abs(L.legval(xs, a) - field.values(xs))))


def root_error_floor(N, field, nsample=6001, kind="pointwise"):
    """The error a PERFECT rootfinder would still make -- i.e. the CONDITIONING of
    the rootfinding problem, driven by the interpolation discrepancy.

    WHAT delta-p IS.  The eigensolver does not see the true field f.  It sees the
    GLL interpolant p_N, and it solves p_N(x) = 0 essentially exactly (measured:
    4e-16 against high-precision roots of p_N itself, scratch/06).  So the error
    against the TRUE root x* of f is not the eigensolver's; it is

        delta-p (x)  =  p_N(x) - f(x)          the INTERPOLATION discrepancy

    evaluated where it matters, at the root.  First-order perturbation of a simple
    root (this is exactly the classical root condition number 1/|p'|):

        p_N(x* + dx) = 0,   p_N = f + delta-p,   f(x*) = 0
        =>  0 = f(x*) + delta-p(x*) + (f' + delta-p')(x*) dx + O(dx^2)
        =>  dx  ~=  - delta-p(x*) / p_N'(x*)
        =>  |dx|  <=  |delta-p(x*)| / |p_N'(x*)|

    So: NUMERATOR = how wrong the polynomial is AT the root;
        DENOMINATOR = how steeply it crosses -- a shallow crossing amplifies the
        same discrepancy into a large root displacement.  1/|p'(x*)| IS the
        condition number of the root.

    kind="pointwise"  |delta-p(x*)| / |p_N'(x*)|     -- sharp; matches the measured
                      root error to 1.00x across the test set (scratch/06)
    kind="max"        ||delta-p||_inf / |p_N'(x*)|   -- a rigorous upper bound, but
                      loose by 1-47x because the sup-norm is usually attained far
                      from the root (e.g. at the element ends)

    DOUBLE ROOTS.  If p'(x*) = 0 this formula diverges, correctly signalling that
    the linear theory fails: a double root splits as sqrt(delta-p / (p''/2)), so a
    1e-16 discrepancy becomes a ~1e-8 displacement.  Measured on the 'tangent'
    field: the double root at 0.25 is returned as a pair 0.24999999/0.25000001,
    displaced by ~1e-8 = sqrt(eps).  That is not a bug in either route -- it is the
    conditioning of a double root, and no rootfinder can do better in double
    precision.  Reported as inf here.
    """
    from numpy.polynomial import legendre as L
    x, _ = gll(N)
    a = nodal_to_modal(N) @ field.values(x)
    ex = field.exact()
    if len(ex) == 0:
        return 0.0
    if kind == "max":
        dp = np.full(len(ex), interp_error(N, field, nsample))
    else:
        dp = np.abs(L.legval(ex, a) - field.values(ex))
    der = np.abs(L.legval(ex, L.legder(a)))
    with np.errstate(divide="ignore", invalid="ignore"):
        fl = dp / der
    return float(np.max(fl))


# --------------------------------------------------------------------------
# Reference rootfinders (for contrast, not for production)
# --------------------------------------------------------------------------
def ref_roots_companion(N, f):
    """REFERENCE: fit the monomial form and use numpy.roots (monomial companion).

    This is the method Trefethen's Ch.18 opens by warning against: the monomial
    basis makes rootfinding ill-conditioned (Wilkinson).  Included to show the
    contrast, not as a candidate.
    Eigenproblem: dense companion, (N)x(N), upper Hessenberg, LAPACK dgeev, O(N^3).
    """
    x, _ = gll(N)
    c = np.polyfit(x, f, N)
    r = np.roots(c)
    r = r[np.abs(r.imag) < 1e-8].real
    return np.sort(r[(r >= -1.0) & (r <= 1.0)])


def ref_roots_brent(N, f, nsample=2001):
    """REFERENCE: sample the interpolant densely, bracket sign changes, Brent.

    This is 'what everyone does' -- and it is exactly the subgrid-sampling logic
    marching cubes uses, so its failure mode (two roots inside one subinterval
    cancel and are both missed) is the one this project is trying to remove.
    No eigenproblem.  Cost O(nsample * N) for the sampling (barycentric eval).
    """
    from scipy.optimize import brentq
    x, _ = gll(N)
    b = bary_weights(x)
    xs = np.linspace(-1, 1, nsample)

    def ev(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        num = np.zeros_like(t)
        den = np.zeros_like(t)
        out = np.full(t.shape, np.nan)
        for j in range(len(x)):
            d = t - x[j]
            hit = np.abs(d) < 1e-14
            out[hit] = f[j]
            with np.errstate(divide="ignore", invalid="ignore"):
                q = b[j] / d
            q[hit] = 0.0
            num += q * f[j]
            den += q
        r = num / den
        r[~np.isnan(out)] = out[~np.isnan(out)]
        return r

    fs = ev(xs)
    roots = []
    for i in range(nsample - 1):
        if fs[i] == 0.0:
            roots.append(xs[i])
        elif fs[i] * fs[i + 1] < 0:
            roots.append(brentq(lambda t: float(ev(t)[0]), xs[i], xs[i + 1],
                                xtol=1e-15, rtol=8.9e-16))
    if fs[-1] == 0.0:
        roots.append(xs[-1])
    return np.sort(np.array(roots))


# --------------------------------------------------------------------------
# Timing
# --------------------------------------------------------------------------
def timeit_us(fn, reps=None, min_time=0.05):
    """Median microseconds per call of fn().  Auto-sizes reps to ~min_time seconds."""
    fn()                                      # warm up (JIT-free, but caches/LAPACK workspaces)
    if reps is None:
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        reps = max(5, min(20000, int(min_time / max(dt, 1e-9))))
    samples = []
    batch = max(1, reps // 7)
    for _ in range(7):
        t0 = time.perf_counter()
        for _ in range(batch):
            fn()
        samples.append((time.perf_counter() - t0) / batch)
    return float(np.median(samples)) * 1e6
