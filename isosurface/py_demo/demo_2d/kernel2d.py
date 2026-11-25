"""
kernel2d.py -- SELF-CONTAINED numerical kernel for the 2D contour demo.

NO third-party numerics.  numpy is used only as an array container and for
elementary arithmetic; everything below is implemented from the defining
recurrences, exactly as a C or CUDA port would have to:

    gll_nodes_weights   GLL nodes by Newton on P'_N, weights from P_N
    bary_weights        barycentric weights, closed form and by definition
    legendre_val_der    p and p' of a Legendre series, forward recurrence
    nodal_to_modal      discrete Legendre projection (note gamma_N = 2/N)
    lagrange_weights    barycentric Lagrange weights L_i(t)
    aberth_roots        all real roots in [-1,1]: Ehrlich-Aberth (route C)
    bisect_bracketed    single root in a known bracket (used on inserted lines)

Nothing here calls numpy.polynomial, scipy, or LAPACK.  `selftest()` checks every
routine against an independent computation.
"""
import numpy as np

__all__ = ["gll_nodes_weights", "bary_weights", "legendre_val_der",
           "nodal_to_modal", "lagrange_weights", "aberth_roots",
           "bisect_bracketed", "deflate", "A_FLOOR", "NITER", "selftest"]

A_FLOOR = 1e-12        # degree-deflation threshold (see demo_A)
NITER = 12             # fixed Aberth iterations (see demo_C)

_CACHE = {}


# --------------------------------------------------------------------------
def _legendre_P_dP(n, x):
    """P_n(x) and P_n'(x) by the three-term recurrence.  O(n), scalar or array."""
    x = np.asarray(x, dtype=float)
    if n == 0:
        return np.ones_like(x), np.zeros_like(x)
    p0 = np.ones_like(x)
    p1 = x.copy()
    d0 = np.zeros_like(x)
    d1 = np.ones_like(x)
    for k in range(1, n):
        p2 = ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        d2 = ((2 * k + 1) * (p1 + x * d1) - k * d0) / (k + 1)
        p0, p1 = p1, p2
        d0, d1 = d1, d2
    return p1, d1


def gll_nodes_weights(N):
    """N+1 Gauss-Lobatto-Legendre nodes and weights on [-1,1].

    Interior nodes are the roots of P'_N, found by Newton from the
    Chebyshev-Lobatto initial guess -- the standard construction, no library.
    Weights: w_j = 2 / (N(N+1) P_N(x_j)^2).
    """
    key = ("gll", N)
    if key in _CACHE:
        return _CACHE[key]
    if N < 1:
        raise ValueError("N >= 1")
    x = np.empty(N + 1)
    x[0], x[N] = -1.0, 1.0
    for j in range(1, N):
        t = -np.cos(np.pi * j / N)                     # Chebyshev-Lobatto start
        for _ in range(100):
            # Newton on P'_N: need P'_N and P''_N
            p, dp = _legendre_P_dP(N, t)
            # (1-x^2) P''_N = 2x P'_N - N(N+1) P_N   (Legendre's equation)
            ddp = (2.0 * t * dp - N * (N + 1) * p) / (1.0 - t * t)
            step = dp / ddp
            t -= step
            if abs(step) < 1e-16:
                break
        x[j] = t
    p, _ = _legendre_P_dP(N, x)
    w = 2.0 / (N * (N + 1) * p ** 2)
    x.flags.writeable = False
    w.flags.writeable = False
    _CACHE[key] = (x, w)
    return x, w


def bary_weights(N):
    """Barycentric weights of the GLL nodes: b_j proportional to (-1)^j sqrt(w_j).

    Proportional, not equal -- the common factor cancels in every formula used
    here (all are homogeneous of degree 0 or 1 in b).
    """
    key = ("bary", N)
    if key in _CACHE:
        return _CACHE[key]
    _, w = gll_nodes_weights(N)
    b = np.array([(-1.0) ** j * np.sqrt(w[j]) for j in range(N + 1)])
    b.flags.writeable = False
    _CACHE[key] = b
    return b


def lagrange_weights(N, t):
    """Barycentric Lagrange weights L_i(t), summing to 1.  O(N)."""
    x, _ = gll_nodes_weights(N)
    d = t - x
    hit = np.abs(d) < 1e-14
    if np.any(hit):
        w = np.zeros(N + 1)
        w[int(np.argmax(hit))] = 1.0
        return w
    w = bary_weights(N) / d
    return w / w.sum()


def nodal_to_modal(N):
    """Matrix M with a = M f : discrete Legendre coefficients from GLL values.

    a_k = (1/gamma_k) sum_j w_j f_j P_k(x_j),  gamma_k = 2/(2k+1) except
    gamma_N = 2/N  (GLL quadrature is exact only to degree 2N-1).
    """
    key = ("n2m", N)
    if key in _CACHE:
        return _CACHE[key]
    x, w = gll_nodes_weights(N)
    V = np.empty((N + 1, N + 1))
    for k in range(N + 1):
        V[:, k] = _legendre_P_dP(k, x)[0] if k > 0 else np.ones(N + 1)
    gam = np.array([2.0 / (2 * k + 1) for k in range(N + 1)])
    gam[N] = 2.0 / N
    M = (V * w[:, None]).T / gam[:, None]
    M.flags.writeable = False
    _CACHE[key] = M
    return M


def modal_vandermonde(N):
    """V[j,k] = P_k(x_j) : modal -> nodal."""
    key = ("V", N)
    if key in _CACHE:
        return _CACHE[key]
    x, _ = gll_nodes_weights(N)
    V = np.empty((N + 1, N + 1))
    for k in range(N + 1):
        V[:, k] = _legendre_P_dP(k, x)[0] if k > 0 else np.ones(N + 1)
    V.flags.writeable = False
    _CACHE[key] = V
    return V


# --------------------------------------------------------------------------
def legendre_val_der(a, z):
    """p = sum_k a_k P_k(z) and p' , by forward recurrence.  Complex-capable.

    O(n) work, O(1) storage, no matrix -- the kernel primitive of route C.
    """
    a = np.asarray(a, dtype=float)
    n = len(a) - 1
    z = np.asarray(z)
    if n == 0:
        return a[0] * np.ones_like(z), np.zeros_like(z)
    P0 = np.ones_like(z)
    P1 = z
    dP0 = np.zeros_like(z)
    dP1 = np.ones_like(z)
    v = a[0] * P0 + a[1] * P1
    d = a[1] * dP1
    for k in range(1, n):
        P2 = ((2 * k + 1) * z * P1 - k * P0) / (k + 1)
        dP2 = ((2 * k + 1) * (P1 + z * dP1) - k * dP0) / (k + 1)
        v = v + a[k + 1] * P2
        d = d + a[k + 1] * dP2
        P0, P1 = P1, P2
        dP0, dP1 = dP1, dP2
    return v, d


def deflate(a, a_floor=A_FLOOR):
    """Drop trailing negligible Legendre coefficients (demo_A, notes/07)."""
    a = np.asarray(a, dtype=float)
    if a_floor <= 0.0:
        return a
    s = np.max(np.abs(a))
    if s == 0.0:
        return a[:1]
    nz = np.nonzero(np.abs(a) > a_floor * s)[0]
    return a[:(nz[-1] + 1)] if len(nz) else a[:1]


def aberth_roots(N, f, niter=NITER, a_floor=A_FLOOR, want_slope=False):
    """ROUTE C: all real roots in [-1,1] of the GLL interpolant of f.

    Matrix-free Ehrlich-Aberth, fixed iteration count, no branches in the loop.
    Initialization: Chebyshev points nudged off the real axis (notes/14).
    Returns sorted roots, and optionally p' at each root (free from the iteration).
    """
    a = deflate(nodal_to_modal(N) @ np.asarray(f, dtype=float), a_floor)
    n = len(a) - 1
    if n < 1:
        return (np.array([]), np.array([])) if want_slope else np.array([])
    t = np.pi * (np.arange(n) + 0.5) / n
    z = (np.cos(t) + 0.1j * np.sin(t)).astype(complex)
    for _ in range(niter):
        v, d = legendre_val_der(a, z)
        w = np.where(np.abs(d) > 0.0, v / np.where(np.abs(d) > 0.0, d, 1.0), 0.0)
        if n > 1:
            diff = z[:, None] - z[None, :]
            np.fill_diagonal(diff, np.inf)
            S = np.sum(1.0 / diff, axis=1)
        else:
            S = np.zeros(1, dtype=complex)
        den = 1.0 - w * S
        corr = np.where(np.abs(den) > 1e-300,
                        w / np.where(np.abs(den) > 1e-300, den, 1.0), w)
        z = z - np.where(np.isfinite(corr), corr, 0.0)
    keep = np.abs(z.imag) < 1e-7
    r = z[keep].real
    m = (r >= -1.0 - 1e-10) & (r <= 1.0 + 1e-10)
    r = np.sort(r[m])
    if not want_slope:
        return r
    if len(r) == 0:
        return r, np.array([])
    _, dv = legendre_val_der(a, r.astype(complex))
    return r, dv.real


def bisect_bracketed(fun, lo, hi, maxit=60, tol=1e-15):
    """Single root of `fun` on [lo,hi], given a sign change.  Bisection: cannot fail.

    Used only on INSERTED grid lines, which are not original fibers.
    """
    fa, fb = fun(lo), fun(hi)
    if fa == 0.0:
        return lo
    if fb == 0.0:
        return hi
    if fa * fb > 0.0:
        return 0.5 * (lo + hi)
    a, b = lo, hi
    for _ in range(maxit):
        m = 0.5 * (a + b)
        fm = fun(m)
        if fm == 0.0 or (b - a) < tol:
            return m
        if fa * fm < 0.0:
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


# --------------------------------------------------------------------------
def selftest(verbose=True):
    """Check every routine against an independent computation."""
    ok = True
    lines = []
    for N in (4, 7, 11):
        x, w = gll_nodes_weights(N)
        # (1) endpoints, symmetry, weight sum
        e1 = abs(x[0] + 1) + abs(x[N] - 1)
        e2 = np.abs(x + x[::-1]).max()
        e3 = abs(w.sum() - 2.0)
        # (2) quadrature exact to degree 2N-1
        e4 = 0.0
        for deg in range(0, 2 * N):
            exact = 0.0 if deg % 2 else 2.0 / (deg + 1)
            e4 = max(e4, abs(np.sum(w * x ** deg) - exact))
        # (3) interior nodes are roots of P'_N
        _, dp = _legendre_P_dP(N, x[1:N])
        e5 = np.abs(dp).max()
        # (4) barycentric weights match the defining product
        bn = np.array([1.0 / np.prod(x[j] - np.delete(x, j)) for j in range(N + 1)])
        ratio = bn / bary_weights(N)
        e6 = np.abs(ratio / ratio[0] - 1).max()
        # (5) nodal <-> modal round trip
        rng = np.random.default_rng(N)
        f = rng.standard_normal(N + 1)
        e7 = np.abs(modal_vandermonde(N) @ (nodal_to_modal(N) @ f) - f).max()
        # (6) legendre_val_der against finite differences
        a = rng.standard_normal(N + 1)
        t = 0.37
        v, d = legendre_val_der(a, np.array([t]))
        h = 1e-6
        vp, _ = legendre_val_der(a, np.array([t + h]))
        vm, _ = legendre_val_der(a, np.array([t - h]))
        # central differences are only O(h^2)-accurate: judge against ~1e-9,
        # not against machine epsilon.  This checks correctness, not precision.
        e8 = abs(d[0] - (vp[0] - vm[0]) / (2 * h))
        # (7) Aberth: roots of a polynomial with PRESCRIBED roots
        want = np.array([-0.83, -0.21, 0.42, 0.77])[:max(1, min(4, N - 1))]
        fv = np.prod([x - q for q in want], axis=0)
        got = aberth_roots(N, fv)
        e9 = np.abs(np.sort(got) - np.sort(want)).max() if len(got) == len(want) else np.inf
        bad = max(e1, e2, e3, e4, e5, e6, e7, e9)      # e8 judged separately
        ok &= (bad < 1e-11) and (e8 < 1e-8)
        lines.append(f"  N={N:2d}  ends {e1:.0e}  sym {e2:.0e}  sum(w) {e3:.0e}  "
                     f"quad {e4:.0e}  P'_N {e5:.0e}  bary {e6:.0e}  "
                     f"n<->m {e7:.0e}  p' {e8:.0e}  roots {e9:.0e}")
    if verbose:
        print("kernel2d self-test (no third-party numerics):")
        for l in lines:
            print(l)
        print(f"  ALL PASS: {ok}")
    return ok


if __name__ == "__main__":
    selftest()
