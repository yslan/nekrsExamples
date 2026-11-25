#!/usr/bin/env python3
"""
run_demo.py -- ROUTE B demo driver (nodal barycentric arrowhead pencil).

Source mesh is always a SINGLE reference element [-1,1] with GLL nodes of order N.

Usage
    python3 run_demo.py                       # all studies, orders 3..15
    python3 run_demo.py --study convergence   # one study
    python3 run_demo.py -N 8                  # single-order detail dump
    python3 run_demo.py -N 8 --field sine3 --level 0.5
    python3 run_demo.py --orders 5,7,9,11     # custom order list
    python3 run_demo.py --no-fig              # skip figure generation

Studies
    convergence  root error vs N for 8 analytic test fields (+ figure)
    cost         us/solve vs N, nnz, flop model (+ figure)
    degeneracy   what happens as the top mode a_N -> 0   <-- route B's strength
    identity     verify det(A - lam B) = -p(lam) pointwise
    reference    vs monomial companion and dense-sample+Brent
    all          every study (default)

See route_b.py for the full math, the eigenproblem characterisation, and the
complexity analysis.
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)

from gllbase import gll, TESTS, match_roots, nodal_to_modal, gll_bary_weights  # noqa: E402
import studies                                                    # noqa: E402
from route_b import roots_route_b, arrowhead_pencil, pencil_nnz   # noqa: E402

LABEL = "route B (nodal barycentric pencil)"


def flop_model_b(N):
    """dgghrd ~(20/3)n^3 + dhgeqz ~30 n^3 (eigenvalues only), n = N+2."""
    n = N + 2.0
    return (20.0 / 3.0) * n ** 3 + 30.0 * n ** 3


def study_identity(orders=(4, 6, 8, 11, 15), verbose=True):
    """Check det(A - lam*B) = -kappa * p(lam) pointwise: the pencil IS the interpolant.

    NOTE on the constant kappa.  gll_bary_weights() returns the closed form
    (-1)^j sqrt(w_j), which is the true barycentric weight set only UP TO A COMMON
    SCALE kappa(N) (measured: 2.77 at N=4, -3241 at N=15).  The pencil is
    homogeneous of degree 1 in b, so every root is unaffected and the scale is
    irrelevant in production -- but the determinant identity carries it, so the
    ratio det/p is the CONSTANT -kappa rather than -1.  The meaningful test is
    therefore that the ratio does not vary with lambda; that is what is reported.
    """
    from numpy.polynomial import legendre as L
    rng = np.random.default_rng(7)
    rows = []
    for N in orders:
        x, _ = gll(N)
        f = rng.standard_normal(N + 1)
        a = nodal_to_modal(N) @ f
        A, Bm = arrowhead_pencil(N, f)
        for lam in (-0.83, -0.21, 0.0, 0.47, 0.96, 1.6):
            d = np.linalg.det(A - lam * Bm)
            p = L.legval(lam, a)
            rows.append(dict(N=N, lam=lam, det=d, p=p,
                             ratio=(d / p if p != 0 else np.nan)))
        rn = np.array([r["ratio"] for r in rows if r["N"] == N])
        for r in rows:
            if r["N"] == N:
                r["norm"] = r["ratio"] / np.median(rn)      # scale-free check
    if verbose:
        print(f"\n=== {LABEL}: determinant identity  det(A - lam B) = -p(lam) ===\n")
        print("    N     lambda        det(A-lam B)          p(lambda)          det/p"
              "          normalized")
        print("  " + "-" * 88)
        for r in rows:
            print(f"  {r['N']:3d}   {r['lam']:+6.2f}   {r['det']:+.10e}   "
                  f"{r['p']:+.10e}   {r['ratio']:+.6e}   {r['norm']:+.12f}")
        nn = np.array([r["norm"] for r in rows])
        print(f"\n  max |normalized - 1| = {np.abs(nn - 1).max():.3e}")
        print("  det/p is CONSTANT in lambda (to rounding) at every N => the pencil's")
        print("  characteristic polynomial is p up to the scale kappa(N) carried by the")
        print("  closed-form barycentric weights.  Roots are invariant to kappa.")
    return rows


def detail(N, field, level=None):
    """Single-order dump: pencil, spectrum, infinite eigenvalues, roots."""
    t = TESTS[field]
    lv = t.level if level is None else level
    x, w = gll(N)
    f = t.raw(x) - lv
    b = gll_bary_weights(N)
    print(f"\n=== ROUTE B detail: N={N}, field='{field}' ({t.desc}), level={lv:g} ===\n")
    print(f"  GLL nodes   x_j = {np.array2string(np.asarray(x), precision=6, max_line_width=100)}")
    print(f"  bary wts    b_j = {np.array2string(np.asarray(b), precision=6, max_line_width=100)}")
    print(f"     spread max|b|/min|b| = {np.abs(b).max()/np.abs(b).min():.3f}"
          "   (O(N) on GLL; O(4^N/N) on equispaced -- see route_b.py)")
    print(f"  nodal f_j - c   = {np.array2string(f, precision=6, max_line_width=100)}")

    a = nodal_to_modal(N) @ f
    print(f"\n  (for reference only -- route B never forms these)")
    print(f"  leading Legendre coefficient a_N = {a[N]:+.6e}   "
          f"|a_N|/max|a| = {abs(a[N])/np.max(np.abs(a)):.3e}")

    A, Bm = arrowhead_pencil(N, f)
    nzA, nzB, dn = pencil_nnz(N)
    print(f"\n  pencil (A,B): {N+2}x{N+2}   nnz(A)={nzA} / {dn} "
          f"({100*nzA/dn:.1f}% dense, ARROWHEAD), nnz(B)={nzB} (diagonal, singular)")
    if N <= 8:
        print("\n  A =")
        print("  " + np.array2string(A, precision=4, suppress_small=True,
                                     max_line_width=140).replace("\n", "\n  "))
        print("  diag(B) =", np.array2string(np.diag(Bm), precision=1))

    r, info = roots_route_b(N, f, return_info=True)
    print(f"\n  eigenvalues: {info['n_finite']} finite, {info['n_infinite']} infinite")
    print(f"     (floor is 1 infinite from the singular B entry; extras = degree deficiency)")
    for v in sorted(info["spectrum"], key=lambda z: (abs(z.imag) > 1e-8, z.real)):
        tag = "  <- real root in [-1,1]" if (abs(v.imag) < 1e-8 and abs(v.real) <= 1) else ""
        print(f"    {v.real:+.12f} {v.imag:+.12f}i{tag}")

    ex = t.exact() if level is None else None
    print(f"\n  roots found ({len(r)}): {np.array2string(r, precision=12)}")
    if ex is not None:
        m = match_roots(r, ex)
        print(f"  exact roots  ({len(ex)}): {np.array2string(ex, precision=12)}")
        print(f"  max error = {m['max_err']:.3e}   missed = {m['n_missed']}   "
              f"spurious = {m['n_spurious']}")


def main():
    ap = argparse.ArgumentParser(description="Route B demo: nodal barycentric arrowhead pencil",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("--study", default="all",
                    choices=["all", "convergence", "cost", "degeneracy",
                             "identity", "reference"])
    ap.add_argument("-N", type=int, default=None, help="single-order detail dump")
    ap.add_argument("--field", default="sine3", choices=sorted(TESTS.keys()))
    ap.add_argument("--level", type=float, default=None, help="override the isolevel")
    ap.add_argument("--orders", default=None, help="comma list, e.g. 5,7,9")
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--figdir", default=os.path.join(HERE, "..", "figs"))
    args = ap.parse_args()

    if args.N is not None:
        detail(args.N, args.field, args.level)
        return

    orders = ([int(s) for s in args.orders.split(",")] if args.orders
              else studies.ORDERS)
    route = roots_route_b
    os.makedirs(args.figdir, exist_ok=True)

    print("=" * 80)
    print(f"  {LABEL}")
    print("  single GLL reference element on [-1,1], orders " +
          ",".join(map(str, orders)))
    print("=" * 80)

    if args.study in ("all", "identity"):
        study_identity()

    if args.study in ("all", "convergence"):
        res = studies.study_convergence(route, LABEL, orders)
        if not args.no_fig:
            p = studies.plot_convergence(res, LABEL,
                                         os.path.join(args.figdir, "B_convergence.png"),
                                         orders)
            print(f"\n  figure -> {os.path.relpath(p, HERE)}")

    if args.study in ("all", "cost"):
        rows = studies.study_cost(route, LABEL, orders,
                                  nnz_fn=lambda N: (pencil_nnz(N)[0], pencil_nnz(N)[2]),
                                  flop_model=flop_model_b)
        print("\n  eigensolver: LAPACK dggev on a GENERALIZED problem, (N+2)x(N+2)")
        print("               A = arrowhead (nnz = 3N+3), B = diag(0,1,..,1) SINGULAR")
        print("               = dggbal + dgghrd (Hessenberg-triangular) + dhgeqz (QZ)")
        print("               infinite eigenvalues filtered by the RELATIVE test")
        print("               |beta| > tol*|alpha| -- see route_b.py")
        if not args.no_fig:
            p = studies.plot_cost(rows, LABEL, os.path.join(args.figdir, "B_cost.png"),
                                  "(20/3+30)(N+2)^3")
            print(f"  figure -> {os.path.relpath(p, HERE)}")

    if args.study in ("all", "degeneracy"):
        studies.study_degeneracy(route, LABEL)
        print("\n  Route B never divides by a_N.  Degree deficiency shows up as an EXTRA")
        print("  infinite eigenvalue (n_infinite rises above its floor of 1), which is")
        print("  information, not a failure.  Compare demo_A.")

    if args.study in ("all", "reference"):
        studies.study_reference(route, LABEL)


if __name__ == "__main__":
    main()
