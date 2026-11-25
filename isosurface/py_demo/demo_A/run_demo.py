#!/usr/bin/env python3
"""
run_demo.py -- ROUTE A demo driver (modal Legendre comrade matrix + degree deflation).

Deflation is ON by default (a_floor = 1e-12); --plain turns it off to expose the
original failure mode.  See route_a.py for why, and the 'deflation' study below.

Source mesh is always a SINGLE reference element [-1,1] with GLL nodes of order N.

Usage
    python3 run_demo.py                       # all studies, orders 3..15
    python3 run_demo.py --study convergence   # one study
    python3 run_demo.py -N 8                  # single-order detail dump
    python3 run_demo.py -N 8 --field sine3 --level 0.5
    python3 run_demo.py --orders 5,7,9,11     # custom order list
    python3 run_demo.py --no-fig              # skip figure generation
    python3 run_demo.py --plain               # DISABLE deflation (original route A)
    python3 run_demo.py --a-floor 1e-10       # custom deflation threshold

Studies
    convergence  root error vs N for 7 analytic test fields (+ figure)
    cost         us/solve vs N, nnz, flop model (+ figure)
    degeneracy   what happens as the top mode a_N -> 0   <-- deflation's purpose
    deflation    threshold sweep: plain vs deflated, and how often it fires
    reference    vs monomial companion and dense-sample+Brent
    all          every study (default)

See route_a.py for the full math, the eigenproblem characterisation, and the
complexity analysis.
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)

from gllbase import gll, TESTS, match_roots, nodal_to_modal      # noqa: E402
import studies                                                    # noqa: E402
from route_a import (roots_route_a, comrade_matrix, comrade_nnz,  # noqa: E402
                     modal_coeffs, deflate_degree, A_FLOOR_DEFAULT)

LABEL = "route A + deflation (modal Legendre comrade)"


def flop_model_a(N):
    """dgehrd (10/3)N^3 + dhseqr ~10 N^3, plus 2(N+1)^2 for the nodal->modal apply."""
    return (10.0 / 3.0) * N ** 3 + 10.0 * N ** 3 + 2.0 * (N + 1) ** 2


def study_deflation(a_floor, orders=(3, 5, 6, 8, 11, 15), verbose=True):
    """Plain route A vs deflated route A vs route B, plus how often deflation fires."""
    sys.path.insert(0, os.path.join(HERE, "..", "demo_B"))
    from route_b import roots_route_b
    from gllbase import TESTS, match_roots

    rows = []
    for k in TESTS:
        t = TESTS[k]
        for N in orders:
            x, _ = gll(N)
            f = t.values(x)
            try:
                rp = roots_route_a(N, f, a_floor=0.0)
                ep = match_roots(rp, t.exact())["max_err"]
            except Exception:
                ep = np.inf
            rd, info = roots_route_a(N, f, a_floor=a_floor, return_info=True)
            ed = match_roots(rd, t.exact())["max_err"]
            eb = match_roots(roots_route_b(N, f), t.exact())["max_err"]
            rows.append(dict(field=k, N=N, plain=ep, defl=ed, b=eb,
                             fired=info["deflated"], deg=info["degree"],
                             ratio=info["a_N_ratio"]))
    if verbose:
        print(f"\n=== route A: deflation (a_floor = {a_floor:g}) ===\n")
        print("  Deflation drops trailing Legendre coefficients below a_floor*max|a|")
        print("  BEFORE building the comrade matrix, so the division in the last row is")
        print("  by a genuine leading coefficient instead of by roundoff.\n")
        print("  field       N   |a_N|/max   deg  fired |  plain A     deflated    route B")
        print("  " + "-" * 76)
        for r in rows:
            mark = " * " if (np.isfinite(r["plain"]) and r["defl"] < r["plain"] * 1e-3) else "   "
            print(f"  {r['field']:<10s} {r['N']:2d}   {r['ratio']:9.2e}  {r['deg']:3d}  "
                  f"{'yes' if r['fired'] else ' - ':>5s} | {r['plain']:9.2e} {mark}"
                  f"{r['defl']:9.2e}   {r['b']:9.2e}")
        nf = sum(r["fired"] for r in rows)
        better = sum(1 for r in rows
                     if np.isfinite(r["plain"]) and r["defl"] < r["plain"] * 1e-3)
        print(f"\n  deflation fired on {nf}/{len(rows)} cases; "
              f"materially improved the answer on {better} (* rows)")
        print("  threshold justification: see route_a.py docstring (20,000-draw sweep)")
    return rows


def detail(N, field, level=None, a_floor=A_FLOOR_DEFAULT):
    """Single-order dump: coefficients, deflation decision, matrix, spectrum, roots."""
    t = TESTS[field]
    lv = t.level if level is None else level
    x, w = gll(N)
    f = t.raw(x) - lv
    print(f"\n=== ROUTE A detail: N={N}, field='{field}' ({t.desc}), level={lv:g},"
          f" a_floor={a_floor:g} ===\n")
    print(f"  GLL nodes   x_j = {np.array2string(np.asarray(x), precision=6, max_line_width=100)}")
    print(f"  GLL weights w_j = {np.array2string(np.asarray(w), precision=6, max_line_width=100)}")
    print(f"  nodal f_j - c   = {np.array2string(f, precision=6, max_line_width=100)}")

    a = modal_coeffs(N, f)
    print(f"\n  Legendre coefficients a_k (a = M f):")
    for k, ak in enumerate(a):
        bar = "#" * int(np.clip(np.log10(abs(ak) + 1e-300) + 17, 0, 40))
        print(f"    a_{k:<2d} = {ak:+.12e}   {bar}")
    ratio = abs(a[N]) / np.max(np.abs(a))
    print(f"\n  leading coefficient a_N = {a[N]:+.6e}   |a_N|/max|a| = {ratio:.3e}")

    ad = deflate_degree(a, a_floor)
    d = len(ad) - 1
    if d < N:
        print(f"  DEFLATION FIRED: |a_N|/max|a| <= {a_floor:g} -> true degree {d}, "
              f"solving a {d}x{d} problem instead of {N}x{N}")
        print(f"    (plain route A would divide by {a[N]:+.3e} here)")
        try:
            rp = roots_route_a(N, f, a_floor=0.0)
            print(f"    plain route A would return: {np.array2string(rp, precision=8)}")
        except Exception as e:
            print(f"    plain route A would FAIL: {str(e)[:60]}")
    else:
        print(f"  deflation did not fire (degree {d} = N); dividing by a genuine a_N")

    C = comrade_matrix(ad)
    print(f"\n  comrade matrix C: {d}x{d}, nnz={np.count_nonzero(C)} / {d*d} "
          f"({100*np.count_nonzero(C)/max(d*d,1):.1f}% dense), tridiagonal + full last row")
    if d <= 10:
        print("  " + np.array2string(C, precision=4, suppress_small=True,
                                     max_line_width=140).replace("\n", "\n  "))
    print(f"  cond(C) = {np.linalg.cond(C):.3e}")

    r, lam = roots_route_a(N, f, return_all=True, a_floor=a_floor)
    print(f"\n  spectrum ({len(lam)} eigenvalues):")
    for v in sorted(lam, key=lambda z: (abs(z.imag) > 1e-8, z.real)):
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
    ap = argparse.ArgumentParser(description="Route A demo: modal Legendre comrade matrix",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("--study", default="all",
                    choices=["all", "convergence", "cost", "degeneracy", "deflation",
                             "reference"])
    ap.add_argument("-N", type=int, default=None, help="single-order detail dump")
    ap.add_argument("--field", default="sine3", choices=sorted(TESTS.keys()))
    ap.add_argument("--level", type=float, default=None, help="override the isolevel")
    ap.add_argument("--orders", default=None, help="comma list, e.g. 5,7,9")
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--figdir", default=os.path.join(HERE, "..", "figs"))
    ap.add_argument("--plain", action="store_true",
                    help="disable degree deflation (original, fragile route A)")
    ap.add_argument("--a-floor", type=float, default=None,
                    help=f"deflation threshold (default {A_FLOOR_DEFAULT:g})")
    args = ap.parse_args()

    a_floor = 0.0 if args.plain else (args.a_floor if args.a_floor is not None
                                      else A_FLOOR_DEFAULT)

    if args.N is not None:
        detail(args.N, args.field, args.level, a_floor)
        return

    orders = ([int(s) for s in args.orders.split(",")] if args.orders
              else studies.ORDERS)
    route = lambda N, f: roots_route_a(N, f, a_floor=a_floor)
    os.makedirs(args.figdir, exist_ok=True)

    label = LABEL if a_floor > 0 else "route A PLAIN (no deflation)"
    print("=" * 80)
    print(f"  {label}   [a_floor = {a_floor:g}]")
    print("  single GLL reference element on [-1,1], orders " +
          ",".join(map(str, orders)))
    print("=" * 80)

    if args.study in ("all", "convergence"):
        res = studies.study_convergence(route, label, orders)
        if not args.no_fig:
            p = studies.plot_convergence(res, label,
                                         os.path.join(args.figdir, "A_convergence.png"),
                                         orders)
            print(f"\n  figure -> {os.path.relpath(p, HERE)}")

    if args.study in ("all", "cost"):
        rows = studies.study_cost(route, label, orders,
                                  nnz_fn=lambda N: (comrade_nnz(N)[0], comrade_nnz(N)[1]),
                                  flop_model=flop_model_a)
        print("\n  eigensolver: LAPACK dgeev on an NxN non-symmetric matrix")
        print("               (tridiagonal + full last row, nnz = 3N-3)")
        print("               = dgebal (balance) + dgehrd (Hessenberg) + dhseqr (Francis QR)")
        print("               balancing is essential -- see route_a.py")
        if not args.no_fig:
            p = studies.plot_cost(rows, label, os.path.join(args.figdir, "A_cost.png"),
                                  "(10/3+10) N^3")
            print(f"  figure -> {os.path.relpath(p, HERE)}")

    if args.study in ("all", "deflation"):
        study_deflation(a_floor if a_floor > 0 else A_FLOOR_DEFAULT)

    if args.study in ("all", "degeneracy"):
        studies.study_degeneracy(route, label)
        if a_floor > 0:
            print("\n  With deflation ON the sweep is flat: the degree is trimmed before the")
            print("  comrade matrix is built, so the division is never by roundoff.")
            print("  Re-run with --plain to see the original failure.")
        else:
            print("\n  PLAIN route A divides by a_N (the spill-over row) and degrades as it")
            print("  shrinks.  Deflation (the default) removes this; re-run without --plain.")

    if args.study in ("all", "reference"):
        studies.study_reference(route, label)


if __name__ == "__main__":
    main()
