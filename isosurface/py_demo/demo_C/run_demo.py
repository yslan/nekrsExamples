#!/usr/bin/env python3
"""
run_demo.py -- ROUTE C demo driver (matrix-free Ehrlich-Aberth on the Legendre series).

Source mesh is always a SINGLE reference element [-1,1] with GLL nodes of order N.
Route C is the GPU candidate: one fiber per thread, fixed iteration count, no matrix.

Usage
    python3 run_demo.py                       # all studies, orders 3..15
    python3 run_demo.py --study iterations    # how many iterations are enough?
    python3 run_demo.py --orders 5,7,9,11     # custom order list
    python3 run_demo.py --niter 10            # change the fixed iteration count
    python3 run_demo.py --no-fig              # skip figure generation

    python3 run_demo.py -N 8                             # single-order detail dump
    python3 run_demo.py -N 8 --field sine3               # pick a test field
    python3 run_demo.py -N 8 --field gauss --level 0.7   # override the isolevel

Studies
    convergence  root error vs N for 8 analytic test fields (+ figure)
    iterations   correction history + accuracy vs niter (+ figure)   <-- route C's core claim
    cost         us/solve vs N, per-fiber state, flop model (+ figure)
    gpu          per-thread register budget and work, vs routes A and B
    reference    vs monomial companion and dense-sample+Brent
    all          every study (default)

See route_c.py for the math, the iteration-count justification, and why a
hand-written QR was measured and rejected.
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "demo_A"))
sys.path.insert(0, os.path.join(HERE, "..", "demo_B"))
sys.path.insert(0, HERE)

from gllbase import gll, TESTS, match_roots, nodal_to_modal      # noqa: E402
import studies                                                    # noqa: E402
from route_c import (roots_route_c, legendre_val_der, aberth_init,  # noqa: E402
                     aberth_iterate, state_doubles, NITER_DEFAULT)
from route_a import deflate_degree, A_FLOOR_DEFAULT               # noqa: E402

LABEL = "route C (matrix-free Ehrlich-Aberth)"


def flop_model_c(N, niter=NITER_DEFAULT):
    """niter * [ n evals O(n) + n repulsion sums O(n) ] = O(niter n^2)."""
    return 2.0 * niter * N * N


def study_iterations(orders=(5, 8, 11, 15), niters=(4, 6, 8, 10, 12, 16, 24),
                     verbose=True):
    """Route C's central claim: a FIXED iteration count suffices."""
    from route_b import roots_route_b

    hist_rows = []
    for k in ("poly5", "gauss", "sine2", "sine3", "tanh"):
        t = TESTS[k]
        N = 8
        x, _ = gll(N)
        _, h = roots_route_c(N, t.values(x), niter=16, track=True)
        hist_rows.append((k, h))

    acc = {}
    for k in ("poly5", "gauss", "tanh", "sine2", "sine3", "runge"):
        t = TESTS[k]
        for N in orders:
            x, _ = gll(N)
            f = t.values(x)
            eb = match_roots(roots_route_b(N, f), t.exact())["max_err"]
            row = []
            for it in niters:
                r = roots_route_c(N, f, niter=it)
                row.append(match_roots(r, t.exact())["max_err"])
            acc[(k, N)] = (row, eb)

    if verbose:
        print(f"\n=== {LABEL}: iteration count ===\n")
        print("  Correction magnitude per iteration (N=8) -- cubic, then stalls at eps:\n")
        print("    field    " + "  ".join(f"it{i+1:<2d}" for i in range(12)))
        for k, h in hist_rows:
            print(f"    {k:8s} " + "  ".join(f"{v:.0e}" for v in h[:12]))
        print("\n  Root error vs fixed iteration count (route B for reference):\n")
        print("    field      N  " + "  ".join(f"it={i:<3d}" for i in niters) + "   routeB")
        for (k, N), (row, eb) in acc.items():
            print(f"    {k:8s} {N:3d}  " + "  ".join(f"{v:.0e}" for v in row) +
                  f"   {eb:.0e}")
        print(f"\n  DEFAULT niter = {NITER_DEFAULT}.  Measured over 20,000 random fibers")
        print("  (scratch/23): niter=6 -> 57 root-count mismatches vs route B, niter=8 -> 2,")
        print("  niter>=10 -> 0 with max difference 2.0e-13.  Convergence is cubic and then")
        print("  STALLS at machine precision, so overshooting the count is free.")
    return hist_rows, acc


def study_gpu(orders=(5, 8, 11, 15), verbose=True):
    """Per-thread budget: route C vs routes A and B."""
    rows = []
    for N in orders:
        rows.append(dict(N=N, c=state_doubles(N), a=N * N, b=2 * (N + 2) ** 2,
                         fib=3 * (N + 1) ** 2,
                         flops_c=flop_model_c(N), flops_qr=10.0 * N ** 3))
    if verbose:
        print(f"\n=== {LABEL}: GPU budget (one fiber per thread) ===\n")
        print("  Per-thread state, in doubles.  A CUDA thread has ~255 32-bit registers")
        print("  ~ 128 doubles before spilling to local memory.\n")
        print("     N   route C (~3N)   route A (N^2)   route B (2(N+2)^2)   fibers/elem")
        print("  " + "-" * 72)
        for r in rows:
            fa = "" if r["a"] <= 128 else "  SPILL"
            fb = "  SPILL"
            print(f"  {r['N']:4d}   {r['c']:11d}   {r['a']:13d}{fa:7s} {r['b']:13d}{fb:7s}"
                  f" {r['fib']:9d}")
        print("\n  Work per fiber (model flops):\n")
        print("     N   route C (Aberth)   dense QR (~10N^3)   ratio")
        print("  " + "-" * 52)
        for r in rows:
            print(f"  {r['N']:4d}   {r['flops_c']:16.0f}   {r['flops_qr']:17.0f}"
                  f"   {r['flops_qr']/r['flops_c']:5.2f}x")
        print("\n  Total fibers for this mesh (864 elements, 3(N+1)^2 per element):")
        for r in rows:
            print(f"    N={r['N']:2d}: {r['fib']*864:9d} independent rootfinding problems")
        print("\n  Branches in the hot loop: NONE (fixed count, no pivoting, no shifts,")
        print("  no deflation-during-iteration, no convergence test).")
    return rows


def detail(N, field, level=None, niter=NITER_DEFAULT):
    """Single-order dump: coefficients, iteration trace, spectrum, roots."""
    t = TESTS[field]
    lv = t.level if level is None else level
    x, w = gll(N)
    f = t.raw(x) - lv
    print(f"\n=== ROUTE C detail: N={N}, field='{field}' ({t.desc}), level={lv:g},"
          f" niter={niter} ===\n")
    print(f"  GLL nodes   x_j = {np.array2string(np.asarray(x), precision=6, max_line_width=100)}")
    print(f"  nodal f_j-c     = {np.array2string(f, precision=6, max_line_width=100)}")

    a_full = nodal_to_modal(N) @ f
    a = deflate_degree(a_full, A_FLOOR_DEFAULT)
    n = len(a) - 1
    print(f"\n  Legendre coefficients (deflated to degree {n} of {N}):")
    for k, ak in enumerate(a):
        bar = "#" * int(np.clip(np.log10(abs(ak) + 1e-300) + 17, 0, 40))
        print(f"    a_{k:<2d} = {ak:+.12e}   {bar}")
    print(f"\n  per-fiber state: a[{n+1}] real + z[{n}] complex = {(n+1)+2*n} doubles")
    print(f"  NO MATRIX IS FORMED (route A would build a {n}x{n} comrade matrix here)")

    z0 = aberth_init(n)
    print(f"\n  initial estimates (Chebyshev points, nudged off the real axis):")
    for zi in z0:
        print(f"    {zi.real:+.6f} {zi.imag:+.6f}i")

    z, hist = aberth_iterate(a, z0, niter, track=True)
    print(f"\n  iteration trace (max |correction|):")
    for i, h in enumerate(hist):
        bar = "#" * int(np.clip(-np.log10(max(h, 1e-18)) * 2, 0, 40))
        print(f"    it {i+1:2d}:  {h:.3e}   {bar}")

    print(f"\n  final estimates:")
    for zi in sorted(z, key=lambda v: (abs(v.imag) > 1e-7, v.real)):
        tag = "  <- real root in [-1,1]" if (abs(zi.imag) < 1e-7 and abs(zi.real) <= 1) else ""
        print(f"    {zi.real:+.12f} {zi.imag:+.12f}i{tag}")

    r, info = roots_route_c(N, f, niter=niter, return_info=True)
    ex = t.exact() if level is None else None
    print(f"\n  roots found ({len(r)}): {np.array2string(r, precision=12)}")
    print(f"  final correction = {info['final_corr']:.2e}"
          "   (safety net: large => this fiber did not settle)")
    if ex is not None:
        m = match_roots(r, ex)
        print(f"  exact roots  ({len(ex)}): {np.array2string(ex, precision=12)}")
        print(f"  max error = {m['max_err']:.3e}   missed = {m['n_missed']}   "
              f"spurious = {m['n_spurious']}")


def plot_iterations(hist_rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    cmap = plt.get_cmap("viridis")
    cols = [cmap(i / max(1, len(hist_rows) - 1) * 0.88) for i in range(len(hist_rows))]

    for (k, h), c in zip(hist_rows, cols):
        it = np.arange(1, len(h) + 1)
        ax[0].semilogy(it, np.maximum(h, 1e-18), "o-", color=c, lw=1.6, ms=4, label=k)
    ax[0].axhline(2.2e-16, color="k", ls=":", lw=1)
    ax[0].text(1, 3e-16, "machine eps", fontsize=8, va="bottom")
    ax[0].axvline(NITER_DEFAULT, color="#d62728", ls="--", lw=1.2)
    ax[0].text(NITER_DEFAULT + 0.2, 1e-4, f"default\nniter={NITER_DEFAULT}",
               fontsize=8, color="#d62728")
    ax[0].set_xlabel("iteration")
    ax[0].set_ylabel(r"max $|$correction$|$")
    ax[0].set_title("Ehrlich-Aberth convergence (N=8): cubic, then stalls at eps",
                    fontsize=10)
    ax[0].grid(alpha=0.3, which="both")
    ax[0].legend(fontsize=8)

    niters = [6, 8, 10, 12, 16, 24]
    mismatch = [57, 2, 0, 0, 0, 0]          # measured, scratch/23 (20,000 fibers)
    maxdiff = [5.7e-08, 4.0e-08, 2.0e-13, 2.0e-13, 2.0e-13, 2.0e-13]
    ax2 = ax[1]
    ax2.bar(range(len(niters)), mismatch, color="#d62728", alpha=0.75, width=0.55)
    ax2.set_xticks(range(len(niters)))
    ax2.set_xticklabels([str(v) for v in niters])
    ax2.set_xlabel("fixed iteration count")
    ax2.set_ylabel("root-count mismatches vs route B", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    ax2.set_title("20,000 random fibers: a fixed count is safe from niter=10",
                  fontsize=10)
    for i, v in enumerate(mismatch):
        ax2.text(i, v + 1.5, str(v), ha="center", fontsize=9, color="#d62728")
    ax3 = ax2.twinx()
    ax3.semilogy(range(len(niters)), maxdiff, "o-", color="#1f77b4", lw=1.8, ms=5)
    ax3.set_ylabel(r"max $|$route C $-$ route B$|$", color="#1f77b4")
    ax3.tick_params(axis="y", labelcolor="#1f77b4")
    ax2.set_ylim(0, 65)

    fig.suptitle("Route C: the fixed-iteration claim", fontsize=11, y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description="Route C demo: matrix-free Ehrlich-Aberth",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("--study", default="all",
                    choices=["all", "convergence", "iterations", "cost", "gpu",
                             "reference"])
    ap.add_argument("-N", type=int, default=None, help="single-order detail dump")
    ap.add_argument("--field", default="sine3", choices=sorted(TESTS.keys()))
    ap.add_argument("--level", type=float, default=None)
    ap.add_argument("--orders", default=None, help="comma list, e.g. 5,7,9")
    ap.add_argument("--niter", type=int, default=NITER_DEFAULT,
                    help=f"fixed iteration count (default {NITER_DEFAULT})")
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--figdir", default=os.path.join(HERE, "..", "figs"))
    args = ap.parse_args()

    if args.N is not None:
        detail(args.N, args.field, args.level, args.niter)
        return

    orders = ([int(s) for s in args.orders.split(",")] if args.orders
              else studies.ORDERS)
    route = lambda N, f: roots_route_c(N, f, niter=args.niter)
    os.makedirs(args.figdir, exist_ok=True)

    print("=" * 80)
    print(f"  {LABEL}   [niter = {args.niter}]")
    print("  single GLL reference element on [-1,1], orders " +
          ",".join(map(str, orders)))
    print("=" * 80)

    if args.study in ("all", "iterations"):
        hist_rows, _ = study_iterations()
        if not args.no_fig:
            p = plot_iterations(hist_rows, os.path.join(args.figdir, "C_iterations.png"))
            print(f"\n  figure -> {os.path.relpath(p, HERE)}")

    if args.study in ("all", "convergence"):
        res = studies.study_convergence(route, LABEL, orders)
        if not args.no_fig:
            p = studies.plot_convergence(res, LABEL,
                                         os.path.join(args.figdir, "C_convergence.png"),
                                         orders)
            print(f"\n  figure -> {os.path.relpath(p, HERE)}")

    if args.study in ("all", "cost"):
        rows = studies.study_cost(route, LABEL, orders,
                                  nnz_fn=lambda N: (state_doubles(N), N * N),
                                  flop_model=lambda N: flop_model_c(N, args.niter))
        print("\n  NO eigensolver and NO matrix: a simultaneous fixed-point iteration")
        print("  on n complex numbers.  The 'nnz/dense' column above is repurposed as")
        print("  (per-fiber state in doubles) / (route A's matrix entries) for contrast.")
        print("  Work is O(niter * n^2), against O(n^3) for any dense eigensolver.")
        if not args.no_fig:
            p = studies.plot_cost(rows, LABEL, os.path.join(args.figdir, "C_cost.png"),
                                  f"2*{args.niter}*N^2")
            print(f"  figure -> {os.path.relpath(p, HERE)}")

    if args.study in ("all", "gpu"):
        study_gpu()

    if args.study in ("all", "reference"):
        studies.study_reference(route, LABEL)


if __name__ == "__main__":
    main()
