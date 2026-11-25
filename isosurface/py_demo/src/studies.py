"""
studies.py -- the study/benchmark/plot engine shared by demo_A and demo_B.

Each demo passes in its own route callable and a label; everything below is
route-agnostic so the two demos are strictly comparable (same fields, same
metric, same timing protocol, same figure style).

Studies
    convergence   error vs N for every test field        -> table + figure
    cost          us/solve vs N, with nnz and flop model -> table + figure
    degeneracy    behaviour as the top mode a_N -> 0     -> table
    identity      det(A - lam B) == +/- p(lam)           -> table  (route B only)
    reference     contrast vs monomial companion + Brent -> table
"""
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gllbase import (gll, TESTS, match_roots, timeit_us,          # noqa: E402
                     ref_roots_companion, ref_roots_brent, nodal_to_modal,
                     interp_error, root_error_floor)

ORDERS = list(range(3, 16))
FIELD_ORDER = ["poly5", "tangent", "gauss", "tanh", "runge", "sine2", "sine3", "sine5"]


# --------------------------------------------------------------------------
def _fmt(e):
    if e is None:
        return "  --  "
    if not np.isfinite(e):
        return " FAIL "
    if e == 0.0:
        return " 0.0  "
    return f"{e:.1e}"


def study_convergence(route, label, orders=ORDERS, keys=None, verbose=True):
    """max root error vs N for each test field.  Returns dict[key] -> array."""
    keys = keys or FIELD_ORDER
    res = {}
    for k in keys:
        t = TESTS[k]
        errs, miss, spur, floor = [], [], [], []
        for N in orders:
            x, _ = gll(N)
            try:
                r = route(N, t.values(x))
                m = match_roots(r, t.exact())
                errs.append(m["max_err"])
                miss.append(m["n_missed"])
                spur.append(m["n_spurious"])
            except Exception:
                errs.append(np.inf)
                miss.append(-1)
                spur.append(-1)
            floor.append(root_error_floor(N, t))
        res[k] = dict(err=np.array(errs), missed=np.array(miss),
                      spurious=np.array(spur), floor=np.array(floor),
                      nroots=len(t.exact()))
    if verbose:
        print(f"\n=== {label}: convergence -- max |root error| vs N ===\n")
        hdr = "  field      nr | " + " ".join(f"N={N:<5d}" for N in orders)
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for k in keys:
            row = " ".join(f"{_fmt(e):7s}" for e in res[k]["err"])
            print(f"  {k:<10s} {res[k]['nroots']:2d} | {row}")
        print("\n  root-displacement floor  |delta-p(x*)| / |p'(x*)|  -- the CONDITIONING of")
        print("  the problem: how far the true root moves because the interpolant is not the")
        print("  field.  numerator = interpolation discrepancy AT the root; denominator =")
        print("  crossing steepness (1/|p'| is the root condition number).  A route sitting")
        print("  ON this floor is exact -- the remaining error is interpolation, not solving.")
        print("  (inf = double root, where p'(x*)=0 and the linear theory fails; see")
        print("  gllbase.root_error_floor.)\n")
        for k in keys:
            row = " ".join(f"{_fmt(e):7s}" for e in res[k]["floor"])
            print(f"  {k:<10s} {res[k]['nroots']:2d} | {row}")
        print("\n  (missed / spurious roots, blank = clean)")
        for k in keys:
            ms = res[k]["missed"]
            sp = res[k]["spurious"]
            if np.any(ms != 0) or np.any(sp != 0):
                cells = " ".join(f"{m}/{s}    " if (m or s) else "       "
                                 for m, s in zip(ms, sp))
                print(f"  {k:<10s}    | {cells}")
    return res


# --------------------------------------------------------------------------
def study_cost(route, label, orders=ORDERS, nnz_fn=None, flop_model=None,
               verbose=True):
    """Wall time per solve vs N (median of 7 batches), plus nnz and a flop model."""
    rng = np.random.default_rng(12345)
    rows = []
    for N in orders:
        f = rng.standard_normal(N + 1)
        us = timeit_us(lambda N=N, f=f: route(N, f))
        nnz = nnz_fn(N) if nnz_fn else None
        fl = flop_model(N) if flop_model else None
        rows.append(dict(N=N, us=us, nnz=nnz, flops=fl))
    if verbose:
        print(f"\n=== {label}: cost per solve (1 core, median of 7 batches) ===\n")
        print("  WARNING: this measures PYTHON dispatch, not the algorithm.  A route built")
        print("  from many small numpy calls (route C: ~600-1600 per solve) is penalised")
        print("  100x against one that makes a single LAPACK call (routes A and B).  Check")
        print("  the implied flop rate: anything under ~100 Mflop/s on a modern core is")
        print("  overhead.  For a real cost see demo_C/aberth_bench.c and notes/12.\n")
        print("    N    us/solve   solves/s    nnz / dense      model flops")
        print("  " + "-" * 58)
        for r in rows:
            nz = (f"{r['nnz'][0]:4d} /{r['nnz'][1]:5d} "
                  f"({100*r['nnz'][0]/r['nnz'][1]:4.1f}%)") if r["nnz"] else " " * 20
            fl = f"{r['flops']:9.0f}" if r["flops"] else ""
            print(f"  {r['N']:3d}   {r['us']:8.1f}   {1e6/r['us']:8.0f}   {nz}  {fl}")
    return rows


# --------------------------------------------------------------------------
def study_degeneracy(route, label, N=7, eps_list=(1e-2, 1e-4, 1e-6, 1e-8,
                                                  1e-10, 1e-12, 1e-14, 0.0),
                     verbose=True):
    """Drive the top Legendre mode a_N to zero and watch what the route does.

    A well-resolved SEM element HAS a tiny a_N -- that is what spectral
    convergence means -- so this is the normal regime, not a corner case.
    """
    from numpy.polynomial import legendre as L
    x, _ = gll(N)
    V = np.array([[L.legval(xj, [0] * k + [1]) for k in range(N + 1)] for xj in x])
    base = np.array([0.3, -0.2, 0.5, 0.1, -0.4, 0.25, 0.15, 0.0])[:N + 1]

    # reference: the exact roots of the degree-(N-1) polynomial (a_N = 0 exactly)
    a0 = base.copy()
    a0[N] = 0.0
    ref = None
    try:
        ref = np.sort(ref_roots_companion(N, V @ a0))
    except Exception:
        pass

    rows = []
    for eps in eps_list:
        a = base.copy()
        a[N] = eps
        f = V @ a
        try:
            r = route(N, f)
            err = (np.abs(r - ref).max() if (ref is not None and len(r) == len(ref))
                   else None)
            rows.append(dict(eps=eps, ok=True, roots=r, err=err))
        except Exception as e:
            rows.append(dict(eps=eps, ok=False, roots=None, err=None,
                             msg=str(e)[:60]))
    if verbose:
        print(f"\n=== {label}: degeneracy -- top mode a_N -> 0  (N={N}) ===\n")
        print(f"  reference roots (a_N = 0 exactly): {np.round(ref, 12) if ref is not None else '--'}\n")
        print("      a_N      status   #roots   roots in [-1,1]                  err vs ref")
        print("  " + "-" * 76)
        for r in rows:
            if not r["ok"]:
                print(f"  {r['eps']:9.1e}   FAILED   --       {r['msg']}")
            else:
                rs = ", ".join(f"{v:+.9f}" for v in r["roots"])
                print(f"  {r['eps']:9.1e}   ok       {len(r['roots']):2d}       "
                      f"{rs:<32s} {_fmt(r['err'])}")
    return rows


# --------------------------------------------------------------------------
def study_reference(route, label, orders=(5, 8, 11), keys=None, verbose=True):
    """Contrast against two standard alternatives on the same data."""
    keys = keys or FIELD_ORDER
    rows = []
    for k in keys:
        t = TESTS[k]
        for N in orders:
            x, _ = gll(N)
            f = t.values(x)
            out = {"field": k, "N": N}
            for name, fn in (("route", route),
                             ("companion", ref_roots_companion),
                             ("brent", ref_roots_brent)):
                try:
                    r = fn(N, f)
                    m = match_roots(r, t.exact())
                    out[name] = m["max_err"]
                    out[name + "_miss"] = m["n_missed"]
                    out[name + "_spur"] = m["n_spurious"]
                except Exception:
                    out[name] = np.inf
                    out[name + "_miss"] = -1
                    out[name + "_spur"] = -1
            rows.append(out)
    if verbose:
        print(f"\n=== {label}: vs reference methods -- max |root error| (miss/spurious) ===\n")
        print("  field       N  |  this route          monomial companion   dense-sample+Brent")
        print("  " + "-" * 78)
        for r in rows:
            def cell(n):
                return f"{_fmt(r[n]):7s} {r[n+'_miss']}/{r[n+'_spur']}"
            print(f"  {r['field']:<10s} {r['N']:2d}  |  {cell('route'):<19s} "
                  f"{cell('companion'):<20s} {cell('brent')}")
        print("\n  'brent' = sample the interpolant on 2001 points, bracket sign changes.")
        print("  That is the subgrid logic marching cubes uses; its misses are the")
        print("  motivation for this whole project.")
    return rows


# --------------------------------------------------------------------------
def plot_convergence(res, label, path, orders=ORDERS):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = [k for k in FIELD_ORDER if k in res]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    cmap = plt.get_cmap("viridis")
    cols = [cmap(i / max(1, len(keys) - 1) * 0.88) for i in range(len(keys))]

    for k, c in zip(keys, cols):
        e = np.array(res[k]["err"], dtype=float)
        e = np.where(np.isfinite(e), e, np.nan)
        e = np.maximum(e, 1e-17)                       # log axis floor
        ax[0].semilogy(orders, e, "o-", color=c, lw=1.6, ms=4,
                       label=f"{k} ({res[k]['nroots']}r)")
        fl = np.maximum(np.array(res[k].get("floor", []), dtype=float), 1e-17)
        if len(fl) == len(orders):
            ax[0].semilogy(orders, fl, "-", color=c, lw=0.9, alpha=0.35)
    ax[0].axhline(2.2e-16, color="k", ls=":", lw=1)
    ax[0].text(orders[0], 2.6e-16, "machine eps", fontsize=8, va="bottom")
    ax[0].set_xlabel("polynomial order $N$")
    ax[0].set_ylabel(r"max $|x_{\rm found}-x_{\rm exact}|$")
    ax[0].set_title(f"{label}: root error vs order\n"
                    r"(faint = conditioning floor $|\delta p(x^*)|/|p'(x^*)|$)",
                    fontsize=10)
    ax[0].grid(alpha=0.3, which="both")
    ax[0].legend(fontsize=8, ncol=2)
    ax[0].set_ylim(1e-17, 1e1)

    for k, c in zip(keys, cols):
        m = np.array(res[k]["missed"], dtype=float)
        s = np.array(res[k]["spurious"], dtype=float)
        ax[1].plot(orders, m, "o-", color=c, lw=1.6, ms=4, label=f"{k} missed")
        ax[1].plot(orders, s, "x--", color=c, lw=1.0, ms=5, alpha=0.65)
    ax[1].set_xlabel("polynomial order $N$")
    ax[1].set_ylabel("count")
    ax[1].set_title("missed (o, solid) / spurious (x, dashed) roots")
    ax[1].grid(alpha=0.3)
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].legend(fontsize=8, ncol=2)

    fig.suptitle("single GLL reference element, all real roots of $f(x)-c$ on $[-1,1]$",
                 fontsize=9, y=1.005)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_cost(rows, label, path, flop_label="model flops"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    N = np.array([r["N"] for r in rows])
    us = np.array([r["us"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.4))

    ax[0].plot(N, us, "o-", color="#1f77b4", lw=1.8, ms=5)
    ax[0].set_xlabel("polynomial order $N$")
    ax[0].set_ylabel(r"$\mu$s per solve")
    ax[0].set_title(f"{label}: measured cost")
    ax[0].grid(alpha=0.3)
    ax[0].set_ylim(bottom=0)
    ax[0].axvspan(5, 10, color="orange", alpha=0.12)
    ax[0].text(7.5, ax[0].get_ylim()[1] * 0.06, "typical SEM $N$",
               ha="center", fontsize=8, color="darkorange")

    if rows[0].get("flops"):
        fl = np.array([r["flops"] for r in rows], dtype=float)
        ax[1].loglog(N, fl, "s-", color="#d62728", lw=1.6, ms=5, label=flop_label)
        ref = fl[0] * (N / N[0]) ** 3
        ax[1].loglog(N, ref, "k:", lw=1.2, label=r"$O(N^3)$")
        ax[1].set_xlabel("polynomial order $N$")
        ax[1].set_ylabel("flops")
        ax[1].set_title("flop model (asymptotic, not measured)")
        ax[1].grid(alpha=0.3, which="both")
        ax[1].legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
def plot_testfields(path, N_show=8, keys=None):
    """BASELINE FIGURE: every test field, its isolevel, its exact roots, and the
    GLL interpolant at one representative order.

    This is the 'what are we actually solving' figure -- read it before any
    convergence plot.  Each panel shows
        solid grey   the true field f(x)
        dashed       the isolevel c
        coloured     the order-N GLL interpolant of f - c, shifted back by +c
        markers o    the exact roots of f - c
        markers x    the GLL nodes (on the isolevel line)
    so a panel where the coloured curve departs from the grey one is a field the
    element cannot resolve at that order -- and that, not the eigensolver, is what
    limits the root accuracy.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from numpy.polynomial import legendre as L

    keys = keys or FIELD_ORDER
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.9 * nrow),
                             squeeze=False)
    xs = np.linspace(-1, 1, 1200)

    for i, k in enumerate(keys):
        ax = axes[i // ncol][i % ncol]
        t = TESTS[k]
        x, _ = gll(N_show)
        a = nodal_to_modal(N_show) @ t.values(x)

        ax.plot(xs, t.raw(xs), "-", color="0.55", lw=2.2, label="$f(x)$")
        ax.plot(xs, L.legval(xs, a) + t.level, "-", color="#1f77b4", lw=1.3,
                label=f"GLL interp, $N$={N_show}")
        ax.axhline(t.level, color="#d62728", ls="--", lw=1.0,
                   label=f"level $c$={t.level:g}")
        ex = t.exact()
        ax.plot(ex, np.full_like(ex, t.level), "o", color="#d62728", ms=6,
                mfc="none", mew=1.6, label=f"{len(ex)} exact roots")
        ax.plot(x, np.full(len(x), t.level), "x", color="k", ms=4, alpha=0.5,
                label="GLL nodes")

        ie = interp_error(N_show, t)
        ax.set_title(f"{k}   ({len(ex)} roots,  interp err {ie:.1e})", fontsize=9)
        ax.set_xlim(-1.02, 1.02)
        ax.grid(alpha=0.25)
        ax.tick_params(labelsize=8)
        if i == 0:
            ax.legend(fontsize=6.5, loc="best")

    for j in range(len(keys), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    fig.suptitle("Baseline: the test fields on the single GLL reference element "
                 f"$[-1,1]$  (interpolant shown at $N$={N_show})",
                 fontsize=11, y=1.00)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_testfield_resolution(path, keys=None, orders=(5, 8, 11, 15)):
    """BASELINE FIGURE 2: how the interpolant approaches each field as N grows.

    Left  -- max |interpolant - field| vs N (the spectral-convergence picture).
    Right -- the same for the ROOT DISPLACEMENT FLOOR |delta-p(x*)|/|p'(x*)|, which
             is the quantity a rootfinder is actually judged against: the
             interpolation discrepancy at the root, divided by the crossing slope.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = keys or FIELD_ORDER
    orders_full = list(range(3, 16))
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    cmap = plt.get_cmap("viridis")
    cols = [cmap(i / max(1, len(keys) - 1) * 0.88) for i in range(len(keys))]

    for k, c in zip(keys, cols):
        t = TESTS[k]
        ie = [max(interp_error(N, t), 1e-17) for N in orders_full]
        fl = [max(root_error_floor(N, t), 1e-17) for N in orders_full]
        ax[0].semilogy(orders_full, ie, "o-", color=c, lw=1.5, ms=4, label=k)
        ax[1].semilogy(orders_full, fl, "o-", color=c, lw=1.5, ms=4, label=k)

    for a_, ttl, yl in ((ax[0], "interpolation error  $\\max|p_N - f|$",
                         r"$\max|p_N-f|$"),
                        (ax[1], "root displacement floor  $|\\delta p(x^*)|/|p'(x^*)|$",
                         r"floor on $|x_{\rm found}-x_{\rm exact}|$")):
        a_.axhline(2.2e-16, color="k", ls=":", lw=1)
        a_.set_xlabel("polynomial order $N$")
        a_.set_ylabel(yl)
        a_.set_title(ttl, fontsize=10)
        a_.grid(alpha=0.3, which="both")
        a_.legend(fontsize=8, ncol=2)
        a_.axvspan(5, 10, color="orange", alpha=0.10)

    fig.suptitle("Baseline: what the GLL interpolant can deliver, before any rootfinding "
                 "(shaded = typical SEM $N$)", fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
def plot_route_comparison(path, orders=(3, 5, 7, 9, 11, 13, 15)):
    """PAPER FIGURE: all three routes on the same axes, plus the conditioning floor.

    Left  -- max root error vs N for routes A/B/C on a representative field set;
             the three curves should lie on top of each other AND on the floor.
    Right -- the degeneracy sweep that separates them: plain A vs A+deflation vs
             B vs C as the top Legendre mode is driven to zero.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import sys as _sys, os as _os
    base = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
    for d in ("demo_A", "demo_B", "demo_C"):
        _sys.path.insert(0, _os.path.join(base, d))
    from route_a import roots_route_a
    from route_b import roots_route_b
    from route_c import roots_route_c
    from numpy.polynomial import legendre as L

    fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.7))

    keys = ["gauss", "sine2", "sine3"]
    styles = {"A": ("o-", "#1f77b4"), "B": ("s--", "#2ca02c"), "C": ("^:", "#d62728")}
    routes = {"A": roots_route_a, "B": roots_route_b, "C": roots_route_c}
    for ki, k in enumerate(keys):
        t = TESTS[k]
        alpha = 1.0 - 0.28 * ki
        for name, fn in routes.items():
            errs = []
            for N in orders:
                x, _ = gll(N)
                try:
                    errs.append(max(match_roots(fn(N, t.values(x)),
                                                t.exact())["max_err"], 1e-17))
                except Exception:
                    errs.append(np.nan)
            m, c = styles[name]
            ax[0].semilogy(orders, errs, m, color=c, lw=1.5, ms=5, alpha=alpha,
                           label=f"{k}: route {name}" if ki == 0 else None)
        fl = [max(root_error_floor(N, t), 1e-17) for N in orders]
        ax[0].semilogy(orders, fl, "-", color="0.35", lw=2.4, alpha=0.30,
                       label="conditioning floor" if ki == 0 else None)
    ax[0].axhline(2.2e-16, color="k", ls=":", lw=1)
    ax[0].set_xlabel("polynomial order $N$")
    ax[0].set_ylabel(r"max $|x_{\rm found}-x_{\rm exact}|$")
    ax[0].set_title("All three routes sit on the conditioning floor\n"
                    "(gauss, sine2, sine3; fading = successive fields)", fontsize=10)
    ax[0].grid(alpha=0.3, which="both")
    ax[0].legend(fontsize=7.5, ncol=2, loc="lower left")

    # --- degeneracy sweep -------------------------------------------------
    N = 7
    V = np.array([[L.legval(xj, [0] * k + [1]) for k in range(N + 1)]
                  for xj in gll(N)[0]])
    base_a = np.array([0.3, -0.2, 0.5, 0.1, -0.4, 0.25, 0.15, 0.0])
    eps = np.array([1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16])
    ref = None
    curves = {"A (plain)": [], "A + deflation": [], "B": [], "C": []}
    for e in eps:
        a = base_a.copy()
        a[N] = e
        f = V @ a
        r_b = roots_route_b(N, f)
        if ref is None:
            a0 = base_a.copy()
            a0[N] = 0.0
            ref = np.sort(roots_route_b(N, V @ a0))
        def err(r):
            r = np.sort(r)
            return (np.abs(r - ref).max() if len(r) == len(ref) else np.nan)
        try:
            curves["A (plain)"].append(err(roots_route_a(N, f, a_floor=0.0)))
        except Exception:
            curves["A (plain)"].append(np.nan)
        curves["A + deflation"].append(err(roots_route_a(N, f)))
        curves["B"].append(err(r_b))
        curves["C"].append(err(roots_route_c(N, f)))
    cstyle = {"A (plain)": ("o-", "#9467bd"), "A + deflation": ("o-", "#1f77b4"),
              "B": ("s--", "#2ca02c"), "C": ("^:", "#d62728")}
    for name, vals in curves.items():
        m, c = cstyle[name]
        ax[1].loglog(eps, np.maximum(np.array(vals, dtype=float), 1e-17), m,
                     color=c, lw=1.7, ms=5, label=name)
    ax[1].set_xlabel(r"leading Legendre coefficient $a_N$")
    ax[1].set_ylabel(r"error vs the $a_N\!=\!0$ limit")
    ax[1].set_title("Degeneracy sweep ($N=7$): only plain route A degrades",
                    fontsize=10)
    ax[1].grid(alpha=0.3, which="both")
    ax[1].legend(fontsize=8.5)
    ax[1].invert_xaxis()

    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_cost_reality(path):
    """PAPER FIGURE: the three cost regimes -- python scalar, python batched, compiled C.

    Numbers are measured; the compiled column comes from demo_C/aberth_bench.c.
    This figure exists to stop the python microbenchmark being read as a cost model.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    N = np.array([5, 8, 11, 15])
    py_scalar = np.array([484.6, 672.5, 912.5, 1179.7])     # scratch/26
    py_batch = np.array([9.16, 16.28, 30.58, 58.81])        # scratch/26
    compiled = np.array([3.10, 7.71, 14.80, 27.63])         # demo_C/aberth_bench.c
    route_a = np.array([27.5, 29.5, 34.8, 44.8])            # scratch/26

    fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.5))

    w = 0.26
    idx = np.arange(len(N))
    ax[0].bar(idx - w, py_scalar, w, label="route C, python scalar", color="#d62728", alpha=.85)
    ax[0].bar(idx, py_batch, w, label="route C, python batched", color="#ff7f0e", alpha=.85)
    ax[0].bar(idx + w, compiled, w, label="route C, compiled C", color="#2ca02c", alpha=.85)
    ax[0].plot(idx, route_a, "ko--", lw=1.4, ms=5, label="route A (one LAPACK call)")
    ax[0].set_yscale("log")
    ax[0].set_xticks(idx)
    ax[0].set_xticklabels([f"N={v}" for v in N])
    ax[0].set_ylabel(r"$\mu$s per fiber")
    ax[0].set_title("The same algorithm, three implementations\n"
                    "(the python microbenchmark measures dispatch, not arithmetic)",
                    fontsize=10)
    ax[0].grid(alpha=0.3, axis="y", which="both")
    ax[0].legend(fontsize=8)
    for i, (a_, b_) in enumerate(zip(py_scalar, compiled)):
        ax[0].text(i, a_ * 1.35, f"{a_/b_:.0f}$\\times$", ha="center", fontsize=8,
                   color="#d62728")

    real_flops = np.array([14640, 38112, 72600, 135720])    # scratch/26, counted
    rate_compiled = real_flops / compiled / 1e3
    rate_pyscalar = real_flops / py_scalar / 1e3
    rate_pybatch = real_flops / py_batch / 1e3
    ax[1].semilogy(N, rate_compiled, "o-", color="#2ca02c", lw=1.8, ms=6, label="compiled C")
    ax[1].semilogy(N, rate_pybatch, "s-", color="#ff7f0e", lw=1.6, ms=5, label="python batched")
    ax[1].semilogy(N, rate_pyscalar, "^-", color="#d62728", lw=1.6, ms=5, label="python scalar")
    ax[1].axhspan(1.0, 10.0, color="green", alpha=0.08)
    ax[1].text(5.2, 3.0, "realistic single-core range", fontsize=8, color="darkgreen")
    ax[1].set_xlabel("polynomial order $N$")
    ax[1].set_ylabel("achieved Gflop/s (real flops)")
    ax[1].set_title("Implied flop rate exposes the artifact", fontsize=10)
    ax[1].grid(alpha=0.3, which="both")
    ax[1].legend(fontsize=8.5)

    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path
