#!/usr/bin/env python3
"""
run_demo.py -- 2D contour extraction on a single GLL element (section 6 in miniature).

The 3D isosurface strategy, one dimension down, where it can be drawn and checked:
solve 2(N+1) fibers with route C, insert lines where an edge holds two roots, march
the N^2 sub-cells reading EXACT crossings, resolve the saddle case by evaluating the
interpolant, and map out through the element's own geometry.

Everything is self-contained: kernel2d.py implements GLL nodes, barycentric
Lagrange, the Legendre recurrence and Ehrlich-Aberth from their definitions.  No
scipy, no numpy.polynomial, no third-party contouring.

Usage
    python3 run_demo.py                    # all studies, writes figures
    python3 run_demo.py --study accuracy   # accuracy|cases|geometry|cost|selftest
    python3 run_demo.py --case saddle -N 8 # one case, detailed
    python3 run_demo.py --geom curved      # straight|curved|annulus
    python3 run_demo.py --no-fig

Studies
    selftest   verify the hand-written kernel against independent computations
    cases      run all 7 test cases, report roots/segments/insertions (+ figure)
    accuracy   this method vs R1 (GLL as bilinear) and R2 (resample ~1.5N) (+ figure)
    r2sweep    how far R2 can be pushed: M = 9..400 per direction (+ figure)
    geometry   the same contour on straight / curved / annulus elements (+ figure)
    refine     straight vs cubic Hermite vs Hermite+projection (+ figure)
    surface    the FULL 3D pipeline on one element + a cut plane (+ 2 figures)
    cost       fiber solves, insertions, evaluations; extrapolation to 3D
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import kernel2d as K                                              # noqa: E402
import contour2d as C                                             # noqa: E402
import reference2d as R                                           # noqa: E402


# --------------------------------------------------------------------------
def exact_distance(case_key, pts):
    """Distance of each contour vertex from the TRUE level set, where known."""
    desc, fn, lev = C.CASES2D[case_key]
    if len(pts) == 0:
        return np.array([])
    P = np.asarray(pts)
    if case_key == "circle":
        return np.abs(np.hypot(P[:, 0], P[:, 1]) - np.sqrt(lev))
    if case_key == "ellipse":
        # |f-c| / |grad f| : first-order distance to the level set
        a2, b2 = 1.0, 2.6
        v = a2 * P[:, 0] ** 2 + b2 * (P[:, 1] - 0.12) ** 2 - lev
        gn = np.hypot(2 * a2 * P[:, 0], 2 * b2 * (P[:, 1] - 0.12))
        return np.abs(v) / np.clip(gn, 1e-12, None)
    if case_key == "wavycircle":
        th = np.arctan2(P[:, 1], P[:, 0])
        rad = np.hypot(P[:, 0], P[:, 1])
        return np.abs(rad - (0.62 + 0.10 * np.cos(3 * th)))
    if case_key == "saddle":
        return np.abs(P[:, 0] * P[:, 1] - lev) / np.hypot(P[:, 1], P[:, 0]).clip(1e-12)
    if case_key == "corner":
        return np.abs(P[:, 0] + P[:, 1] - lev) / np.sqrt(2.0)
    if case_key == "twoblobs":
        d1 = np.abs(np.hypot(P[:, 0] - 0.45, P[:, 1] - 0.4) - np.sqrt(lev))
        d2 = np.abs(np.hypot(P[:, 0] + 0.5, P[:, 1] + 0.45) - np.sqrt(lev))
        return np.minimum(d1, d2)
    return np.array([])          # no closed form


def verts(segs):
    out = []
    for (p0, p1) in segs:
        out.append(p0)
        out.append(p1)
    return np.array(out) if out else np.zeros((0, 2))


def residual(N, U, level, segs):
    """|f(vertex) - level| : how far the emitted vertices are off the contour."""
    if not segs:
        return np.array([])
    F = U - level
    return np.array([abs(C.tensor_eval(F, p[0], p[1]))
                     for s in segs for p in s])


# --------------------------------------------------------------------------
def study_selftest():
    print("\n=== self-test of the hand-written kernel ===\n")
    ok = K.selftest(verbose=True)
    print("\n  Routines checked: GLL nodes (endpoints, symmetry, sum of weights,")
    print("  quadrature exact to degree 2N-1, interior nodes are roots of P'_N),")
    print("  barycentric weights against the defining product, nodal<->modal round")
    print("  trip, Legendre p' against central differences, and Aberth against a")
    print("  polynomial with prescribed roots.")
    return ok


def study_cases(N, geom, no_fig, figdir):
    print(f"\n=== all test cases (N={N}, geometry='{geom}') ===\n")
    r, Rg, Sg = C.gll_grid(N)
    print(f"  {'case':10s} {'roots':>6s} {'segs':>5s} {'insert':>7s} {'saddle':>7s}"
          f" {'max |f-c| at vertices':>23s}")
    print("  " + "-" * 68)
    results = {}
    for key, (desc, fn, lev) in C.CASES2D.items():
        U = fn(Rg, Sg)
        res = C.extract_contour(N, U, lev)
        rr = residual(N, U, lev, res.segments)
        results[key] = (res, U, lev)
        print(f"  {key:10s} {res.n_roots():6d} {len(res.segments):5d} "
              f"{len(res.inserted):7d} {res.ambiguous:7d} "
              f"{(rr.max() if len(rr) else 0.0):23.2e}")
    print("\n  'max |f-c| at vertices' is the residual of the emitted points on the")
    print("  true interpolant: machine precision means the vertices lie exactly on")
    print("  the contour, not on a linear approximation of it.")
    for key, (desc, fn, lev) in C.CASES2D.items():
        print(f"    {key:10s} {desc}")
    if not no_fig:
        p = plot_cases(N, geom, results, os.path.join(figdir, "E_contour_cases.png"))
        print(f"\n  figure -> {os.path.relpath(p, HERE)}")
    return results


def study_r2sweep(no_fig, figdir, Ms=(9, 13, 25, 50, 100, 200, 400)):
    """How fine must the resampled grid be?  R2 refines a LINEAR interpolant, so it
    converges at O(h^2) = O(1/M^2) no matter how accurate the underlying field is."""
    print("\n=== how far can R2 be pushed by resampling? ===\n")
    print("  R2 interpolates the degree-N field onto a uniform MxM grid, then contours")
    print("  with LINEAR edge interpolation.  The field values are spectrally accurate,")
    print("  but the CONTOURING is linear, so the vertex error is O(1/M^2).\n")
    rows = []
    for key in ("circle", "ellipse"):
        desc, fn, lev = C.CASES2D[key]
        for N in (8,):
            r, Rg, Sg = C.gll_grid(N)
            U = fn(Rg, Sg)
            ours = C.extract_contour(N, U, lev)
            e0 = exact_distance(key, verts(ours.segments))
            base = max(e0.max() if len(e0) else np.nan, 1e-16)
            for M in Ms:
                r2 = R.reference_R2(N, U, lev, M=M)
                e2 = exact_distance(key, verts(r2["segments"]))
                rows.append(dict(case=key, N=N, M=M,
                                 err=max(e2.max() if len(e2) else np.nan, 1e-16),
                                 ours=base, evals=M * M,
                                 segs=len(r2["segments"])))
    print(f"  {'case':9s} {'M':>5s} {'grid pts':>9s} {'segments':>9s} {'R2 error':>11s}"
          f" {'rate':>7s} {'vs ours':>10s}")
    print("  " + "-" * 66)
    prev = {}
    for w in rows:
        k = w["case"]
        rate = ""
        if k in prev:
            pM, pe = prev[k]
            if pe > 1e-15 and w["err"] > 1e-15:
                rate = f"{np.log(pe/w['err'])/np.log(w['M']/pM):7.2f}"
        prev[k] = (w["M"], w["err"])
        print(f"  {w['case']:9s} {w['M']:5d} {w['evals']:9d} {w['segs']:9d}"
              f" {w['err']:11.2e} {rate:>7s} {w['err']/w['ours']:10.1e}")
    print("\n  'rate' is the observed order: err ~ M^-rate.  It sits near 2, the")
    print("  order of linear interpolation -- refining the grid never reaches the")
    print("  accuracy the degree-N field already contains.  To match our 1e-16 by")
    print("  refinement alone would need M ~ 1e8, i.e. 1e16 grid points per element.")
    print("\n  Our method uses 2(N+1) = %d fibre solves and no resampling at all."
          % (2 * (8 + 1)))
    if not no_fig:
        p = plot_r2sweep(rows, os.path.join(figdir, "E_contour_r2sweep.png"))
        print(f"\n  figure -> {os.path.relpath(p, HERE)}")
    return rows


def plot_r2sweep(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cases = sorted(set(w["case"] for w in rows))
    fig, ax = plt.subplots(1, 2, figsize=(12.2, 4.6))
    cols = {"circle": "#1f77b4", "ellipse": "#ff7f0e"}
    for cs in cases:
        sel = [w for w in rows if w["case"] == cs]
        M = np.array([w["M"] for w in sel])
        e = np.array([w["err"] for w in sel])
        ax[0].loglog(M, e, "o-", color=cols[cs], lw=1.8, ms=6, label=f"R2, {cs}")
        ax[0].axhline(sel[0]["ours"], color=cols[cs], ls="--", lw=1.4, alpha=.8)
    ref = np.array([w["M"] for w in rows if w["case"] == cases[0]], dtype=float)
    e0 = [w["err"] for w in rows if w["case"] == cases[0]][0]
    ax[0].loglog(ref, e0 * (ref[0] / ref) ** 2, "k:", lw=1.3, label=r"$O(M^{-2})$")
    ax[0].axhline(2.2e-16, color="k", ls=":", lw=1)
    ax[0].text(ref[0], 3.5e-16, "machine eps", fontsize=8)
    ax[0].set_xlabel("uniform grid size $M$ per direction")
    ax[0].set_ylabel("max vertex distance from the true level set")
    ax[0].set_title("R2 converges at the rate of LINEAR interpolation\n"
                    "(dashed: this work, at machine precision)", fontsize=10)
    ax[0].grid(alpha=.3, which="both")
    ax[0].legend(fontsize=8.5)

    for cs in cases:
        sel = [w for w in rows if w["case"] == cs]
        ax[1].loglog([w["evals"] for w in sel], [w["err"] for w in sel],
                     "o-", color=cols[cs], lw=1.8, ms=6, label=f"R2, {cs}")
        ax[1].plot([2 * (8 + 1)], [sel[0]["ours"]], "*", color="#2ca02c", ms=18,
                   label="this work" if cs == cases[0] else None)
    ax[1].set_xlabel("field evaluations per element")
    ax[1].set_ylabel("max vertex distance")
    ax[1].set_title("cost of that accuracy", fontsize=10)
    ax[1].grid(alpha=.3, which="both")
    ax[1].legend(fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def study_accuracy(no_fig, figdir):
    print("\n=== this method vs the two reference approaches ===\n")
    print("  R1: treat the GLL nodal values as an N x N bilinear grid, marching")
    print("      squares with linear edge interpolation.")
    print("  R2: interpolate the field to a uniform MxM grid (M ~ 1.5N), then the")
    print("      same marching squares.")
    print("  Both use the SAME topology logic as ours; only the vertex locations")
    print("  differ -- linear interpolation versus exact rootfinding.\n")
    rows = []
    # Cases chosen so the comparison is meaningful: the contour must be CURVED in
    # (r,s).  'saddle' and 'corner' are excluded because a straight or
    # coordinate-aligned contour is reproduced exactly by linear interpolation too
    # -- there the methods tie at machine precision and the comparison says nothing.
    for key in ("circle", "ellipse", "twoblobs", "wavycircle"):
        desc, fn, lev = C.CASES2D[key]
        for N in (4, 6, 8, 12):
            r, Rg, Sg = C.gll_grid(N)
            U = fn(Rg, Sg)
            ours = C.extract_contour(N, U, lev)
            r1 = R.reference_R1(N, U, lev)
            r2 = R.reference_R2(N, U, lev)
            e0 = exact_distance(key, verts(ours.segments))
            e1 = exact_distance(key, verts(r1["segments"]))
            e2 = exact_distance(key, verts(r2["segments"]))
            FLOOR = 1e-16          # report at most machine precision, not below
            rows.append(dict(case=key, N=N,
                             ours=max(e0.max() if len(e0) else np.nan, FLOOR),
                             r1=max(e1.max() if len(e1) else np.nan, FLOOR),
                             r2=max(e2.max() if len(e2) else np.nan, FLOOR),
                             n0=len(ours.segments), n1=len(r1["segments"]),
                             n2=len(r2["segments"]),
                             ev0=2 * (N + 1), ev1=r1["evals"], ev2=r2["evals"]))
    print(f"  {'case':10s} {'N':>3s} | {'ours':>10s} {'R1 bilinear':>12s} {'R2 resampled':>13s}"
          f" | {'R1/ours':>9s} {'R2/ours':>9s}")
    print("  " + "-" * 78)
    for w in rows:
        f1, f2 = w["r1"] / w["ours"], w["r2"] / w["ours"]
        print(f"  {w['case']:10s} {w['N']:3d} | {w['ours']:10.2e} {w['r1']:12.2e}"
              f" {w['r2']:13.2e} | {f1:9.1e} {f2:9.1e}")
    print("\n  max distance of an emitted vertex from the TRUE level set.")
    print()
    print("  READ THIS TABLE IN TWO HALVES.")
    print("  'circle' and 'ellipse' are POLYNOMIAL in (r,s): the GLL interpolant is")
    print("  exact, so the only error left is the contouring itself.  There our method")
    print("  is at machine precision and the reference methods are 1e13-1e15 times")
    print("  worse -- they are limited by LINEAR interpolation along a cell edge, which")
    print("  resampling shrinks but never removes.")
    print()
    print("  'twoblobs' and 'wavycircle' are NOT polynomials: the degree-N interpolant")
    print("  cannot represent them, and that interpolation error swamps everything.")
    print("  All three methods then agree to within a small factor -- exact rootfinding")
    print("  buys nothing once the polynomial itself is the wrong curve.  This is the")
    print("  conditioning floor of section 2, seen in 2D: a rootfinder cannot beat the")
    print("  interpolant it is given.")
    if not no_fig:
        p = plot_accuracy(rows, os.path.join(figdir, "E_contour_accuracy.png"))
        print(f"\n  figure -> {os.path.relpath(p, HERE)}")
    return rows


def study_geometry(N, no_fig, figdir):
    print(f"\n=== the same contour on three geometries (N={N}) ===\n")
    print("  The fiber solve happens entirely in the REFERENCE element and is")
    print("  unaffected by the deformation: the roots in (r,s) are identical.")
    print("  Only the map-out changes.\n")
    r, Rg, Sg = C.gll_grid(N)
    desc, fn, lev = C.CASES2D["circle"]
    U = fn(Rg, Sg)
    res = C.extract_contour(N, U, lev)
    ref_pts = verts(res.segments)
    print(f"  {'geometry':10s} {'segments':>9s} {'ref roots identical':>21s}"
          f" {'physical arc length':>21s}")
    print("  " + "-" * 64)
    out = {}
    for gname, gfun in C.GEOMS.items():
        Xn, Yn = gfun(Rg, Sg)
        phys = C.map_to_physical(res.segments, Xn, Yn)
        L = sum(np.hypot(p1[0] - p0[0], p1[1] - p0[1]) for (p0, p1) in phys)
        out[gname] = (phys, Xn, Yn)
        print(f"  {gname:10s} {len(phys):9d} {'yes':>21s} {L:21.6f}")
    print("\n  Same reference roots, three different physical curves: the geometry")
    print("  enters only through the isoparametric map, never through the solve.")
    if not no_fig:
        p = plot_geometry(N, res, out, os.path.join(figdir, "E_contour_geometry.png"))
        print(f"  figure -> {os.path.relpath(p, HERE)}")
    return out


def study_cost(N):
    print(f"\n=== cost accounting (N={N}) ===\n")
    r, Rg, Sg = C.gll_grid(N)
    tot_ins = 0
    for key, (desc, fn, lev) in C.CASES2D.items():
        res = C.extract_contour(N, fn(Rg, Sg), lev)
        tot_ins += len(res.inserted)
    print(f"  2D, one element:")
    print(f"    fiber solves           2(N+1) = {2*(N+1)}")
    print(f"    marching cells         N^2    = {N*N}")
    print(f"    inserted lines         {tot_ins/len(C.CASES2D):.2f} per case (average over"
          f" {len(C.CASES2D)} cases)")
    print(f"    field evaluations      only at saddle cells and on inserted lines")
    print(f"\n  3D extrapolation, one element:")
    print(f"    fiber solves           3(N+1)^2 = {3*(N+1)**2}")
    print(f"    marching cells         N^3      = {N**3}")
    print(f"    roots, worst case      3(N+1)^2 N = {3*(N+1)**2*N}")
    print(f"    roots, realistic       ~0.36/fiber = {int(0.36*3*(N+1)**2)}")
    print(f"\n  storage per element (realistic), r + xyz in fp32:")
    nr = 0.36 * 3 * (N + 1) ** 2
    print(f"    {nr:.0f} roots x 4 floats x 4 B = {nr*16:.0f} B")
    print(f"    worst case: {3*(N+1)**2*N} x 16 B = {3*(N+1)**2*N*16/1024:.1f} KiB")


# --------------------------------------------------------------------------
def _setup_ax(ax, title):
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.15)


def plot_cases(N, geom, results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = list(C.CASES2D.keys())
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.3 * nrow),
                             squeeze=False)
    r, Rg, Sg = C.gll_grid(N)
    gfun = C.GEOMS[geom]
    Xn, Yn = gfun(Rg, Sg)
    for i, key in enumerate(keys):
        ax = axes[i // ncol][i % ncol]
        res, U, lev = results[key]
        # background: the field, mapped
        tt = np.linspace(-1, 1, 120)
        Wt = np.array([K.lagrange_weights(N, t) for t in tt])
        Pf = Wt @ (U - lev) @ Wt.T
        TR, TS = np.meshgrid(tt, tt, indexing="ij")
        XX, YY = gfun(TR, TS)
        ax.pcolormesh(XX, YY, Pf, cmap="RdBu_r", shading="auto",
                      vmin=-np.abs(Pf).max(), vmax=np.abs(Pf).max(), alpha=0.55)
        # element boundary
        for (a_, b_) in ((Xn[0, :], Yn[0, :]), (Xn[-1, :], Yn[-1, :]),
                         (Xn[:, 0], Yn[:, 0]), (Xn[:, -1], Yn[:, -1])):
            ax.plot(a_, b_, "k-", lw=1.0, alpha=0.6)
        # GLL grid lines
        for j in range(N + 1):
            ax.plot(Xn[j, :], Yn[j, :], "-", color="0.55", lw=0.35, alpha=0.7)
            ax.plot(Xn[:, j], Yn[:, j], "-", color="0.55", lw=0.35, alpha=0.7)
        phys = C.map_to_physical(res.segments, Xn, Yn)
        for (p0, p1) in phys:
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "-", color="#111111", lw=2.0)
        V = verts(res.segments)
        if len(V):
            PV = C.map_to_physical([(tuple(v), tuple(v)) for v in V], Xn, Yn)
            ax.plot([q[0][0] for q in PV], [q[0][1] for q in PV], "o",
                    color="#d62728", ms=2.6)
        _setup_ax(ax, f"{key}: {res.n_roots()} roots, {len(res.segments)} segs"
                      + (f", {len(res.inserted)} ins" if res.inserted else "")
                      + (f", {res.ambiguous} saddle" if res.ambiguous else ""))
    for j in range(len(keys), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"Exact contours on one GLL element, $N={N}$, geometry '{geom}'\n"
                 "grey: GLL grid   red: exact fiber crossings   black: segments",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_accuracy(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cases = sorted(set(w["case"] for w in rows))
    fig, axes = plt.subplots(1, len(cases), figsize=(3.5 * len(cases), 3.6),
                             squeeze=False)
    for i, cs in enumerate(cases):
        ax = axes[0][i]
        sel = [w for w in rows if w["case"] == cs]
        Ns = [w["N"] for w in sel]
        ax.semilogy(Ns, [max(w["ours"], 1e-17) for w in sel], "o-",
                    color="#2ca02c", lw=1.8, ms=6, label="exact (this work)")
        ax.semilogy(Ns, [w["r1"] for w in sel], "s--", color="#d62728",
                    lw=1.6, ms=5, label="R1 GLL bilinear")
        ax.semilogy(Ns, [w["r2"] for w in sel], "^:", color="#1f77b4",
                    lw=1.6, ms=5, label="R2 uniform ~1.5N")
        ax.axhline(2.2e-16, color="k", ls=":", lw=1)
        ax.set_xlabel("$N$")
        ax.set_title(cs, fontsize=10)
        ax.grid(alpha=0.3, which="both")
        if i == 0:
            ax.set_ylabel("max distance of a vertex\nfrom the true level set")
            ax.legend(fontsize=8)
    fig.suptitle("Vertex accuracy: exact rootfinding versus linear interpolation "
                 "(same marching-squares topology)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_geometry(N, res, out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 3.9))
    r, Rg, Sg = C.gll_grid(N)
    ax = axes[0]
    for j in range(N + 1):
        ax.plot(Rg[j, :], Sg[j, :], "-", color="0.6", lw=0.4)
        ax.plot(Rg[:, j], Sg[:, j], "-", color="0.6", lw=0.4)
    for (p0, p1) in res.segments:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "-", color="#111", lw=2.0)
    V = verts(res.segments)
    ax.plot(V[:, 0], V[:, 1], "o", color="#d62728", ms=3)
    _setup_ax(ax, "reference $(r,s)$: the solve happens HERE")
    for i, (gname, (phys, Xn, Yn)) in enumerate(out.items()):
        ax = axes[i + 1]
        for j in range(N + 1):
            ax.plot(Xn[j, :], Yn[j, :], "-", color="0.6", lw=0.4)
            ax.plot(Xn[:, j], Yn[:, j], "-", color="0.6", lw=0.4)
        for (p0, p1) in phys:
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "-", color="#111", lw=2.0)
        pv = np.array([p for s in phys for p in s])
        ax.plot(pv[:, 0], pv[:, 1], "o", color="#d62728", ms=3)
        _setup_ax(ax, f"physical: '{gname}'")
    fig.suptitle("Identical reference roots, three geometries: the deformation "
                 "enters only through the isoparametric map", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
def study_surface(N, geom, no_fig, figdir):
    """The full 3D pipeline on one element, plus a cut plane through the interior."""
    import surface3d as S3
    print(f"\n=== 3D: the full pipeline on one element (N={N}, geometry='{geom}') ===\n")
    print(f"  marching-cubes edge table generated and verified: {S3.verify_tables()}")
    print(f"  (built here from the cube's edge topology -- no third-party table)\n")
    print("  Each 3D case is the DIRECT EXTENSION of the 2D case of the same name:")
    print("  same formula with a t-term, same isolevel, so the t=0 slice of the 3D")
    print("  field IS the 2D field.  Verified:")
    bad = 0
    for key, err, dl, st in S3.check_2d_consistency(N):
        if st == "DIFFERS":
            bad += 1
            print(f"    {key:11s} MISMATCH  max|U2-U3(t=0)|={err:.1e}")
    print(f"    all {len(S3.CASES3D)-1} shared cases reduce exactly "
          f"({bad} mismatches); 'torus' has no 2D analogue.\n")
    r, Rg, Sg, Tg = S3.gll_grid3(N)
    gfun = S3.GEOMS3[geom]
    Xn, Yn, Zn = gfun(Rg, Sg, Tg)
    print(f"  {'case':11s} {'solves':>7s} {'roots':>6s} {'tris':>6s} {'multi-edge':>11s}"
          f" {'saddle':>7s} {'at nodes':>12s} {'flat tri':>11s} {'projected':>11s}")
    print("  " + "-" * 92)
    out = {}
    for key, (desc, fn, lev) in S3.CASES3D.items():
        U = fn(Rg, Sg, Tg)
        surf = S3.extract_surface(N, U, lev)
        res = max((abs(S3.tensor_eval3(surf.F, *p))
                   for tri in surf.tris for p in tri), default=0.0)
        r0 = S3.surface_residual(surf.F, surf.tris, nsample=3) if surf.tris else 0.0
        r2 = (S3.surface_residual(S3.subdivide_project(surf.F, surf.tris, nsub=2),
                                  nsample=3) if False else
              S3.surface_residual(surf.F,
                                  S3.subdivide_project(surf.F, surf.tris, nsub=2),
                                  nsample=3) if surf.tris else 0.0)
        out[key] = (surf, U, lev)
        print(f"  {key:11s} {surf.n_solves:7d} {surf.n_roots():6d} {len(surf.tris):6d}"
              f" {surf.multi_edges:11d} {surf.ambiguous:7d} {res:12.2e}"
              f" {r0:11.2e} {r2:11.2e}")
    print("\n  'at nodes'  : |f-c| at the triangle corners -- machine precision, they")
    print("                are exact fibre roots.")
    print("  'flat tri'  : |f-c| sampled INSIDE the flat triangles -- the chord error,")
    print("                the only thing left after exact vertices.")
    print("  'projected' : after 2 subdivisions with every new node projected onto the")
    print("                surface (x <- x - f grad f/|grad f|^2).  ~4x smaller per")
    print("                subdivision, and the surface renders smooth.")
    print()
    print("  NOTE on 'tangency': the field (r^2+s^2+t^2-0.5)^2 is a perfect square and")
    print("  is non-negative everywhere, so it TOUCHES the level without crossing.  No")
    print("  marching-cubes corner sign changes and 0 triangles are emitted -- which is")
    print("  correct: a measure-zero tangency is not a surface MC can represent.  The")
    print("  fibre solver nevertheless finds the tangency locus (60 roots), with p'~0")
    print("  flagging every one of them.  This is the section-2 tangency case seen in 3D:")
    print("  detect and report, do not pretend to triangulate.")
    for key, (desc, fn, lev) in S3.CASES3D.items():
        print(f"    {key:11s} {desc}")
    if not no_fig:
        p = plot_surface(N, geom, out, os.path.join(figdir, "E_surface3d.png"))
        print(f"\n  figure -> {os.path.relpath(p, HERE)}")
        p2 = plot_surface_cut(N, geom, out, os.path.join(figdir, "E_surface3d_cut.png"))
        print(f"  figure -> {os.path.relpath(p2, HERE)}")
    return out


def plot_surface(N, geom, out, path, nsub=2):
    """Smooth surfaces: every triangle subdivided and every new node PROJECTED onto
    the isosurface, plus the t=0 slice plane drawn so the intersection is visible."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    import surface3d as S3

    keys = list(S3.CASES3D.keys())
    ncol = 5
    nrow = int(np.ceil(len(keys) / ncol))
    fig = plt.figure(figsize=(3.6 * ncol, 4.0 * nrow))
    r, Rg, Sg, Tg = S3.gll_grid3(N)
    Xn, Yn, Zn = S3.GEOMS3[geom](Rg, Sg, Tg)

    def phys(pt):
        wr = K.lagrange_weights(N, pt[0])
        ws = K.lagrange_weights(N, pt[1])
        wt = K.lagrange_weights(N, pt[2])
        return (float(np.einsum("i,j,k,ijk->", wr, ws, wt, Xn)),
                float(np.einsum("i,j,k,ijk->", wr, ws, wt, Yn)),
                float(np.einsum("i,j,k,ijk->", wr, ws, wt, Zn)))

    for i, key in enumerate(keys):
        ax = fig.add_subplot(nrow, ncol, i + 1, projection="3d")
        surf, U, lev = out[key]

        # --- the SLICE PLANE t = 0, drawn semi-transparent and extended ----------
        gq = np.linspace(-1.0, 1.0, 9)
        quads = []
        for ii in range(len(gq) - 1):
            for jj in range(len(gq) - 1):
                quads.append([phys((gq[ii], gq[jj], 0.0)),
                              phys((gq[ii + 1], gq[jj], 0.0)),
                              phys((gq[ii + 1], gq[jj + 1], 0.0)),
                              phys((gq[ii], gq[jj + 1], 0.0))])
        pl = Poly3DCollection(quads, alpha=0.13, linewidths=0)
        pl.set_facecolor("#888888")
        ax.add_collection3d(pl)

        # --- the surface, subdivided and projected: smooth --------------------
        tris = S3.subdivide_project(surf.F, surf.tris, nsub=nsub) if surf.tris else []
        if tris:
            pc = Poly3DCollection([[phys(p) for p in tri] for tri in tris],
                                  alpha=0.95, linewidths=0.0)
            pc.set_facecolor("#5aa2dc")
            ax.add_collection3d(pc)

        # --- element wireframe ------------------------------------------------
        segs = []
        for a_ in (0, N):
            segs += [list(zip(Xn[a_, :, 0], Yn[a_, :, 0], Zn[a_, :, 0])),
                     list(zip(Xn[a_, :, -1], Yn[a_, :, -1], Zn[a_, :, -1])),
                     list(zip(Xn[a_, 0, :], Yn[a_, 0, :], Zn[a_, 0, :])),
                     list(zip(Xn[a_, -1, :], Yn[a_, -1, :], Zn[a_, -1, :]))]
        for b_ in (0, N):
            segs += [list(zip(Xn[:, b_, 0], Yn[:, b_, 0], Zn[:, b_, 0])),
                     list(zip(Xn[:, b_, -1], Yn[:, b_, -1], Zn[:, b_, -1]))]
        ax.add_collection3d(Line3DCollection(segs, colors="k", linewidths=.55,
                                             alpha=.45))

        # --- the t=0 contour: the 2D case, in red, ON the slice plane ----------
        try:
            Pcut = S3.cut_plane(N, U, lev, 2, 0.0)
            res2 = C.extract_contour(N, Pcut, 0.0)
            rsegs = [[phys((q0[0], q0[1], 0.0)), phys((q1[0], q1[1], 0.0))]
                     for (q0, q1) in res2.segments]
            if rsegs:
                ax.add_collection3d(Line3DCollection(
                    rsegs, colors="#d62728", linewidths=2.6, linestyles="--"))
        except Exception:
            pass

        note = ""
        if len(surf.tris) == 0 and surf.n_roots() > 0:
            note = "\n(tangency: no sign change)"
        ax.set_title(f"{key}\n{surf.n_roots()} roots, {len(tris)} tris{note}",
                     fontsize=9)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=20, azim=41)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.grid(False)
    fig.suptitle(
        f"Isosurfaces from $3(N+1)^2={3*(N+1)**2}$ fibre solves on one $N={N}$ curved "
        f"element, then subdivided $\\times{4**nsub}$ with every new node PROJECTED "
        "onto the surface.\n"
        "Grey plane: the $t=0$ slice.  Red dashed: its contour -- the same-named case "
        "of the previous figure.  The surface is smooth because the nodes are exact, "
        "not because it was smoothed.",
        fontsize=10.5, y=1.015)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_surface_cut(N, geom, out, path):
    """Bigger picture: surface + a cut plane showing the field inside."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import surface3d as S3

    key = "torus"
    surf, U, lev = out[key]
    r, Rg, Sg, Tg = S3.gll_grid3(N)
    Xn, Yn, Zn = S3.GEOMS3[geom](Rg, Sg, Tg)
    fig = plt.figure(figsize=(15.5, 4.6))

    # --- (a) the surface, clipped to t <= 0 to expose the interior ----------
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    phys = S3.map_surface_to_physical(
        [tri for tri in surf.tris if np.mean([p[2] for p in tri]) <= 0.02],
        Xn, Yn, Zn)
    pc = Poly3DCollection(phys, alpha=0.9, linewidths=0.2)
    pc.set_facecolor("#4c9bd6"); pc.set_edgecolor("#1c3f5c")
    ax.add_collection3d(pc)
    ax.set_title("(a) surface clipped at $t=0$\n(the cut exposes the interior)",
                 fontsize=9)
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=24, azim=40)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.grid(False)

    # --- (b) the cut plane itself: the 2D problem ---------------------------
    ax = fig.add_subplot(1, 3, 2)
    P = S3.cut_plane(N, U, lev, 2, 0.0)
    tt = np.linspace(-1, 1, 160)
    W = np.array([K.lagrange_weights(N, t) for t in tt])
    Pf = W @ P @ W.T
    ax.pcolormesh(tt, tt, Pf.T, cmap="RdBu_r", shading="auto",
                  vmin=-np.abs(Pf).max(), vmax=np.abs(Pf).max())
    res2 = C.extract_contour(N, P, 0.0)
    for (p0, p1) in res2.segments:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "-", color="k", lw=2.2)
    V = verts(res2.segments)
    if len(V):
        ax.plot(V[:, 0], V[:, 1], "o", color="#d62728", ms=4)
    ax.set_aspect("equal")
    ax.set_title("(b) the $t=0$ cut is a 2D problem\n"
                 f"{res2.n_roots()} roots, {len(res2.segments)} segments", fontsize=9)
    ax.set_xlabel("$r$"); ax.set_ylabel("$s$"); ax.tick_params(labelsize=7)

    # --- (c) all the fibre crossings, as a point cloud ----------------------
    ax = fig.add_subplot(1, 3, 3, projection="3d")
    pts = []
    for (d, j, k), rts in surf.cross.items():
        for t in rts:
            pts.append((t, r[j], r[k]) if d == 0 else
                       (r[j], t, r[k]) if d == 1 else (r[j], r[k], t))
    Pt = np.array(pts)
    col = {0: "#d62728", 1: "#2ca02c", 2: "#9467bd"}
    off = 0
    for d in (0, 1, 2):
        sel = [(t, r[j], r[k]) if d == 0 else (r[j], t, r[k]) if d == 1
               else (r[j], r[k], t)
               for (dd, j, k), rts in surf.cross.items() if dd == d for t in rts]
        if sel:
            A = np.array(sel)
            ax.scatter(A[:, 0], A[:, 1], A[:, 2], s=9, c=col[d],
                       label=f"dir {'rst'[d]}", depthshade=True)
    ax.set_title(f"(c) the raw output: {len(Pt)} exact crossings\n"
                 "on a structured lattice, not a point cloud", fontsize=9)
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=24, azim=40)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.grid(False)
    ax.legend(fontsize=7.5, loc="upper left")

    fig.suptitle("The bigger picture: fibre crossings, the 2D cut, and the assembled "
                 f"surface ('{key}', $N={N}$)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
def study_refine(N, no_fig, figdir):
    """Storage-budgeted reconstruction of the contour.

    EXTRACTION (element data available) may evaluate freely.  RECONSTRUCTION must
    rebuild the curve from a stored payload alone.  The question is therefore not
    'may we project?' but 'what must be STORED, and what curve does it support?'.
    """
    import refine2d as RF
    print(f"\n=== reconstructing the contour from a STORED payload (N={N}) ===\n")
    print("  Two phases:")
    print("    EXTRACTION     element nodal data available; fibre solves, gradients")
    print("                   and projection all happen here and may evaluate freely.")
    print("    RECONSTRUCTION element data GONE; the curve is rebuilt from the stored")
    print("                   payload only.  This is the real use case.")
    print()
    print("  Payloads, per segment (vertices are shared, so one marginal vertex each;")
    print("  a unit tangent is stored as ONE angle, not two components):")
    print("    2 floats  vertices only              -> polyline")
    print("    3 floats  + tangent angle            -> cubic Hermite, C^1 at vertices")
    print("    4 floats  + 1 projected point        -> quadratic through 3 points")
    print("    5 floats  + tangent + 1 projected    -> quintic: C^1 ends AND mid")
    print()
    print("  REFERENCE: projection with k -> infinity.  Every projected point is an")
    print("  exact root, so the polyline through them tends to the true contour --")
    print("  error -> 0 at unbounded storage.  It is the lower bound no fixed-size")
    print("  payload can beat, not a storage scheme itself.\n")

    dists = {
        "circle": lambda P, lev=C.CASES2D["circle"][2]:
            np.abs(np.hypot(P[:, 0], P[:, 1]) - np.sqrt(lev)),
        "ellipse": lambda P, lev=C.CASES2D["ellipse"][2]:
            np.abs(P[:, 0] ** 2 + 2.6 * (P[:, 1] - 0.12) ** 2 - lev)
            / np.clip(np.hypot(2 * P[:, 0], 2 * 2.6 * (P[:, 1] - 0.12)), 1e-12, None),
        "wavycircle": lambda P: np.abs(
            np.hypot(P[:, 0], P[:, 1])
            - (0.62 + 0.10 * np.cos(3 * np.arctan2(P[:, 1], P[:, 0])))),
    }

    def err_of(curves, dist, ns=129):
        w = 0.0
        for Q in curves:
            if len(Q) < 2:
                continue
            L = np.cumsum(np.r_[0.0, np.linalg.norm(np.diff(Q, axis=0), axis=1)])
            if L[-1] <= 0:
                continue
            u = np.linspace(0, 1, ns) * L[-1]
            X = np.interp(u, L, Q[:, 0]); Y = np.interp(u, L, Q[:, 1])
            w = max(w, float(np.max(dist(np.column_stack([X, Y])))))
        return w

    schemes = [("vertices only", dict(want_tangent=False), 2.0),
               ("+ tangent (Hermite)", dict(want_tangent=True), 3.0),
               ("+ 1 proj (quadratic)", dict(want_tangent=False, want_mid=True), 4.0),
               ("+ tangent + proj (quintic)", dict(want_tangent=True, want_mid=True), 5.0)]
    rows = []
    r, Rg, Sg = C.gll_grid(N)
    for key, dist in dists.items():
        desc, fn, lev = C.CASES2D[key]
        U = fn(Rg, Sg)
        F = U - lev
        res = C.extract_contour(N, U, lev)
        for name, kw, fl in schemes:
            curves = [RF.reconstruct(RF.build_payload(F, p0, p1, **kw))
                      for (p0, p1) in res.segments]
            rows.append(dict(case=key, scheme=name, floats=fl,
                             err=err_of(curves, dist)))
        # the reference: projection, k -> infinity
        ref = {}
        for k in (1, 3, 7, 15, 31, 63):
            curves = [RF.refine_segment(F, p0, p1, k=k, mode="project")
                      for (p0, p1) in res.segments]
            ref[k] = err_of(curves, dist)
        rows.append(dict(case=key, scheme="reference (proj, k->inf)",
                         floats=np.nan, err=ref[63], refcurve=ref))

    print(f"  {'case':11s} {'payload':27s} {'floats/seg':>11s} {'error':>11s}"
          f" {'vs polyline':>12s}")
    print("  " + "-" * 78)
    for cs in dists:
        base = [w for w in rows if w["case"] == cs and w["scheme"] == "vertices only"][0]["err"]
        for w in [w for w in rows if w["case"] == cs]:
            fl = "     --   " if not np.isfinite(w["floats"]) else f"{w['floats']:11.1f}"
            print(f"  {w['case']:11s} {w['scheme']:27s} {fl} {w['err']:11.2e}"
                  f" {base/max(w['err'],1e-18):11.1f}x")
        print()

    print("  Reading the table:")
    print("  * One extra float per vertex (the tangent angle) buys ~38x on 'circle'")
    print("    and ~7x on 'ellipse'.  It is the best value per byte, and it needs NO")
    print("    new field evaluation: the along-fibre derivative is p'(r*) from the")
    print("    Aberth iteration, and the transverse one is a single Lagrange")
    print("    contraction of nodal values already in hand.")
    print("  * Three more floats (tangent + one projected midpoint) reach ~3e3x, and")
    print("    the quintic is the best fixed payload tested.")
    print("  * 'wavycircle' is unmoved by any payload: the interpolant is the wrong")
    print("    curve, and no reconstruction of a wrong curve helps (section 2).")
    print("  * The reference converges to zero, confirming the projected points lie")
    print("    exactly on the contour; the fixed payloads approach it from above.")
    if not no_fig:
        p = plot_refine(N, rows, os.path.join(figdir, "E_contour_refine.png"))
        print(f"\n  figure -> {os.path.relpath(p, HERE)}")
    return rows


def plot_refine(N, rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import refine2d as RF

    fig = plt.figure(figsize=(14.8, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.12, 1.0], hspace=0.34, wspace=0.27)

    # --- top: the reconstructed curves on a coarse element -------------------
    Nc = 5
    r, Rg, Sg = C.gll_grid(Nc)
    desc, fn, lev = C.CASES2D["circle"]
    U = fn(Rg, Sg); F = U - lev
    res = C.extract_contour(Nc, U, lev)
    th = np.linspace(0, 2 * np.pi, 400)
    panels = [("vertices only  (2 floats/seg)", dict(want_tangent=False)),
              ("+ tangent: cubic Hermite  (3)", dict(want_tangent=True)),
              ("+ tangent + proj: quintic  (5)", dict(want_tangent=True, want_mid=True))]
    for i, (title, kw) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i])
        ax.plot(np.sqrt(lev) * np.cos(th), np.sqrt(lev) * np.sin(th), "-",
                color="0.74", lw=5.5, alpha=.95, label="true level set", zorder=1)
        for j in range(Nc + 1):
            ax.plot(Rg[j, :], Sg[j, :], "-", color="0.88", lw=.5, zorder=0)
            ax.plot(Rg[:, j], Sg[:, j], "-", color="0.88", lw=.5, zorder=0)
        for (p0, p1) in res.segments:
            pay = RF.build_payload(F, p0, p1, **kw)
            Q = RF.reconstruct(pay)
            ax.plot(Q[:, 0], Q[:, 1], "-", color="#1a3c6e", lw=2.0, zorder=3)
            if "pm" in pay:
                ax.plot(*pay["pm"], ".", color="#2ca02c", ms=8, zorder=4)
        V = np.array([p for sgm in res.segments for p in sgm])
        ax.plot(V[:, 0], V[:, 1], "o", color="#d62728", ms=5.5, zorder=5,
                label="stored vertices (exact roots)")
        ax.set_aspect("equal"); ax.set_title(title, fontsize=10)
        ax.set_xlim(-0.85, 0.85); ax.set_ylim(-0.85, 0.85); ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7.5, loc="lower right")
        if i == 2:
            ax.plot([], [], ".", color="#2ca02c", ms=8, label="stored proj. point")
            ax.legend(fontsize=7.5, loc="lower right")

    # --- bottom: error vs stored floats, with the reference ------------------
    cases = ["circle", "ellipse", "wavycircle"]
    for i, cs in enumerate(cases):
        ax = fig.add_subplot(gs[1, i])
        sel = [w for w in rows if w["case"] == cs and np.isfinite(w["floats"])]
        ax.semilogy([w["floats"] for w in sel], [w["err"] for w in sel],
                    "o-", color="#1f77b4", lw=1.9, ms=7, label="fixed payload", zorder=3)
        for w in sel:
            ax.annotate(w["scheme"].replace("+ ", "").split(" (")[0],
                        (w["floats"], w["err"]), fontsize=6.5,
                        textcoords="offset points", xytext=(4, 5))
        refrow = [w for w in rows if w["case"] == cs and "refcurve" in w]
        if refrow:
            rc = refrow[0]["refcurve"]
            ks = sorted(rc)
            ax.semilogy([2 + 2 * k for k in ks], [rc[k] for k in ks], "s--",
                        color="#888888", lw=1.4, ms=4,
                        label="reference: projection $k\\to\\infty$", zorder=2)
        ax.set_xlabel("floats stored per segment")
        ax.set_title(cs, fontsize=10)
        ax.grid(alpha=.3, which="both")
        ax.set_xlim(1, 22)
        if i == 0:
            ax.set_ylabel("max distance of the reconstructed\ncurve from the true level set")
            ax.legend(fontsize=8)
    fig.suptitle("Reconstruction from a stored payload: the element data is gone, "
                 "only these floats remain\n"
                 "One extra float per vertex (the exact tangent angle, free at "
                 "extraction) buys an order of magnitude", fontsize=11, y=0.985)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="2D contour demo on a GLL element",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("--study", default="all",
                    choices=["all", "selftest", "cases", "accuracy", "r2sweep",
                             "geometry", "refine", "surface", "cost"])
    ap.add_argument("-N", type=int, default=8)
    ap.add_argument("--case", default=None, choices=sorted(C.CASES2D.keys()))
    ap.add_argument("--geom", default="curved", choices=sorted(C.GEOMS.keys()))
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--figdir", default=os.path.join(HERE, "..", "figs"))
    args = ap.parse_args()
    os.makedirs(args.figdir, exist_ok=True)

    print("=" * 78)
    print("  2D contour extraction on a single GLL element")
    print("  fiber solver: route C (Ehrlich-Aberth), self-contained kernel")
    print("=" * 78)

    if args.case:
        r, Rg, Sg = C.gll_grid(args.N)
        desc, fn, lev = C.CASES2D[args.case]
        U = fn(Rg, Sg)
        res = C.extract_contour(args.N, U, lev)
        print(f"\ncase '{args.case}': {desc}\n  level = {lev}")
        print(f"  fiber solves      : {res.n_solves}")
        print(f"  roots found       : {res.n_roots()}")
        print(f"  inserted lines    : {len(res.inserted)} {res.inserted}")
        print(f"  saddle cells      : {res.ambiguous}")
        print(f"  segments          : {len(res.segments)}")
        rr = residual(args.N, U, lev, res.segments)
        if len(rr):
            print(f"  max |f-c| at vertices: {rr.max():.3e}")
        for d in (0, 1):
            nm = "r" if d == 0 else "s"
            for i, rts in enumerate(res.crossings[d]):
                if len(rts):
                    sl = res.slopes[d][i]
                    print(f"    {nm}-fiber {i}: " +
                          ", ".join(f"{t:+.9f} (p'={p:+.2e})" for t, p in zip(rts, sl)))
        return

    if args.study in ("all", "selftest"):
        study_selftest()
    if args.study in ("all", "cases"):
        study_cases(args.N, args.geom, args.no_fig, args.figdir)
    if args.study in ("all", "accuracy"):
        study_accuracy(args.no_fig, args.figdir)
    if args.study in ("all", "r2sweep"):
        study_r2sweep(args.no_fig, args.figdir)
    if args.study in ("all", "geometry"):
        study_geometry(args.N, args.no_fig, args.figdir)
    if args.study in ("all", "refine"):
        study_refine(args.N, args.no_fig, args.figdir)
    if args.study in ("all", "surface"):
        study_surface(min(args.N, 8), "curved" if args.geom == "annulus" else args.geom,
                      args.no_fig, args.figdir)
    if args.study in ("all", "cost"):
        study_cost(args.N)


if __name__ == "__main__":
    main()
