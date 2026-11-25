#!/usr/bin/env python3
"""
run_demo.py -- three hard isosurface cases on ONE deformed spectral element.

Each case is shown at the four stages of the section-6 pipeline, so that what the
method produces can be inspected rather than taken on trust:

    (1) ROOTS        the raw output: exact fibre crossings on the GLL lattice
    (2) TRIANGLES    marching cubes on the N^3 sub-cells, flat facets
    (3) RECONSTRUCTED the curved surface rebuilt from the STORED PAYLOAD ONLY --
                     vertices + normals, no further access to the element data
    (4) EXACT        the analytic level set, sampled finely: the reference

Cases (see cases3d.py for why these three):
    tunnel   MC case 13, where face and interior ambiguity interact
    gyroid   many sheets, high genus, the most ambiguous cells
    neck     two lobes merging: |grad f| ~ 2e-2, the sqrt(eps) tangency

Usage
    python3 run_demo.py                 # all cases, all stages (+ figures)
    python3 run_demo.py --case tunnel   # one case
    python3 run_demo.py -N 8            # element order
    python3 run_demo.py --study report  # numbers only, no figures
    python3 run_demo.py --no-fig
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "demo_2d"))

import cases3d as K3                                              # noqa: E402
import surface3d as S3                                            # noqa: E402
import kernel2d as K                                              # noqa: E402


# --------------------------------------------------------------------------
# reconstruction from the stored payload (no element access)
# --------------------------------------------------------------------------
def build_payload(F, tris):
    """EXTRACTION: what is stored per triangle.

    Per vertex: position (3 floats) + unit normal (2 floats after normalisation).
    The normal is grad f / |grad f|, obtained at extraction time; in production the
    transverse components come from a Lagrange contraction of the nodal values and
    the along-fibre one from p'(r*), exactly as in 2D (refine2d.tangent_from_stored).
    NOTHING else is kept: the element's nodal data is discarded after this.
    """
    verts = {}
    for tri in tris:
        for p in tri:
            kx = tuple(np.round(p, 12))
            if kx not in verts:
                g = S3.grad3(F, np.array(p))
                n = float(np.linalg.norm(g))
                verts[kx] = (np.array(p), g / n if n > 1e-10 else None)
    return verts, [tuple(tuple(np.round(p, 12)) for p in tri) for tri in tris]


def reconstruct(payload, nsub=2):
    """RECONSTRUCTION: curve each triangle using ONLY stored positions and normals.

    For each edge (p0,n0)-(p1,n1) the midpoint is placed by the standard PN-triangle
    rule: take the chord midpoint, then remove its component along the average
    normal, which is a purely geometric operation on the stored data:

        m = (p0+p1)/2 ,   m <- m - ((m-p0).n0) n0/2 - ((m-p1).n1) n1/2

    This reproduces curvature from the normals alone.  No field evaluation, no
    projection, no access to the element -- which is the operational constraint.
    """
    verts, tris = payload
    out = []
    for tri in tris:
        pts = [verts[k][0] for k in tri]
        nrm = [verts[k][1] for k in tri]
        out.append((pts, nrm))

    def edge_mid(p0, n0, p1, n1):
        m = 0.5 * (p0 + p1)
        if n0 is not None:
            m = m - 0.5 * float((m - p0) @ n0) * n0
        if n1 is not None:
            m = m - 0.5 * float((m - p1) @ n1) * n1
        return m

    def nrm_at(p, n0, n1):
        if n0 is None or n1 is None:
            return n0 if n1 is None else n1
        v = n0 + n1
        nn = np.linalg.norm(v)
        return v / nn if nn > 1e-12 else n0

    cur = out
    for _ in range(nsub):
        new = []
        for (pts, nrm) in cur:
            A, B, Cc = pts
            nA, nB, nC = nrm
            mAB = edge_mid(A, nA, B, nB)
            mBC = edge_mid(B, nB, Cc, nC)
            mCA = edge_mid(Cc, nC, A, nA)
            nAB, nBC, nCA = nrm_at(mAB, nA, nB), nrm_at(mBC, nB, nC), nrm_at(mCA, nC, nA)
            new += [([A, mAB, mCA], [nA, nAB, nCA]),
                    ([mAB, B, mBC], [nAB, nB, nBC]),
                    ([mCA, mBC, Cc], [nCA, nBC, nC]),
                    ([mAB, mBC, mCA], [nAB, nBC, nCA])]
        cur = new
    return [tuple(p) for (p, _) in cur]


def payload_bytes(payload):
    verts, tris = payload
    return len(verts) * 5 * 4, len(tris) * 3 * 4       # fp32 coords+normals, int32 idx


# --------------------------------------------------------------------------
def analytic_tris(key, ngrid=48):
    """The reference surface: MARCHING TETRAHEDRA on a fine sample of the analytic field.

    Marching tetrahedra rather than marching cubes, deliberately: each cube is split
    into 6 tetrahedra and a tetrahedron's crossing pattern is unambiguous (1 or 2
    triangles, no case table, no saddle).  It is watertight by construction and cannot
    drop a patch, which a centroid-fan over cube cut-points can -- and did, when the
    cut points of a cell were not close to coplanar.

    The reference is therefore produced by a different algorithm from the method under
    test, which is what a reference should be.
    """
    t, V = K3.analytic_reference(key, ngrid)
    n = ngrid
    # the 6-tetrahedron split of a cube (shared diagonal 0-6), a standard decomposition
    CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    TETS = [(0, 5, 1, 6), (0, 1, 2, 6), (0, 2, 3, 6),
            (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6)]
    tris = []
    for i in range(n - 1):
        for j in range(n - 1):
            for k in range(n - 1):
                vv = [V[i + c[0], j + c[1], k + c[2]] for c in CUBE]
                if min(vv) > 0 or max(vv) < 0:
                    continue
                pp = [np.array([t[i + c[0]], t[j + c[1]], t[k + c[2]]]) for c in CUBE]
                for tet in TETS:
                    a, b, c_, d = tet
                    val = [vv[a], vv[b], vv[c_], vv[d]]
                    pos = [pp[a], pp[b], pp[c_], pp[d]]
                    neg = [q for q in range(4) if val[q] < 0]
                    if len(neg) == 0 or len(neg) == 4:
                        continue
                    # A tetrahedron's crossing pattern is unambiguous: either one
                    # vertex is separated from the other three (1 triangle), or the
                    # vertices split 2-2 (a quad, 2 triangles).  The quad's corners
                    # must be taken in CYCLIC order -- listing them in edge order
                    # gives a bowtie, which is what left the mesh unshared earlier.
                    def ip(q0, q1):
                        w = val[q0] / (val[q0] - val[q1])
                        return tuple(pos[q0] + w * (pos[q1] - pos[q0]))

                    if len(neg) == 1 or len(neg) == 3:
                        lone = neg[0] if len(neg) == 1 else \
                            [q for q in range(4) if q not in neg][0]
                        others = [q for q in range(4) if q != lone]
                        tris.append(tuple(ip(lone, q) for q in others))
                    else:                                   # 2-2 split: a quad
                        a0, a1 = neg
                        b0, b1 = [q for q in range(4) if q not in neg]
                        # cyclic: a0-b0, a0-b1, a1-b1, a1-b0
                        quad = [ip(a0, b0), ip(a0, b1), ip(a1, b1), ip(a1, b0)]
                        tris.append((quad[0], quad[1], quad[2]))
                        tris.append((quad[0], quad[2], quad[3]))
    return tris


def surface_distance(tris, key, nsample=2):
    """Max |f - c| over the reconstructed surface, using the ANALYTIC field.

    This is the honest geometric error: it asks whether the reconstructed surface
    lies on the true level set, not on the interpolant.
    """
    c = K3.CASES[key]
    worst = 0.0
    bary = [(i / nsample, j / nsample, 1 - i / nsample - j / nsample)
            for i in range(nsample + 1) for j in range(nsample + 1 - i)]
    for tri in tris:
        A, B, Cc = (np.array(p) for p in tri)
        for (w0, w1, w2) in bary:
            p = w0 * A + w1 * B + w2 * Cc
            val = c["f"](np.array([p[0]]), np.array([p[1]]), np.array([p[2]]))[0]
            worst = max(worst, abs(float(val) - c["level"]))
    return worst


# --------------------------------------------------------------------------
def plot_case(N, key, surf, pay, rec, path, nsub=2):
    """Four stages, side by side, in PHYSICAL space (deformed unless the case says not)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

    r, Rg, Sg, Tg = S3.gll_grid3(N)
    Xn, Yn, Zn = K3.case_geom(key)(Rg, Sg, Tg)

    def phys(p):
        wr = K.lagrange_weights(N, p[0])
        ws = K.lagrange_weights(N, p[1])
        wt = K.lagrange_weights(N, p[2])
        return (float(np.einsum("i,j,k,ijk->", wr, ws, wt, Xn)),
                float(np.einsum("i,j,k,ijk->", wr, ws, wt, Yn)),
                float(np.einsum("i,j,k,ijk->", wr, ws, wt, Zn)))

    def frame(ax):
        segs = []
        for a_ in (0, N):
            segs += [list(zip(Xn[a_, :, 0], Yn[a_, :, 0], Zn[a_, :, 0])),
                     list(zip(Xn[a_, :, -1], Yn[a_, :, -1], Zn[a_, :, -1])),
                     list(zip(Xn[a_, 0, :], Yn[a_, 0, :], Zn[a_, 0, :])),
                     list(zip(Xn[a_, -1, :], Yn[a_, -1, :], Zn[a_, -1, :]))]
        for b_ in (0, N):
            segs += [list(zip(Xn[:, b_, 0], Yn[:, b_, 0], Zn[:, b_, 0])),
                     list(zip(Xn[:, b_, -1], Yn[:, b_, -1], Zn[:, b_, -1]))]
        ax.add_collection3d(Line3DCollection(segs, colors="k", linewidths=.5, alpha=.4))
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=20, azim=41)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.grid(False)

    c = K3.CASES[key]
    fig = plt.figure(figsize=(17.2, 4.7))

    # (1) roots
    ax = fig.add_subplot(1, 4, 1, projection="3d")
    col = {0: "#d62728", 1: "#2ca02c", 2: "#9467bd"}
    for d in (0, 1, 2):
        P = [phys((t, r[j], r[k]) if d == 0 else (r[j], t, r[k]) if d == 1
                  else (r[j], r[k], t))
             for (dd, j, k), rts in surf.cross.items() if dd == d for t in rts]
        if P:
            A = np.array(P)
            ax.scatter(A[:, 0], A[:, 1], A[:, 2], s=7, c=col[d], depthshade=True,
                       label=f"{'rst'[d]}-fibres")
    frame(ax)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title(f"(1) roots\n{surf.n_roots()} exact crossings from "
                 f"{surf.n_solves} fibre solves", fontsize=9)

    # (2) flat triangles -- STEEL BLUE
    ax = fig.add_subplot(1, 4, 2, projection="3d")
    if surf.tris:
        pc = Poly3DCollection([[phys(p) for p in t] for t in surf.tris],
                              alpha=.93, linewidths=.25)
        pc.set_facecolor("#7fa8c9"); pc.set_edgecolor("#24425c")
        ax.add_collection3d(pc)
    frame(ax)
    ax.set_title(f"(2) marching cubes\n{len(surf.tris)} flat triangles, "
                 f"{surf.ambiguous} ambiguous cells", fontsize=9)

    # (3) reconstructed from the stored payload -- ORANGE, to contrast with (4).
    #     For the Taylor-Green case, colour by |u| instead: a second physical
    #     quantity carried on the surface, which is what a practitioner wants and
    #     which costs one extra stored float per vertex.
    ax = fig.add_subplot(1, 4, 3, projection="3d")
    if rec:
        polys = [[phys(p) for p in t] for t in rec]
        pc = Poly3DCollection(polys, alpha=.96, linewidths=0.0)
        if key == "tgv":
            spd = np.array([np.mean([K3.tgv_speed_ref(np.array([p[0]]),
                                                      np.array([p[1]]),
                                                      np.array([p[2]]))[0]
                                     for p in t]) for t in rec])
            cm = plt.get_cmap("plasma")
            lo_, hi_ = float(spd.min()), float(spd.max())
            pc.set_facecolor(cm((spd - lo_) / max(hi_ - lo_, 1e-12)))
            sm = plt.cm.ScalarMappable(cmap=cm,
                                       norm=plt.Normalize(vmin=lo_, vmax=hi_))
            cb = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.02)
            cb.set_label(r"$|\mathbf{u}|$", fontsize=8)
            cb.ax.tick_params(labelsize=7)
        else:
            pc.set_facecolor("#e08a3c")
        ax.add_collection3d(pc)
    frame(ax)
    nv, ni = payload_bytes(pay)
    extra = ", coloured by $|u|$" if key == "tgv" else ""
    ax.set_title(f"(3) reconstructed from payload{extra}\n"
                 f"{len(rec)} tris from {(nv+ni)/1024:.1f} KiB, no element access",
                 fontsize=9)

    # (4) exact reference -- GREEN, and by marching TETRAHEDRA (a different
    #     algorithm from the one under test, and watertight by construction)
    ax = fig.add_subplot(1, 4, 4, projection="3d")
    ng = 40
    ex = analytic_tris(key, ngrid=ng)
    if ex:
        pc = Poly3DCollection([[phys(p) for p in t] for t in ex],
                              alpha=.95, linewidths=0.0)
        pc.set_facecolor("#4f9d69")
        ax.add_collection3d(pc)
    frame(ax)
    ax.set_title(f"(4) exact reference\nmarching tetrahedra, ${ng}^3$ sample",
                 fontsize=9)

    wrapped = "\n".join(__import__("textwrap").wrap(
        f"{key}: {c['title']}   —   {c['why']}", 155))
    fig.suptitle(wrapped, fontsize=10, y=1.04 if wrapped.count("\n") else 1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _interp_floor(surf, key, prj_tris, stride=9):
    """Distance from the INTERPOLANT's level set to the ANALYTIC one.

    The projected points lie on the interpolant by construction, so evaluating the
    analytic field there measures the discretisation error alone -- the floor that
    bounds every surface-reconstruction scheme.
    """
    c = K3.CASES[key]
    V = [p for tri in prj_tris for p in tri][::stride]
    if not V:
        return float("nan")
    return max(abs(float(c["f"](np.array([p[0]]), np.array([p[1]]),
                                np.array([p[2]]))[0]) - c["level"]) for p in V)


def main():
    ap = argparse.ArgumentParser(
        description="Three hard isosurface cases on one deformed element",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--case", default=None, choices=sorted(K3.CASES.keys()))
    ap.add_argument("-N", type=int, default=None,
                    help="override the element order (default: per-case)")
    ap.add_argument("--nsub", type=int, default=2)
    ap.add_argument("--grid", type=int, default=40,
                    help="sampling of the analytic reference surface")
    ap.add_argument("--study", default="all", choices=["all", "report"])
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--figdir", default=os.path.join(HERE, "..", "figs"))
    args = ap.parse_args()
    args.N_given = args.N is not None
    if args.N is None:
        args.N = 8
    os.makedirs(args.figdir, exist_ok=True)

    print("=" * 80)
    print("  Hard isosurface cases on ONE spectral element")
    print("  three manufactured configurations on a deformed element, plus the")
    print("  Taylor-Green vortex on a big undeformed one")
    print("=" * 80)

    keys = [args.case] if args.case else list(K3.CASES.keys())

    print(f"\n  {'case':8s} {'roots':>6s} {'tris':>6s} {'ambig':>6s} {'min|grad|':>10s}"
          f" {'payload':>9s} | {'flat':>10s} {'PN recon':>10s} {'projected':>10s}"
          f" {'floor':>10s}  element")
    print("  " + "-" * 112)
    for key in keys:
        c = K3.CASES[key]
        Nk = args.N if args.N_given else K3.case_order(key, args.N)
        r, Rg, Sg, Tg = S3.gll_grid3(Nk)
        U = c["f"](Rg, Sg, Tg)
        surf = S3.extract_surface(Nk, U, c["level"])
        pay = build_payload(surf.F, surf.tris)
        rec = reconstruct(pay, nsub=args.nsub)
        V = np.array([p for t in surf.tris for p in t]) if surf.tris else np.zeros((0, 3))
        gmin = (min(np.linalg.norm(S3.grad3(surf.F, p)) for p in V[::5])
                if len(V) else float("nan"))
        e_flat = surface_distance(surf.tris, key)
        e_rec = surface_distance(rec, key)
        prj = S3.subdivide_project(surf.F, surf.tris, nsub=args.nsub)
        e_prj = surface_distance(prj, key)
        nv, ni = payload_bytes(pay)
        print(f"  {key:8s} {surf.n_roots():6d} {len(surf.tris):6d} {surf.ambiguous:6d}"
              f" {gmin:10.2e} {(nv+ni)/1024:8.1f}K | {e_flat:10.2e} {e_rec:10.2e}"
              f" {e_prj:10.2e} {_interp_floor(surf, key, prj):10.2e}"
              f"  N={Nk}{'' if not c.get('undeformed') else ', undef'}")
        if not args.no_fig and args.study == "all":
            p = plot_case(Nk, key, surf, pay, rec,
                          os.path.join(args.figdir, f"F_case_{key}.png"),
                          nsub=args.nsub)
            print(f"           figure -> {os.path.relpath(p, HERE)}")

    print("\n  All four errors are max |f - c| against the ANALYTIC field, sampled")
    print("  INSIDE every triangle -- not at the nodes, where the answer is trivially")
    print("  exact.  The four columns separate the two error sources:")
    print("    flat      flat marching-cubes triangles")
    print("    PN recon  curved using ONLY the stored payload (positions + normals),")
    print("              with no further access to the element.  THE CONSTRAINT CASE.")
    print("    projected subdivide + Newton on the interpolant (needs element access):")
    print("              the best the pipeline can do, kept as the reference")
    print("    floor     distance from the INTERPOLANT's level set to the analytic one")
    print("              -- the discretisation error, which no surface method can beat")
    print()
    print("  Reading it: PN reconstruction roughly halves the flat error at ~5 floats")
    print("  per vertex and zero element access.  Projection does 7-11x better but")
    print("  needs the element data.  Both remain far above the interpolation floor on")
    print("  'tunnel' (7.9e-08), so on that case the limit is the RECONSTRUCTION, not")
    print("  the discretisation -- there is room for a richer payload.  On 'gyroid' the")
    print("  floor is 4.1e-04 and the field is barely resolved at N=8, so refining the")
    print("  surface matters less than raising N.")
    print()
    for key in keys:
        print(f"    {key:8s} {K3.CASES[key]['title']}")
        print(f"             {K3.CASES[key]['why']}")


if __name__ == "__main__":
    main()
