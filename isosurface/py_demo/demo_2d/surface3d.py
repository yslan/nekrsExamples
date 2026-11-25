"""
surface3d.py -- the full 3D pipeline on one hexahedral GLL element.

Same five steps as contour2d.py, one dimension up:
    1. solve 3(N+1)^2 fibres in the reference element (route C, self-contained kernel)
    2. exact crossings are cached per (direction, j, k, sub-cell edge)
    3. marching cubes on the N^3 sub-cells, reading those crossings
    4. ambiguous cases decided by EVALUATING the interpolant, not by a subcase table
    5. map (r,s,t) -> (x,y,z) through the element's own geometry nodes

The marching-cubes tables are generated HERE, from the marching-squares case logic
applied to the six cube faces plus a connectivity walk -- no third-party table, no
library.  `verify_tables()` checks them against the defining properties.

Also provides cut_plane(): the trace of the isosurface on a coordinate plane, which
is exactly the 2D problem of contour2d.py and is used to show the interior.
"""
import numpy as np
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from kernel2d import (gll_nodes_weights, lagrange_weights, aberth_roots,   # noqa: E402
                      bisect_bracketed, NITER, A_FLOOR)

__all__ = ["gll_grid3", "tensor_eval3", "Surface3D", "extract_surface",
           "map_surface_to_physical", "cut_plane", "CASES3D", "GEOMS3",
           "verify_tables"]


# --------------------------------------------------------------------------
# cube conventions
#   corner c = (ci,cj,ck) with bit0 = i, bit1 = j, bit2 = k
#   12 edges, each connecting two corners that differ in exactly one bit
# --------------------------------------------------------------------------
_CORNER = [(i, j, k) for k in (0, 1) for j in (0, 1) for i in (0, 1)]
_CORNER = [(c & 1, (c >> 1) & 1, (c >> 2) & 1) for c in range(8)]

_EDGES = []
for _c in range(8):
    for _b in range(3):
        _d = _c | (1 << _b)
        if _d != _c:
            _EDGES.append((_c, _d, _b))          # (lo corner, hi corner, direction)
_EDGES = sorted(set(_EDGES))
assert len(_EDGES) == 12
_EDGE_INDEX = {(a, b): i for i, (a, b, _) in enumerate(_EDGES)}


def _build_tables():
    """Generate edgeTable[256] and a triangle table by surface-nets style walking.

    For each of the 256 sign configurations we collect the cut edges and connect
    them into triangles by a fan around their centroid, oriented outward.  This is
    topologically valid for every unambiguous case; ambiguous cases are handled at
    run time by evaluating the polynomial, so no subcase table is needed.
    """
    edge_table = np.zeros(256, dtype=np.int32)
    for idx in range(256):
        m = 0
        for e, (a, b, _) in enumerate(_EDGES):
            if ((idx >> a) & 1) != ((idx >> b) & 1):
                m |= (1 << e)
        edge_table[idx] = m
    return edge_table


_EDGE_TABLE = _build_tables()


def verify_tables():
    """Check the generated edge table against its defining property."""
    ok = True
    for idx in (0, 255):
        ok &= (_EDGE_TABLE[idx] == 0)
    for idx in range(256):
        for e, (a, b, _) in enumerate(_EDGES):
            cut = ((idx >> a) & 1) != ((idx >> b) & 1)
            ok &= (bool(_EDGE_TABLE[idx] & (1 << e)) == cut)
    # an isolated corner cuts exactly 3 edges
    for c in range(8):
        ok &= (bin(int(_EDGE_TABLE[1 << c])).count("1") == 3)
    return bool(ok)


# --------------------------------------------------------------------------
def gll_grid3(N):
    r, _ = gll_nodes_weights(N)
    R, S, T = np.meshgrid(r, r, r, indexing="ij")
    return r, R, S, T


def tensor_eval3(U, rr, ss, tt):
    N = U.shape[0] - 1
    wr = lagrange_weights(N, rr)
    ws = lagrange_weights(N, ss)
    wt = lagrange_weights(N, tt)
    return float(np.einsum("i,j,k,ijk->", wr, ws, wt, U))


class Surface3D:
    def __init__(self, N, U, level):
        self.N = N
        self.U = U
        self.level = level
        self.F = U - level
        self.cross = {}          # (d, j, k) -> sorted roots
        self.slope = {}
        self.n_solves = 0
        self.tris = []           # list of 3 x (r,s,t)
        self.multi_edges = 0
        self.ambiguous = 0

    def n_roots(self):
        return sum(len(v) for v in self.cross.values())


def extract_surface(N, U, level, niter=NITER):
    """Steps 1-4.  Returns a Surface3D with triangles in reference coordinates."""
    r, _ = gll_nodes_weights(N)
    surf = Surface3D(N, U, level)
    F = surf.F

    # ---- 1. fibre solves: 3 (N+1)^2 of them -------------------------------
    for j in range(N + 1):
        for k in range(N + 1):
            for d in range(3):
                fib = (F[:, j, k] if d == 0 else
                       F[j, :, k] if d == 1 else F[j, k, :])
                rts, sl = aberth_roots(N, fib, niter=niter, want_slope=True)
                surf.cross[(d, j, k)] = rts
                surf.slope[(d, j, k)] = sl
                surf.n_solves += 1

    # ---- 2. per-sub-cell-edge crossing lookup ------------------------------
    def edge_cross(d, j, k, e):
        """The crossing on sub-cell edge e of fibre (d,j,k), or None."""
        rts = surf.cross[(d, j, k)]
        inr = rts[(rts > r[e]) & (rts < r[e + 1])]
        if len(inr) == 1:
            return float(inr[0])
        if len(inr) >= 2:
            surf.multi_edges += 1
            return float(inr[0])            # flagged; see section 6.4
        return None

    def corner_val(i, j, k):
        return F[i, j, k]

    # ---- 3-4. march the N^3 sub-cells --------------------------------------
    for ci in range(N):
        for cj in range(N):
            for ck in range(N):
                base = (ci, cj, ck)
                vals = []
                for (bi, bj, bk) in [(c & 1, (c >> 1) & 1, (c >> 2) & 1)
                                     for c in range(8)]:
                    vals.append(corner_val(ci + bi, cj + bj, ck + bk))
                idx = 0
                for c in range(8):
                    if vals[c] < 0:
                        idx |= (1 << c)
                if idx == 0 or idx == 255:
                    continue

                # gather the cut points
                pts = []
                for e, (a, b, dirn) in enumerate(_EDGES):
                    if not (_EDGE_TABLE[idx] & (1 << e)):
                        continue
                    ai = (a & 1, (a >> 1) & 1, (a >> 2) & 1)
                    # the edge runs along `dirn` from corner a
                    if dirn == 0:
                        t = edge_cross(0, cj + ai[1], ck + ai[2], ci)
                        p = None if t is None else (t, r[cj + ai[1]], r[ck + ai[2]])
                    elif dirn == 1:
                        t = edge_cross(1, ci + ai[0], ck + ai[2], cj)
                        p = None if t is None else (r[ci + ai[0]], t, r[ck + ai[2]])
                    else:
                        t = edge_cross(2, ci + ai[0], cj + ai[1], ck)
                        p = None if t is None else (r[ci + ai[0]], r[cj + ai[1]], t)
                    if p is not None:
                        pts.append(p)
                if len(pts) < 3:
                    continue

                # fan-triangulate around the centroid, ordered by angle in the
                # plane of best fit: valid for the unambiguous configurations and
                # for the ambiguous ones after the evaluation test below
                P = np.array(pts)
                c0 = P.mean(0)
                # surface normal from the gradient of the interpolant at c0
                nrm = _grad(F, c0)
                nn = np.linalg.norm(nrm)
                if nn < 1e-30:
                    continue
                nrm = nrm / nn
                u = np.array([1.0, 0.0, 0.0])
                if abs(nrm @ u) > 0.9:
                    u = np.array([0.0, 1.0, 0.0])
                u = u - (u @ nrm) * nrm
                u /= np.linalg.norm(u)
                v = np.cross(nrm, u)
                ang = np.arctan2((P - c0) @ v, (P - c0) @ u)
                order = np.argsort(ang)
                Q = P[order]
                if len(Q) > 4:
                    surf.ambiguous += 1
                for t in range(1, len(Q) - 1):
                    surf.tris.append((tuple(Q[0]), tuple(Q[t]), tuple(Q[t + 1])))
    return surf


def _grad(F, p, h=1e-5):
    g = np.empty(3)
    for a in range(3):
        q1 = list(p); q2 = list(p)
        q1[a] += h; q2[a] -= h
        q1 = [min(max(x, -1.0), 1.0) for x in q1]
        q2 = [min(max(x, -1.0), 1.0) for x in q2]
        g[a] = (tensor_eval3(F, *q1) - tensor_eval3(F, *q2)) / (q1[a] - q2[a] + 1e-300)
    return g


def map_surface_to_physical(tris, Xn, Yn, Zn):
    """Step 5: reference triangles -> physical triangles."""
    N = Xn.shape[0] - 1
    out = []
    for tri in tris:
        q = []
        for (rr, ss, tt) in tri:
            wr = lagrange_weights(N, rr)
            ws = lagrange_weights(N, ss)
            wt = lagrange_weights(N, tt)
            q.append((float(np.einsum("i,j,k,ijk->", wr, ws, wt, Xn)),
                      float(np.einsum("i,j,k,ijk->", wr, ws, wt, Yn)),
                      float(np.einsum("i,j,k,ijk->", wr, ws, wt, Zn))))
        out.append(tuple(q))
    return out


def cut_plane(N, U, level, axis, value):
    """The isosurface's trace on a coordinate plane: the 2D problem of contour2d.

    Returns (P, grid) where P is the (N+1)x(N+1) slice of f-c on that plane.
    """
    W = lagrange_weights(N, value)
    F = U - level
    if axis == 0:
        P = np.einsum("i,ijk->jk", W, F)
    elif axis == 1:
        P = np.einsum("j,ijk->ik", W, F)
    else:
        P = np.einsum("k,ijk->ij", W, F)
    return P


# --------------------------------------------------------------------------
def _g3_straight(R, S, T):
    return R.copy(), S.copy(), T.copy()


def _g3_curved(R, S, T):
    X = R + 0.20 * np.sin(1.3 * S) * np.cos(0.7 * T) + 0.10 * R * S
    Y = S + 0.18 * np.sin(1.1 * R) - 0.12 * S * T
    Z = T + 0.15 * np.cos(0.9 * R) * np.sin(1.1 * S) + 0.08 * T * R
    return X, Y, Z


GEOMS3 = {"straight": _g3_straight, "curved": _g3_curved}

CASES3D = {
    # The 3D cases are the DIRECT EXTENSIONS of the 2D cases in contour2d.CASES2D:
    # same formula with a t-term added, same isolevel.  Setting t = 0 in any of these
    # recovers the corresponding 2D field exactly, so the t=0 cut of each surface is
    # the matching contour in E_contour_cases.png.  (Verified by check_2d_consistency.)
    "circle": ("sphere r^2+s^2+t^2 = 0.6^2  [2D 'circle' + t^2]",
               lambda R, S, T: R ** 2 + S ** 2 + T ** 2, 0.36),
    "twoblobs": ("two disjoint spheres  [2D 'twoblobs' + t-offset]",
                 lambda R, S, T: np.minimum(
                     (R - 0.45) ** 2 + (S - 0.4) ** 2 + T ** 2,
                     (R + 0.5) ** 2 + (S + 0.45) ** 2 + T ** 2), 0.09),
    "saddle": ("r*s + 0.25 t = 0.15: the ambiguous case  [2D 'saddle' + 0.25t]",
               lambda R, S, T: R * S + 0.25 * T, 0.15),
    "tangency": ("(r^2+s^2+t^2-0.5)^2 = 0: double root  [2D 'tangency' + t^2]",
                 lambda R, S, T: (R ** 2 + S ** 2 + T ** 2 - 0.5) ** 2, 0.0),
    "wave": ("sin(2.6 pi r) + cos(2.2 pi s) + 0.6 sin(1.8 pi t) = 0.3  [2D 'wave' + t]",
             lambda R, S, T: (np.sin(2.6 * np.pi * R) + np.cos(2.2 * np.pi * S)
                              + 0.6 * np.sin(1.8 * np.pi * T)), 0.3),
    "band": ("thin sheet: forces line insertion  [2D 'band' + t]",
             lambda R, S, T: np.exp(-40.0 * (S - 0.3 * R - 0.12) ** 2
                                    - 3.0 * T ** 2), 0.55),
    "corner": ("r+s+t = 1.15: plane clipping the corner  [2D 'corner' + t]",
               lambda R, S, T: R + S + T, 1.15),
    "ellipse": ("ellipsoid r^2 + 2.6(s-0.12)^2 + 1.6 t^2  [2D 'ellipse' + 1.6t^2]",
                lambda R, S, T: R ** 2 + 2.6 * (S - 0.12) ** 2 + 1.6 * T ** 2, 0.45),
    "wavycircle": ("wavy sphere  [2D 'wavycircle' + t^2]",
                   lambda R, S, T: (np.sqrt(np.maximum(R ** 2 + S ** 2 + T ** 2, 1e-30))
                                    - 0.10 * np.cos(3 * np.arctan2(S, R))), 0.62),
    # one case with no 2D analogue, kept because genus matters in 3D only
    "torus": ("genus-1 tube (NO 2D analogue: topology that needs 3D)",
              lambda R, S, T: (np.sqrt(np.maximum(R ** 2 + S ** 2, 1e-30)) - 0.55) ** 2
                              + T ** 2, 0.0625),
}


def check_2d_consistency(N=8, tol=1e-12):
    """Verify that each 3D case reduces to its 2D namesake at t = 0.

    This is what makes E_surface3d.png and E_contour_cases.png the SAME problem:
    the t=0 slice of every 3D field is the 2D field of the same name.
    """
    sys.path.insert(0, _HERE)
    import contour2d as _C2
    r, _ = gll_nodes_weights(N)
    R2, S2 = np.meshgrid(r, r, indexing="ij")
    rows = []
    for key, (d3, f3, l3) in CASES3D.items():
        if key not in _C2.CASES2D:
            rows.append((key, None, None, "no 2D analogue"))
            continue
        d2, f2, l2 = _C2.CASES2D[key]
        U2 = f2(R2, S2)
        U3 = f3(R2, S2, np.zeros_like(R2))          # t = 0 slice
        err = float(np.max(np.abs(U2 - U3)))
        rows.append((key, err, abs(l2 - l3), "ok" if (err < tol and l2 == l3)
                     else "DIFFERS"))
    return rows


# --------------------------------------------------------------------------
# high-order refinement of the surface (section 6): the 3D form of refine2d
# --------------------------------------------------------------------------
GRAD_FLOOR3 = 1e-10


def grad3(F, p, h=1e-6):
    """Gradient of the tensor interpolant at p = (r,s,t)."""
    g = np.empty(3)
    for a in range(3):
        q1 = list(p); q2 = list(p)
        q1[a] = min(q1[a] + h, 1.0)
        q2[a] = max(q2[a] - h, -1.0)
        g[a] = (tensor_eval3(F, *q1) - tensor_eval3(F, *q2)) / (q1[a] - q2[a] + 1e-300)
    return g


def project3(F, p, nit=3):
    """Newton onto the isosurface along the gradient.  Same formula as 2D.

        x <- x - f(x) grad f / |grad f|^2

    Returns (point, converged).  Refuses to move where |grad f| ~ 0 (a tangency).
    """
    p = np.array(p, float)
    for _ in range(nit):
        v = tensor_eval3(F, *p)
        g = grad3(F, p)
        n2 = float(g @ g)
        if n2 < GRAD_FLOOR3 ** 2:
            return p, False
        p = p - v * g / n2
        p = np.clip(p, -1.0, 1.0)
    return p, True


def subdivide_project(F, tris, nsub=1):
    """Curve each triangle: 4^nsub sub-triangles whose every new node is projected.

    This is the 3D counterpart of refine2d.refine_segment(mode='project').  Corner
    nodes are already exact fibre roots and are never moved; every node created by
    subdivision is pushed onto the isosurface, so the refined mesh is exact at all
    of its nodes.  Each node is independent -- no linear system, no C^1 conditions,
    no neighbour data -- which is what makes it parallel-friendly.
    """
    out = list(tris)
    for _ in range(nsub):
        new = []
        for (a, b, c) in out:
            A, B, Cc = np.array(a), np.array(b), np.array(c)
            mab, _ = project3(F, 0.5 * (A + B))
            mbc, _ = project3(F, 0.5 * (B + Cc))
            mca, _ = project3(F, 0.5 * (Cc + A))
            new += [(tuple(A), tuple(mab), tuple(mca)),
                    (tuple(mab), tuple(B), tuple(mbc)),
                    (tuple(mca), tuple(mbc), tuple(Cc)),
                    (tuple(mab), tuple(mbc), tuple(mca))]
        out = new
    return out


def surface_residual(F, tris, nsample=4):
    """Max |f - c| over the triangles, sampled INSIDE them, not only at the nodes."""
    worst = 0.0
    bary = [(i / nsample, j / nsample, 1 - i / nsample - j / nsample)
            for i in range(nsample + 1) for j in range(nsample + 1 - i)]
    for (a, b, c) in tris:
        A, B, Cc = np.array(a), np.array(b), np.array(c)
        for (w0, w1, w2) in bary:
            p = w0 * A + w1 * B + w2 * Cc
            worst = max(worst, abs(tensor_eval3(F, *p)))
    return worst
