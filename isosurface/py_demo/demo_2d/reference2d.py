"""
reference2d.py -- the two conventional approaches, for comparison.

Both are what a practitioner would actually do today, and both are implemented
here WITHOUT a third-party contouring library: marching squares with LINEAR edge
interpolation is ~40 lines, and writing it ourselves makes the comparison honest
(same topology logic, the only difference is where the vertices come from).

    R1  GLL -> N^2 linear cells -> marching squares with linear interpolation
        (i.e. treat the GLL nodal values as a bilinear grid: what a standard
         visualisation tool does when handed the element's nodes)

    R2  GLL -> interpolate to a uniform MxM grid (M ~ 1.5N) -> same marching squares
        (i.e. resample the high-order field, then contour the resampled data:
         what a tool does when told to 'refine' a high-order element)

Both produce vertices that are LINEARLY interpolated along cell edges, so their
error is the linear-interpolation error of the crossing, not machine precision.
"""
import numpy as np
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from kernel2d import gll_nodes_weights, lagrange_weights      # noqa: E402

__all__ = ["marching_squares_linear", "reference_R1", "reference_R2",
           "resample_uniform"]


# the same case table as the exact method, so only the VERTEX source differs
_MS = {0: [], 15: [],
       1: [(3, 0)], 14: [(3, 0)],
       2: [(0, 1)], 13: [(0, 1)],
       3: [(3, 1)], 12: [(3, 1)],
       4: [(1, 2)], 11: [(1, 2)],
       6: [(0, 2)], 9: [(0, 2)],
       7: [(3, 2)], 8: [(3, 2)],
       5: "saddle", 10: "saddle"}


def marching_squares_linear(gr, gs, P):
    """Marching squares on values P[i,j] at grid (gr, gs), LINEAR edge interpolation.

    The saddle cases are resolved by the standard bilinear asymptotic decider,
    which is the best a sign-based method can do without the underlying polynomial.
    """
    nr, ns = len(gr), len(gs)
    sg = np.sign(P)
    segs = []
    n_amb = 0

    def lerp(t0, t1, v0, v1):
        return t0 + (t1 - t0) * (v0 / (v0 - v1))

    def epoint(e, i, j):
        if e == 0:
            return (lerp(gr[i], gr[i + 1], P[i, j], P[i + 1, j]), gs[j])
        if e == 1:
            return (gr[i + 1], lerp(gs[j], gs[j + 1], P[i + 1, j], P[i + 1, j + 1]))
        if e == 2:
            return (lerp(gr[i], gr[i + 1], P[i, j + 1], P[i + 1, j + 1]), gs[j + 1])
        if e == 3:
            return (gr[i], lerp(gs[j], gs[j + 1], P[i, j], P[i, j + 1]))
        return None

    for i in range(nr - 1):
        for j in range(ns - 1):
            c = [sg[i, j], sg[i + 1, j], sg[i + 1, j + 1], sg[i, j + 1]]
            idx = sum((1 << k) for k in range(4) if c[k] < 0)
            ent = _MS[idx]
            if ent == []:
                continue
            if ent == "saddle":
                n_amb += 1
                # asymptotic decider: sign of the bilinear saddle value
                a00, a10 = P[i, j], P[i + 1, j]
                a11, a01 = P[i + 1, j + 1], P[i, j + 1]
                den = (a00 - a10 - a01 + a11)
                sad = (a00 * a11 - a10 * a01) / den if den != 0 else 0.0
                if (sad < 0) == (c[0] < 0):
                    pairs = [(3, 0), (1, 2)] if idx == 5 else [(0, 1), (2, 3)]
                else:
                    pairs = [(0, 1), (2, 3)] if idx == 5 else [(3, 0), (1, 2)]
            else:
                pairs = ent
            for (e0, e1) in pairs:
                p0, p1 = epoint(e0, i, j), epoint(e1, i, j)
                if p0 is not None and p1 is not None:
                    segs.append((p0, p1))
    return segs, n_amb


def reference_R1(N, U, level):
    """R1: contour the GLL nodal values directly as an N x N bilinear grid."""
    x, _ = gll_nodes_weights(N)
    segs, namb = marching_squares_linear(x, x, U - level)
    return dict(segments=segs, ambiguous=namb, grid=(x, x),
                cells=N * N, evals=(N + 1) ** 2, label="R1 GLL as bilinear")


def resample_uniform(N, U, M):
    """Interpolate the degree-N tensor field onto a uniform M x M grid."""
    t = np.linspace(-1.0, 1.0, M)
    W = np.array([lagrange_weights(N, ti) for ti in t])       # (M, N+1)
    return t, W @ U @ W.T


def reference_R2(N, U, level, M=None):
    """R2: resample to a uniform M x M grid (default M ~ 1.5N), then contour."""
    if M is None:
        M = int(round(1.5 * N)) + 1
    t, P = resample_uniform(N, U, M)
    segs, namb = marching_squares_linear(t, t, P - level)
    return dict(segments=segs, ambiguous=namb, grid=(t, t),
                cells=(M - 1) ** 2, evals=M * M,
                label=f"R2 uniform {M}x{M}")
