"""
contour2d.py -- 2D contour extraction on a single (possibly curved) GLL element.

This is the 3D isosurface strategy of section 6, one dimension down, where it can
be drawn and checked by eye.  The pipeline is identical:

    1. FIBER SOLVE   2(N+1) one-dimensional rootfinds on the reference element,
                     using route C (matrix-free Ehrlich-Aberth, demo_C).
    2. INSERTION     where a sub-cell edge holds >= 2 roots, insert a line so that
                     every edge carries at most one crossing (marching squares'
                     assumption).  The inserted node's SIGN is known a priori.
    3. MARCHING      marching squares on the N^2 sub-cells, reading exact crossing
                     locations from step 1 instead of interpolating.
    4. AMBIGUITY     the one ambiguous marching-squares case (saddle, case 5/10) is
                     decided by EVALUATING the interpolant at the cell centre --
                     no lookup table, no asymptotic decider.
    5. MAP OUT       reference (r,s) -> physical (x,y) by barycentric interpolation
                     of the element's own geometry nodes.

Everything the topology needs -- crossing location, per-edge count, sign of p' --
comes out of the fiber solve.  No field evaluation is needed between roots: between
two consecutive roots the sign is constant, so the full sign structure follows from
one nodal sign plus the parity of each root.

Exports
    gll_grid, tensor_eval, ContourResult
    solve_fibers, insert_lines, marching_squares, extract_contour
    map_to_physical, CASES2D
"""
import numpy as np
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

# SELF-CONTAINED: every numerical routine is in kernel2d.py, written from the
# defining recurrences.  No numpy.polynomial, no scipy, no LAPACK, and no import
# from the 1D demos -- this directory stands alone and ports directly to C.
from kernel2d import (gll_nodes_weights, bary_weights, nodal_to_modal,   # noqa: E402
                      lagrange_weights as _lw, legendre_val_der,
                      aberth_roots, bisect_bracketed, deflate,
                      A_FLOOR, NITER)

__all__ = ["gll_grid", "tensor_eval", "lagrange_weights", "ContourResult",
           "solve_fibers", "insert_lines", "marching_squares",
           "extract_contour", "map_to_physical", "CASES2D", "GEOMS"]


# --------------------------------------------------------------------------
# reference element utilities
# --------------------------------------------------------------------------
def gll_grid(N):
    """(r, R, S) : 1D nodes and the (N+1)x(N+1) tensor grid, indexing='ij'."""
    r, _ = gll_nodes_weights(N)
    R, S = np.meshgrid(r, r, indexing="ij")
    return r, R, S


lagrange_weights = _lw


def tensor_eval(U, rr, ss):
    """Evaluate the degree-N tensor interpolant of nodal values U at (rr, ss)."""
    N = U.shape[0] - 1
    return float(lagrange_weights(N, rr) @ U @ lagrange_weights(N, ss))


# --------------------------------------------------------------------------
# step 1: fiber solves
# --------------------------------------------------------------------------
class ContourResult:
    """Everything the topology needs, gathered from the fiber solves."""

    def __init__(self, N, U, level):
        self.N = N
        self.U = U
        self.level = level
        self.F = U - level                 # nodal values of f - c
        # crossings[d][i] = sorted roots on fiber i in direction d
        #   d = 0 : r varies, s = r_node[i]      d = 1 : s varies, r = r_node[i]
        self.crossings = [[], []]
        self.slopes = [[], []]             # p'(root), for orientation / tangency
        self.inserted = []                 # (direction, coordinate) of inserted lines
        self.n_solves = 0
        self.n_edge_newton = 0
        self.segments = []                 # list of ((r0,s0),(r1,s1)) in reference coords
        self.ambiguous = 0                 # saddle cells resolved by evaluation

    # ---- diagnostics ----
    def n_roots(self):
        return sum(len(v) for d in self.crossings for v in d)

    def multi_root_edges(self, grid):
        """count sub-cell edges holding >= 2 roots, against a given 1D grid"""
        n = 0
        for d in (0, 1):
            for rts in self.crossings[d]:
                for e in range(len(grid) - 1):
                    if np.sum((rts > grid[e]) & (rts < grid[e + 1])) >= 2:
                        n += 1
        return n


def solve_fibers(N, U, level, niter=NITER, a_floor=A_FLOOR):
    """Step 1.  Solve all 2(N+1) fibers with route C.  Returns a ContourResult.

    Also records p'(root): its SIGN gives the crossing direction and its MAGNITUDE
    flags tangency.  Both are by-products of the Aberth iteration, not extra work.
    """
    res = ContourResult(N, U, level)
    F = res.F
    for d in (0, 1):
        for i in range(N + 1):
            fib = F[:, i] if d == 0 else F[i, :]
            rts, sl = aberth_roots(N, fib, niter=niter, a_floor=a_floor,
                                   want_slope=True)
            res.n_solves += 1
            res.crossings[d].append(rts)
            res.slopes[d].append(sl)
    return res


# --------------------------------------------------------------------------
# step 2: adaptive line insertion
# --------------------------------------------------------------------------
def insert_lines(res, grid=None):
    """Step 2.  Guarantee at most one root per sub-cell edge.

    Where a sub-cell edge of the current grid holds >= 2 roots, insert a grid line
    midway between the first two.  The inserted node's value is obtained by ONE
    tensor contraction (O(N^2) flops for the whole line), and the sign there is
    known a priori -- between two consecutive roots the sign is the opposite of the
    bracketing nodes' common sign -- so no root information is lost.

    Returns the refined 1D grids (gr, gs); the fiber roots are unchanged, since
    inserting a line does not move any root.
    """
    N = res.N
    base, _ = gll_nodes_weights(N)
    grids = [base.copy(), base.copy()]

    for d in (0, 1):
        extra = []
        for rts in res.crossings[d]:
            g_ = grids[d]
            for e in range(len(g_) - 1):
                inr = rts[(rts > g_[e]) & (rts < g_[e + 1])]
                if len(inr) >= 2:
                    extra.append(0.5 * (inr[0] + inr[1]))
        if extra:
            grids[d] = np.unique(np.concatenate([grids[d], np.array(extra)]))
            for m in extra:
                res.inserted.append((d, float(m)))
    return grids[0], grids[1]


def _value_on_grid(res, gr, gs):
    """Nodal values of f-c on the (possibly refined) tensor grid, by contraction."""
    N = res.N
    Wr = np.array([lagrange_weights(N, t) for t in gr])      # (nr, N+1)
    Ws = np.array([lagrange_weights(N, t) for t in gs])      # (ns, N+1)
    return Wr @ res.F @ Ws.T


# --------------------------------------------------------------------------
# steps 3-4: marching squares with exact crossings
# --------------------------------------------------------------------------
#  corner numbering of a cell (i,j):      edge numbering:
#     3 --- 2         s+                     --2--
#     |     |         ^                     3|   |1
#     0 --- 1         +-> r                   --0--
#  corner k sign -> bit k of the case index (bit set when value < 0)
#
#  edge e connects corners (e, (e+1) % 4); edge 0 is the -s edge (r varies), etc.
_MS_EDGES = {
    0:  [], 15: [],
    1:  [(3, 0)], 14: [(3, 0)],
    2:  [(0, 1)], 13: [(0, 1)],
    3:  [(3, 1)], 12: [(3, 1)],
    4:  [(1, 2)], 11: [(1, 2)],
    6:  [(0, 2)],  9: [(0, 2)],
    7:  [(3, 2)],  8: [(3, 2)],
    # 5 and 10 are the AMBIGUOUS saddle cases: resolved by evaluation, see below
    5:  "saddle", 10: "saddle",
}


def marching_squares(res, gr, gs):
    """Steps 3-4.  March the refined lattice, reading exact crossings.

    Edge crossings are looked up from the fiber solves (not interpolated), so every
    emitted vertex lies on the true contour to machine precision.  The two ambiguous
    cases are decided by evaluating the interpolant at the cell centre.
    """
    N = res.N
    P = _value_on_grid(res, gr, gs)
    nr, ns = len(gr), len(gs)
    sign = np.sign(P)
    segs = []
    res.ambiguous = 0

    def edge_point(kind, i, j):
        """Exact crossing on a cell edge.  kind='r': between (i,j)-(i+1,j)."""
        if kind == "r":
            lo, hi, fixed = gr[i], gr[i + 1], gs[j]
            # the r-fiber at s = gs[j]; if gs[j] is an inserted line it is not an
            # original fiber, so fall back to a bracketed solve on the interpolant
            k = np.where(np.abs(gll_nodes_weights(N)[0] - fixed) < 1e-13)[0]
            if len(k):
                rts = res.crossings[0][k[0]]
                inr = rts[(rts >= lo - 1e-12) & (rts <= hi + 1e-12)]
                if len(inr) == 1:
                    return (float(inr[0]), fixed)
            t = _bracketed(res, "r", lo, hi, fixed)
            return (t, fixed)
        else:
            lo, hi, fixed = gs[j], gs[j + 1], gr[i]
            k = np.where(np.abs(gll_nodes_weights(N)[0] - fixed) < 1e-13)[0]
            if len(k):
                rts = res.crossings[1][k[0]]
                inr = rts[(rts >= lo - 1e-12) & (rts <= hi + 1e-12)]
                if len(inr) == 1:
                    return (fixed, float(inr[0]))
            t = _bracketed(res, "s", lo, hi, fixed)
            return (fixed, t)

    for i in range(nr - 1):
        for j in range(ns - 1):
            c = [sign[i, j], sign[i + 1, j], sign[i + 1, j + 1], sign[i, j + 1]]
            idx = sum((1 << k) for k in range(4) if c[k] < 0)
            entry = _MS_EDGES[idx]
            if entry == []:
                continue
            if entry == "saddle":
                # AMBIGUOUS: evaluate the interpolant at the cell centre.  This is
                # what replaces the asymptotic decider / MC33 subcase tables.
                res.ambiguous += 1
                rc = 0.5 * (gr[i] + gr[i + 1])
                sc = 0.5 * (gs[j] + gs[j + 1])
                centre = tensor_eval(res.F, rc, sc)
                # connect so that the centre's sign matches the enclosed region
                if (centre < 0) == (c[0] < 0):
                    pairs = [(3, 0), (1, 2)] if idx == 5 else [(0, 1), (2, 3)]
                else:
                    pairs = [(0, 1), (2, 3)] if idx == 5 else [(3, 0), (1, 2)]
            else:
                pairs = entry

            for (e0, e1) in pairs:
                p0 = _edge_lookup(edge_point, e0, i, j)
                p1 = _edge_lookup(edge_point, e1, i, j)
                if p0 is not None and p1 is not None:
                    segs.append((p0, p1))
    res.segments = segs
    return segs


def _edge_lookup(edge_point, e, i, j):
    """edge index -> (kind, cell offset)."""
    if e == 0:
        return edge_point("r", i, j)
    if e == 1:
        return edge_point("s", i + 1, j)
    if e == 2:
        return edge_point("r", i, j + 1)
    if e == 3:
        return edge_point("s", i, j)
    return None


def _bracketed(res, kind, lo, hi, fixed):
    """The single crossing on [lo,hi] of an INSERTED grid line.

    Inserted lines are not original fibers, so their crossing is not in the fiber
    solve.  There is a guaranteed sign change, so bisection cannot fail -- and this
    is the one place in the whole scheme where one-root-at-a-time solving is the
    right tool, because the bracket is known and exactly one root is wanted.
    """
    def val(t):
        return (tensor_eval(res.F, t, fixed) if kind == "r"
                else tensor_eval(res.F, fixed, t))
    return bisect_bracketed(val, lo, hi)


# --------------------------------------------------------------------------
# step 5: map the reference contour to physical space
# --------------------------------------------------------------------------
def map_to_physical(segments, Xn, Yn):
    """Reference (r,s) -> physical (x,y) via the element's own geometry nodes.

    The geometry is isoparametric: X and Y are stored as nodal values, so this is
    the same barycentric interpolation used for the field.  A root exact in (r,s)
    is therefore exact in (x,y) to the accuracy of the geometry representation --
    no extra error is introduced.
    """
    N = Xn.shape[0] - 1
    out = []
    for (p0, p1) in segments:
        q = []
        for (rr, ss) in (p0, p1):
            wr = lagrange_weights(N, rr)
            ws = lagrange_weights(N, ss)
            q.append((float(wr @ Xn @ ws), float(wr @ Yn @ ws)))
        out.append((q[0], q[1]))
    return out


# --------------------------------------------------------------------------
# the full pipeline
# --------------------------------------------------------------------------
def extract_contour(N, U, level, niter=NITER, insert=True):
    """Run the whole pipeline and return the ContourResult."""
    res = solve_fibers(N, U, level, niter=niter)
    if insert:
        gr, gs = insert_lines(res)
    else:
        gr, _ = gll_nodes_weights(N)
        gs = gr
    res.grid = (gr, gs)
    marching_squares(res, gr, gs)
    return res


# --------------------------------------------------------------------------
# geometries and test cases
# --------------------------------------------------------------------------
def _geom_straight(R, S):
    return R.copy(), S.copy()


def _geom_curved(R, S):
    """A smooth, invertible deformation of the reference square."""
    X = R + 0.22 * np.sin(1.4 * S) * np.cos(0.6 * R) + 0.10 * R * S
    Y = S + 0.18 * np.sin(1.2 * R) - 0.12 * S * S
    return X, Y


def _geom_annulus(R, S):
    """Polar-like map: a curved 'pie slice', strongly non-affine."""
    rad = 1.4 + 0.5 * S
    th = 0.55 * np.pi * (R * 0.5)
    return rad * np.cos(th) - 1.0, rad * np.sin(th) - 0.55


GEOMS = {"straight": _geom_straight,
         "curved": _geom_curved,
         "annulus": _geom_annulus}


CASES2D = {
    # key: (description, field on the REFERENCE square, level)
    "circle": ("single closed loop: r^2+s^2 = 0.6^2",
               lambda R, S: R ** 2 + S ** 2, 0.36),
    "twoblobs": ("two disjoint components (tests component labelling)",
                 lambda R, S: np.minimum((R - 0.45) ** 2 + (S - 0.4) ** 2,
                                         (R + 0.5) ** 2 + (S + 0.45) ** 2), 0.09),
    "saddle": ("hyperbola r*s = 0.15: the ambiguous marching-squares case",
               lambda R, S: R * S, 0.15),
    "tangency": ("level tangent to the field: double root, sqrt(eps) accuracy",
                 lambda R, S: (R ** 2 + S ** 2 - 0.5) ** 2, 0.0),
    "wave": ("oscillatory: several crossings per fiber",
             lambda R, S: np.sin(2.6 * np.pi * R) + np.cos(2.2 * np.pi * S), 0.3),
    "band": ("thin band: two roots inside one sub-cell edge (forces insertion)",
             lambda R, S: np.exp(-40.0 * (S - 0.3 * R - 0.12) ** 2), 0.55),
    "corner": ("contour clipping the element corner (open arcs)",
               lambda R, S: R + S, 1.15),
    "ellipse": ("elongated closed loop, off-centre",
                lambda R, S: R ** 2 + 2.6 * (S - 0.12) ** 2, 0.45),
    "wavycircle": ("r = 0.62 + 0.10 cos(3 theta): curvature varies along the loop",
                   lambda R, S: (np.hypot(R, S)
                                 - 0.10 * np.cos(3 * np.arctan2(S, R))), 0.62),
}
