"""
refine2d.py -- high-order refinement of the extracted contour.

The marching step emits straight segments between exact crossings.  The vertices are
at machine precision, but the segment BETWEEN them is a chord, and that chord is the
only remaining error.  Two ways to remove it, both of which use information the
polynomial already provides -- neither is a fit:

  (H)  CUBIC HERMITE.  The gradient of the interpolant is exact everywhere, and the
       contour is perpendicular to it, so the unit TANGENT at each vertex is exact
       too.  Endpoints + exact tangents = cubic Hermite, C^1 at every vertex.
       No linear system, no smoothing parameter.

  (P)  HERMITE + PROJECTION.  Subdivide the Hermite curve and push each new point
       onto the contour by Newton along the gradient,
             x <- x - f(x) grad f / |grad f|^2 .
       The new points are then exact roots, not approximations, so the curve is
       exact at every node.  The formula is dimension-independent: the same line
       curves a triangle in 3D (section 6).

Projection needs |grad f| bounded away from zero.  At a tangency grad f -> 0 and the
step blows up; `project` refuses to move in that case and the caller falls back to
Hermite, then to the chord.
"""
import numpy as np
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from contour2d import tensor_eval                                 # noqa: E402

__all__ = ["grad2", "unit_tangent", "hermite_points", "project",
           "refine_segment", "refine_contour", "curve_error", "GRAD_FLOOR"]

GRAD_FLOOR = 1e-10          # below this |grad f| the projection is not attempted


def grad2(F, rr, ss, h=1e-6):
    """Gradient of the interpolant.  Central differences on the tensor evaluation.

    (A production kernel would differentiate the Lagrange basis directly, at the same
    O(N) cost and without the step-size choice; central differences keep this demo
    short and are accurate to ~1e-10, far below the errors being compared.)
    """
    return np.array([
        (tensor_eval(F, rr + h, ss) - tensor_eval(F, rr - h, ss)) / (2 * h),
        (tensor_eval(F, rr, ss + h) - tensor_eval(F, rr, ss - h)) / (2 * h)])


def unit_tangent(F, p, ref=None):
    """Exact unit tangent to the contour at p: perpendicular to grad f."""
    g = grad2(F, p[0], p[1])
    n = np.hypot(g[0], g[1])
    if n < GRAD_FLOOR:
        return None
    t = np.array([-g[1], g[0]]) / n
    if ref is not None and t @ ref < 0:
        t = -t
    return t


def hermite_points(p0, t0, p1, t1, u):
    """Cubic Hermite through p0,p1 with unit tangents t0,t1, scaled by the chord."""
    p0 = np.asarray(p0, float)
    p1 = np.asarray(p1, float)
    L = float(np.linalg.norm(p1 - p0))
    u = np.asarray(u, float)[:, None]
    h00 = 2 * u ** 3 - 3 * u ** 2 + 1
    h10 = u ** 3 - 2 * u ** 2 + u
    h01 = -2 * u ** 3 + 3 * u ** 2
    h11 = u ** 3 - u ** 2
    return h00 * p0 + h10 * (L * t0) + h01 * p1 + h11 * (L * t1)


def project(F, p, nit=3):
    """Newton onto the contour along the gradient.  Returns (point, converged)."""
    p = np.array(p, float)
    for _ in range(nit):
        v = tensor_eval(F, p[0], p[1])
        g = grad2(F, p[0], p[1])
        n2 = float(g @ g)
        if n2 < GRAD_FLOOR ** 2:
            return p, False
        p = p - v * g / n2
        p = np.clip(p, -1.0, 1.0)
    return p, True


def refine_segment(F, p0, p1, k=3, mode="project"):
    """One straight segment -> k+2 points along a high-order curve.

    mode = 'straight' | 'hermite' | 'project'
    Returns an array of points including both endpoints (which are never moved:
    they are already exact roots from the fibre solve).
    """
    u = np.linspace(0.0, 1.0, k + 2)
    if mode == "straight":
        return (1 - u)[:, None] * np.asarray(p0) + u[:, None] * np.asarray(p1)

    d = np.asarray(p1, float) - np.asarray(p0, float)
    t0 = unit_tangent(F, p0, ref=d)
    t1 = unit_tangent(F, p1, ref=d)
    if t0 is None or t1 is None:                     # tangency: fall back
        return (1 - u)[:, None] * np.asarray(p0) + u[:, None] * np.asarray(p1)

    P = hermite_points(p0, t0, p1, t1, u)
    if mode == "hermite":
        return P
    for i in range(1, len(P) - 1):                   # endpoints stay exact
        q, ok = project(F, P[i])
        if ok:
            P[i] = q
    return P


def refine_contour(F, segments, k=3, mode="project"):
    """Refine every segment; returns a list of point arrays (one per segment)."""
    return [refine_segment(F, p0, p1, k=k, mode=mode) for (p0, p1) in segments]


def curve_error(F, segments, dist_fn, k=3, mode="project", nsample=17):
    """Max distance of the WHOLE curve from the true level set.

    The curve is sampled at `nsample` points per segment -- including between the
    refinement nodes, so a method cannot score well merely by placing its nodes on
    the contour while bulging in between.
    """
    worst = 0.0
    u = np.linspace(0.0, 1.0, nsample)
    for (p0, p1) in segments:
        P = refine_segment(F, p0, p1, k=k, mode=mode)
        # sample the piecewise-linear curve through the refined nodes densely
        segl = np.cumsum(np.r_[0.0, np.linalg.norm(np.diff(P, axis=0), axis=1)])
        if segl[-1] <= 0:
            continue
        s = u * segl[-1]
        X = np.interp(s, segl, P[:, 0])
        Y = np.interp(s, segl, P[:, 1])
        worst = max(worst, float(np.max(dist_fn(np.column_stack([X, Y])))))
    return worst


# --------------------------------------------------------------------------
# STORAGE-BUDGETED RECONSTRUCTION
#
# Two phases, and the distinction matters:
#
#   EXTRACTION    the element's nodal data is available.  Fibre solves, marching,
#                 gradients, and projection all happen here and may evaluate the
#                 interpolant freely.
#   RECONSTRUCTION  the element data is GONE.  The curve must be rebuilt from a
#                 stored payload alone -- this is the real use case, and it is what
#                 limits how much the curve can be curved.
#
# So the question is not "may we project?" (we may -- at extraction time) but
# "what must be STORED, and what is the best curve that payload supports?".
#
# Payloads, per segment, counting one marginal vertex (vertices are shared):
#   vertices only            2 floats   -> polyline
#   + tangent angle          3 floats   -> cubic Hermite        (C^1 at every vertex)
#   + 1 projected point      4 floats   -> quadratic through 3 points
#   + tangent + 1 projected  5 floats   -> quintic: C^1 ends AND interpolates the mid
#
# The tangent costs ONE float, not two: it is a unit vector, so an angle suffices.
# Projection with k -> infinity is the LOWER BOUND (error -> 0 at unbounded storage)
# and is kept as the reference curve, not as a storage scheme.
# --------------------------------------------------------------------------
def tangent_from_stored(F, p, along, slope):
    """The exact gradient at a contour vertex, WITHOUT evaluating the field anew.

    `along` is 0 if the vertex lies on an r-fibre (s fixed) and 1 for an s-fibre.
    `slope` is p'(root) from the Aberth iteration -- the along-fibre derivative,
    already computed.  The transverse derivative comes from ONE Lagrange contraction
    of the stored nodal values followed by a Legendre derivative: a weighted sum of
    data already in hand, not a new field evaluation.

    This is the routine that makes cubic Hermite affordable under the storage
    constraint: at extraction time it is nearly free, and only ONE float (the angle)
    need be stored.
    """
    from kernel2d import (lagrange_weights, nodal_to_modal, legendre_val_der)
    N = F.shape[0] - 1
    if along == 0:                               # r varies, s = p[1] fixed
        g_line = lagrange_weights(N, p[0]) @ F   # f(r*, s_j) for all j
        a = nodal_to_modal(N) @ g_line
        _, dv = legendre_val_der(a, np.array([p[1] + 0j]))
        g = np.array([slope, dv[0].real])
    else:                                        # s varies, r = p[0] fixed
        g_line = F @ lagrange_weights(N, p[1])   # f(r_i, s*) for all i
        a = nodal_to_modal(N) @ g_line
        _, dv = legendre_val_der(a, np.array([p[0] + 0j]))
        g = np.array([dv[0].real, slope])
    n = float(np.hypot(g[0], g[1]))
    if n < GRAD_FLOOR:
        return None
    return np.array([-g[1], g[0]]) / n


def _bezier(P, m):
    from math import comb
    P = np.asarray(P, float)
    n = len(P) - 1
    u = np.linspace(0.0, 1.0, m)[:, None]
    Q = np.zeros((m, 2))
    for i in range(n + 1):
        Q += comb(n, i) * (u ** i) * ((1 - u) ** (n - i)) * P[i]
    return Q


def reconstruct(payload, m=65):
    """RECONSTRUCTION: rebuild one segment's curve from its stored payload only.

    payload = dict(p0=, p1=, t0=?, t1=?, pm=?)   -- no field access whatsoever.
    """
    p0 = np.asarray(payload["p0"], float)
    p1 = np.asarray(payload["p1"], float)
    t0 = payload.get("t0")
    t1 = payload.get("t1")
    pm = payload.get("pm")
    u = np.linspace(0.0, 1.0, m)[:, None]

    if t0 is None and pm is None:                       # 2 floats: polyline
        return (1 - u) * p0 + u * p1
    if t0 is not None and pm is None:                   # 3 floats: cubic Hermite
        return hermite_points(p0, t0, p1, t1, u.ravel())
    if t0 is None and pm is not None:                   # 4 floats: quadratic
        c = 2 * np.asarray(pm, float) - 0.5 * (p0 + p1)
        return _bezier([p0, c, p1], m)
    # 5 floats: C^1 at both ends AND interpolating the midpoint
    H = hermite_points(p0, t0, p1, t1, u.ravel())
    Hm = hermite_points(p0, t0, p1, t1, np.array([0.5]))[0]
    bump = 16.0 * (u ** 2) * ((1 - u) ** 2)
    return H + bump * (np.asarray(pm, float) - Hm)


def build_payload(F, p0, p1, along0=None, slope0=None, along1=None, slope1=None,
                  want_tangent=True, want_mid=False):
    """EXTRACTION: assemble what must be stored for one segment."""
    pay = {"p0": tuple(p0), "p1": tuple(p1)}
    if want_tangent:
        d = np.asarray(p1, float) - np.asarray(p0, float)
        t0 = unit_tangent(F, p0, ref=d)
        t1 = unit_tangent(F, p1, ref=d)
        if t0 is not None and t1 is not None:
            pay["t0"], pay["t1"] = t0, t1
    if want_mid:
        u = np.array([0.5])
        if "t0" in pay:
            mid = hermite_points(p0, pay["t0"], p1, pay["t1"], u)[0]
        else:
            mid = 0.5 * (np.asarray(p0, float) + np.asarray(p1, float))
        q, ok = project(F, mid)
        if ok:
            pay["pm"] = tuple(q)
    return pay


def payload_floats(pay):
    """Marginal storage per segment: one vertex (shared) + extras.

    A unit tangent is stored as ONE angle, not two components.
    """
    n = 2.0                                   # the marginal vertex
    if "t0" in pay:
        n += 1.0                              # one angle per vertex
    if "pm" in pay:
        n += 2.0
    return n
