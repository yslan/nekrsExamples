"""
cases3d.py -- five distinctive HARD cases for isosurface extraction on one
spectral element, chosen from what the visualization literature treats as the
difficult configurations, plus one (TURBULENT) added from a production failure.

  (1) TUNNEL  f = rs + st + tr,  c = -0.18
      The Marching Cubes "case 13" configuration: the element's corner signs
      alternate, which is exactly the situation where face ambiguity and interior
      ambiguity interact.  Chernyaev's MC33 exists for this case, and resolving it
      from corner signs alone requires the asymptotic decider plus an interior test
      (Nielson-Hamann; Chernyaev; Lewiner et al.).  Here the interpolant is
      available, so the ambiguity is settled by EVALUATING it -- no subcase table.
      Measured: 18 ambiguous cells at N=8, corner signs (+,-,-,-,-,-,-,+).

  (2) GYROID  f = sin r cos s + sin s cos t + sin t cos r,  c = 0
      A triply-periodic minimal surface: several disconnected sheets threading one
      element, high genus, and many ambiguous cells (81 at N=8).  The standard
      stress test for whether a method reproduces topology rather than merely
      producing triangles.

  (3) TURBULENT  f = sum of three incommensurate mode pairs at w = 8,  c = 0.3, N = 16
      The case the other four do NOT cover: a fibre carrying TWO roots inside one
      GLL interval.  Both bracketing nodes then have the SAME sign, so no sign
      test can see the second crossing and marching cubes leaves a hole
      (notes/16).  The four cases above produce ZERO such pairs -- measured -- and
      that is exactly why the defect survived to production, where TGV at Q = 0.1
      loses 5-12% of its roots this way.
      Posed at N = 16, not 8, deliberately.  The 3D interpolation error of this
      field is 1.2 at N=8, 3.4e-03 at N=14 and 2.6e-04 at N=16: below N=16 the
      element does not resolve the field at all, every surface error is swamped
      by discretisation, and the case cannot show whether the multi-root
      handling works.  Raising N does NOT make the pairs go away -- which is the
      substantive point, since it means they are a property of the lattice and
      not of under-resolution.
      Measured at N=16: 3424 roots, 83 multi-root pairs (2.4%), and never more
      than 2 roots in one interval.

  (4) NECK  f = ((r-a)^2+s^2+t^2)((r+a)^2+s^2+t^2),  a = 0.46,  c = 0.045
      Two lobes at the instant they merge.  The isosurface is on the verge of a
      topology change, and the gradient nearly vanishes at the neck
      (min |grad f| = 2.0e-02, against O(1) elsewhere): the sqrt(eps)-conditioned
      tangency of section 2, in three dimensions.  Also the one case where
      gradient projection must be guarded.

  (5) TGV  the Q-criterion of the Taylor-Green initial condition, on one big
      undeformed element.  A real flow field rather than a manufactured one --
      but SMOOTH: at t = 0 the vortices are well resolved, so it produces no
      multi-root pairs either.  That is the point of keeping both it and
      TURBULENT: the same physics, resolved and under-resolved.

Each is posed on a DEFORMED element, except TGV.  The reference surface is the
analytic level set of the same formula, evaluated on a fine grid -- not another
discretisation of the same pipeline.
"""
import numpy as np

__all__ = ["CASES", "GEOM", "geom_map", "analytic_reference",
           "tgv_velocity", "tgv_q", "TGV_BOX"]


def _f_tunnel(R, S, T):
    return R * S + S * T + T * R


def _f_gyroid(R, S, T):
    return (np.sin(2.3 * R) * np.cos(2.3 * S)
            + np.sin(2.3 * S) * np.cos(2.3 * T)
            + np.sin(2.3 * T) * np.cos(2.3 * R))


# Three incommensurate mode pairs.  The frequencies are deliberately not
# rational multiples of each other: a commensurate field puts its extrema on a
# regular grid, which can align with the GLL nodes and hide the very pairs this
# case exists to produce.  w = 8 against N = 8 is roughly the under-resolution of
# the production TGV run (8 elements at N = 7 across [-pi,pi]^3).
_TURB_W = 8.0


def _f_turbulent(R, S, T, w=_TURB_W):
    return (np.sin(w * R) * np.cos(0.92 * w * S)
            + np.sin(1.07 * w * S) * np.cos(0.87 * w * T)
            + np.sin(0.95 * w * T) * np.cos(1.11 * w * R))


def _f_neck(R, S, T, a=0.46):
    return ((R - a) ** 2 + S ** 2 + T ** 2) * ((R + a) ** 2 + S ** 2 + T ** 2)


# --------------------------------------------------------------------------
# Taylor-Green vortex: a real flow field, on ONE big UNDEFORMED element
# --------------------------------------------------------------------------
TGV_BOX = (0.0, 2.0 * np.pi)          # the full TGV cell, [0,2pi]^3


def tgv_velocity(X, Y, Z):
    """The standard Taylor-Green initial condition."""
    u = np.sin(X) * np.cos(Y) * np.cos(Z)
    v = -np.cos(X) * np.sin(Y) * np.cos(Z)
    w = np.zeros_like(X)
    return u, v, w


def tgv_q(X, Y, Z):
    """Q-criterion, Q = 1/2(|Omega|^2 - |S|^2), the standard vortex identifier.

    Computed ANALYTICALLY from the Taylor-Green velocity, so it is a genuine flow
    diagnostic rather than a manufactured field -- the only case here that a
    practitioner would actually ask for.  With w = 0 the velocity gradient is

        du/dx =  cos x cos y cos z     du/dy = -sin x sin y cos z
        du/dz = -sin x cos y sin z     dv/dx =  sin x sin y cos z
        dv/dy = -cos x cos y cos z     dv/dz =  cos x sin y sin z
        dw/*  =  0

    and Q follows from the symmetric/antisymmetric split.
    """
    cx, sx = np.cos(X), np.sin(X)
    cy, sy = np.cos(Y), np.sin(Y)
    cz, sz = np.cos(Z), np.sin(Z)
    g = np.zeros(X.shape + (3, 3))
    g[..., 0, 0] = cx * cy * cz
    g[..., 0, 1] = -sx * sy * cz
    g[..., 0, 2] = -sx * cy * sz
    g[..., 1, 0] = sx * sy * cz
    g[..., 1, 1] = -cx * cy * cz
    g[..., 1, 2] = cx * sy * sz
    Q = np.zeros_like(X)
    for i in range(3):
        for j in range(3):
            Sij = 0.5 * (g[..., i, j] + g[..., j, i])
            Oij = 0.5 * (g[..., i, j] - g[..., j, i])
            Q += 0.5 * (Oij * Oij - Sij * Sij)
    return Q


def _f_tgv(R, S, T):
    """Q-criterion of the Taylor-Green vortex, as a function of reference coords."""
    lo, hi = TGV_BOX
    X = lo + (R + 1.0) * (hi - lo) / 2.0
    Y = lo + (S + 1.0) * (hi - lo) / 2.0
    Z = lo + (T + 1.0) * (hi - lo) / 2.0
    return tgv_q(X, Y, Z)


CASES = {
    "tunnel": dict(
        f=_f_tunnel, level=-0.18,
        title="MC case 13: the tunnel",
        why="alternating corner signs; face and interior ambiguity interact. "
            "Sign-based MC needs MC33's subcase tables; we evaluate the polynomial."),
    "gyroid": dict(
        f=_f_gyroid, level=0.0,
        title="gyroid: many sheets, high genus",
        why="a triply-periodic minimal surface: several disconnected sheets in one "
            "element and the largest number of ambiguous cells of any case here."),
    "tgv": dict(
        f=_f_tgv, level=0.15, undeformed=True, order=16,
        title="Taylor-Green vortex: Q-criterion",
        why="a REAL flow field, not a manufactured one: the Q-criterion of the "
            "standard Taylor-Green initial condition on one big undeformed element "
            "spanning the whole [0,2pi]^3 cell. Eight vortex tubes meet at the cell "
            "centre, so the isosurface has genuine flow topology; at N=16 one "
            "element resolves Q to 3.4e-07."),
    "turbulent": dict(
        f=_f_turbulent, level=0.3, order=16,
        title="turbulent: two roots in one interval",
        why="the multi-root case. A fibre carries two roots inside one GLL "
            "interval, the bracketing nodes share a sign, and a sign-based method "
            "silently drops the pair -- a hole, not a visible error. The other "
            "four cases produce ZERO such pairs, which is why this defect reached "
            "production. At N=16: 83 pairs in 3424 roots (2.4%), never more than "
            "2 in one interval, with the interpolant resolved to 2.6e-04 so the "
            "surface error is not swamped by discretisation."),
    "neck": dict(
        f=_f_neck, level=0.045,
        title="neck: a topology change in progress",
        why="two lobes at the moment of merging; |grad f| ~ 2e-2 at the neck, so "
            "this is the sqrt(eps) tangency of section 2 and the case where "
            "gradient projection must be guarded."),
}


def geom_map(R, S, T):
    """A smooth, invertible deformation of the reference cube.

    Deliberately non-affine in all three directions so that reference fibres are
    curves in physical space and the isoparametric map is doing real work.
    """
    X = R + 0.20 * np.sin(1.3 * S) * np.cos(0.7 * T) + 0.10 * R * S
    Y = S + 0.18 * np.sin(1.1 * R) - 0.12 * S * T
    Z = T + 0.15 * np.cos(0.9 * R) * np.sin(1.1 * S) + 0.08 * T * R
    return X, Y, Z


GEOM = geom_map


def case_order(key, default=8):
    """Element order for a case.  TGV needs a big element; the others do not."""
    return CASES[key].get("order", default)


def case_geom(key):
    """The geometry map for a case: identity where the case asks to be undeformed."""
    if CASES[key].get("undeformed", False):
        return lambda R, S, T: (R.copy(), S.copy(), T.copy())
    return geom_map


def analytic_reference(key, ngrid=64):
    """The TRUE level set of the analytic formula, on a fine uniform grid.

    This is the reference the extracted surface is compared against.  It is not
    produced by the pipeline under test: it is the analytic field sampled finely and
    contoured by a straightforward (linear) marching cubes, so its own error is
    O(h^2) with h = 2/ngrid -- small enough to serve as ground truth for the
    geometric comparison, and independent of everything being measured.
    """
    c = CASES[key]
    t = np.linspace(-1.0, 1.0, ngrid)
    R, S, T = np.meshgrid(t, t, t, indexing="ij")
    return t, c["f"](R, S, T) - c["level"]


def tgv_speed_ref(R, S, T):
    """|u| of the Taylor-Green velocity at reference coordinates -- used to colour
    the isosurface by a second physical quantity, as a practitioner would."""
    lo, hi = TGV_BOX
    X = lo + (R + 1.0) * (hi - lo) / 2.0
    Y = lo + (S + 1.0) * (hi - lo) / 2.0
    Z = lo + (T + 1.0) * (hi - lo) / 2.0
    u, v, w = tgv_velocity(X, Y, Z)
    return np.sqrt(u * u + v * v + w * w)


def describe(key):
    c = CASES[key]
    return (f"{c['title']}  [N={case_order(key)}"
            f"{', undeformed' if c.get('undeformed') else ', deformed'}]")
