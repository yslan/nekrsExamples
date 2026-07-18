"""Pure-python quadrature tests for the GLL support (Plan D) -- no .vts needed.

1. zwgll sanity: n=5 nodes/weights vs known closed-form values; weight sums.
2. Consistency: GLL-weighted average vs trapezoidal average of the analytic
   fields (each rule on its own grid) agree within discretization tolerance.
3. Accuracy vs the EXACT analytic average: err(GLL) <= err(trap); vz is a
   polynomial (degree 6 in x,y), so GLL (exact for degree <= 2n-3) must hit it
   to ~machine precision while the trapezoid rule carries O(h^2) error.

Grids/extents mirror turbPipe.udf::feBox3d / feBox3dGll: {11,7,21} points on
x,y in [-0.3,0.3], z in [2,4].
"""
import sys

import numpy as np

import fe_read

R = fe_read.R  # 0.5
A = 0.3        # x,y half-width
Z0, Z1 = 2.0, 4.0
NX, NY, NZ = 11, 7, 21
TIME = 0.3

ok = True


def report(name, cond, detail=""):
    global ok
    ok = ok and cond
    print(f"{name}: {detail} -> {'PASS' if cond else 'FAIL'}")


# ---------------------------------------------------------------- 1. zwgll sanity
z5, w5 = fe_read.zwgll(5)
z5_ref = np.array([-1.0, -np.sqrt(3.0 / 7.0), 0.0, np.sqrt(3.0 / 7.0), 1.0])
w5_ref = np.array([1.0 / 10, 49.0 / 90, 32.0 / 45, 49.0 / 90, 1.0 / 10])
report("zwgll(5) nodes", np.abs(z5 - z5_ref).max() < 1e-12, f"maxerr={np.abs(z5 - z5_ref).max():.2e}")
report("zwgll(5) weights", np.abs(w5 - w5_ref).max() < 1e-12, f"maxerr={np.abs(w5 - w5_ref).max():.2e}")

sum_ok = True
for n in (2, 3, 7, 11, 21):
    z, w = fe_read.zwgll(n)
    sum_ok = sum_ok and abs(w.sum() - 2.0) < 1e-12 and np.abs(z + z[::-1]).max() < 1e-12
report("zwgll sums/symmetry", sum_ok, "sum(w)=2, z symmetric for n in {2,3,7,11,21}")


# ------------------------------------------------- grid averages of the two fields
def avg_xy(code):
    """Average vz over x,y in [-A,A]^2 using the rule given by code (0/1)."""
    x, wx = fe_read.axis_nodes_weights(NX, code, -A, A)
    y, wy = fe_read.axis_nodes_weights(NY, code, -A, A)
    X, Y = np.meshgrid(x, y, indexing="ij")
    F = fe_read.ref_vz(X, Y, np.zeros_like(X), TIME)
    return float(np.einsum("i,j,ij->", wx, wy, F))


def avg_z(code):
    """Average temp's z-factor cos(2.7 pi z + 0.2) over [Z0,Z1]."""
    z, wz = fe_read.axis_nodes_weights(NZ, code, Z0, Z1)
    return float(np.dot(wz, np.cos(2.7 * np.pi * z + 0.2)))


# exact averages
# E[x^{2k}] over [-A,A] is A^{2k}/(2k+1); vz = 6/5 (1 - (x^2+y^2)^3 / R^6) cos t
m2, m4, m6 = A**2 / 3, A**4 / 5, A**6 / 7
exact_xy = 6.0 / 5.0 * (1.0 - (2 * m6 + 6 * m4 * m2) / R**6) * np.cos(TIME)
c = 2.7 * np.pi
exact_z = (np.sin(c * Z1 + 0.2) - np.sin(c * Z0 + 0.2)) / (c * (Z1 - Z0))

for name, avg, exact in (("vz over xy", avg_xy, exact_xy), ("temp over z", avg_z, exact_z)):
    a_trap, a_gll = avg(0), avg(1)
    e_trap, e_gll = abs(a_trap - exact), abs(a_gll - exact)

    # 2. consistency: both rules approximate the same integral
    report(f"{name} consistency", abs(a_gll - a_trap) < 5e-2,
           f"|gll-trap|={abs(a_gll - a_trap):.2e}")

    # 3. accuracy ordering vs the exact average
    report(f"{name} accuracy", e_gll <= e_trap,
           f"err(gll)={e_gll:.2e} <= err(trap)={e_trap:.2e}")

# vz is polynomial (degree 6) => GLL quadrature is exact
e_poly = abs(avg_xy(1) - exact_xy)
report("vz GLL exactness", e_poly < 1e-12, f"err={e_poly:.2e} (polynomial, degree 6)")

sys.exit(0 if ok else 1)
