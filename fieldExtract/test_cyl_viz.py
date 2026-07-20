"""Visualize the cylinder demo (Plan G).

Usage:  python test_cyl_viz.py [field]     # field defaults to "vz" (also "temp")

Two panels:
  left  -- raw r-theta point cloud in the x-y plane at a mid-z slice, colored by
           the sampled field (cyldemo*.vts, Cartesian x=r cos t, y=r sin t);
  right -- azimuthally-averaged (r,z) profile (cyldemo_avgy_g*.vts, x=radius, y=0):
           pcolormesh of <field>_theta over (r, z).
Confirms the cylinder sampling + through-plugin theta average by eye. Falls back
to a synthetic analytic grid when the .vts files do not exist yet.
Writes test_cyl_viz.png; shows the figure if a display is available.
"""
import glob
import os
import sys

import numpy as np

import matplotlib

if not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt

import fe_read

# must match turbPipe_cyl.udf
NR, NT, NZ = 12, 48, 11
R0, R1, Z0, Z1 = 0.05, 0.45, 2.0, 4.0

FIELD = sys.argv[1] if len(sys.argv) > 1 else "vz"


def load_cloud():
    """(X, Y, F) of a mid-z slice of the raw cyl point cloud."""
    files = sorted(glob.glob("cyldemo[0-9]*.vts"))
    if files:
        d = fe_read.load(files[-1])
        nr, nt, nz = d["dims"]
        X = d["X"].reshape((nz, nt, nr))[nz // 2]
        Y = d["Y"].reshape((nz, nt, nr))[nz // 2]
        Z = d["Z"].reshape((nz, nt, nr))[nz // 2]
        F = np.asarray(d[FIELD], dtype=np.float64).reshape((nz, nt, nr))[nz // 2]
        return X, Y, Z, F, d["time"], files[-1]
    # synthetic fallback
    r = np.linspace(R0, R1, NR)
    th = np.arange(NT) * (2 * np.pi / NT)
    Rr, Th = np.meshgrid(r, th)  # (nt, nr)
    X, Y = Rr * np.cos(Th), Rr * np.sin(Th)
    Z = np.full_like(X, 0.5 * (Z0 + Z1))
    F = fe_read.ref_vz(X, Y, Z, 0.0) if FIELD == "vz" else fe_read.ref_T(X, Y, Z)
    return X, Y, Z, F, 0.0, "(synthetic)"


def load_profile():
    """(R, Z, F) of the azimuthally-averaged (r,z) profile."""
    files = sorted(glob.glob("cyldemo_avgy_g[0-9]*.vts"))
    if files:
        d = fe_read.load(files[-1])
        nr, ny, nz = d["dims"]
        R = d["X"].reshape((nz, ny, nr))[:, 0, :]  # x = radius
        Z = d["Z"].reshape((nz, ny, nr))[:, 0, :]
        F = np.asarray(d[FIELD], dtype=np.float64).reshape((nz, ny, nr))[:, 0, :]
        return R, Z, F, d["time"], files[-1]
    # synthetic fallback: analytic theta-average on the (r,z) grid
    r = np.linspace(R0, R1, NR)
    z = np.linspace(Z0, Z1, NZ)
    Z2, R2 = np.meshgrid(z, r, indexing="ij")  # (nz, nr)
    if FIELD == "vz":
        F = fe_read.ref_vz(R2, np.zeros_like(R2), Z2, 0.0)  # axisymmetric
    else:
        F = np.zeros_like(R2)  # <temp>_theta = 0
    return R2, Z2, F, 0.0, "(synthetic)"


Xc, Yc, Zc, Fc, tc, csrc = load_cloud()
Rp, Zp, Fp, tp, psrc = load_profile()

# relerr vs the analytic reference (like test_viz.py): the raw cloud vs the field
# itself, the profile vs its analytic azimuthal average (<vz>_theta=vz(r), <temp>_theta=0).
cloud_ref = fe_read.ref_vz(Xc, Yc, Zc, tc) if FIELD == "vz" else fe_read.ref_T(Xc, Yc, Zc)
err_c = fe_read.relerr(Fc, cloud_ref)
prof_ref = fe_read.ref_vz(Rp, np.zeros_like(Rp), Zp, tp) if FIELD == "vz" else np.zeros_like(Rp)
err_p = fe_read.relerr(Fp, prof_ref)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.5, 5.0))

sc = ax0.scatter(Xc, Yc, c=Fc, cmap="jet", s=28, edgecolors="k", linewidths=0.3)
ax0.set_aspect("equal")
fig.colorbar(sc, ax=ax0).set_label(FIELD)
ax0.set_xlabel("x")
ax0.set_ylabel("y")
ax0.set_title(f"{csrc}: raw r-theta cloud (mid-z), {FIELD}  t={tc:g}\nrelerr={err_c:.1e}")

pc = ax1.pcolormesh(Rp, Zp, Fp, shading="gouraud", cmap="jet")
ax1.scatter(Rp, Zp, s=6, c="k", alpha=0.4)
fig.colorbar(pc, ax=ax1).set_label(f"<{FIELD}>_theta")
ax1.set_xlabel("r")
ax1.set_ylabel("z")
ax1.set_title(f"{psrc}: azimuthal avg (r,z), t={tp:g}\nrelerr={err_p:.1e}")

fig.tight_layout()
fig.savefig("test_cyl_viz.png", dpi=150)
print("wrote test_cyl_viz.png")
if matplotlib.get_backend().lower() != "agg":
    plt.show()
