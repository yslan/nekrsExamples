"""Visualize the ADIOS2 .bp backend output (Plan F).

Usage:  python test_bp_viz.py [field]     # field defaults to "vz" (also "temp")

Loads the latest step of box3d.bp (feBox3d->process(..., "adios")) via the
adios2 python bindings and draws two panels of a mid-y x-z slice:
  left  -- the .bp field;
  right -- the paired box3d0000N.vts field (same interpolation), as a visual A/B.
Confirms the BP path opens and carries the right data. Falls back to a synthetic
analytic slice when box3d.bp does not exist yet (before the case has been run).
Writes test_bp_viz.png; shows the figure if a display is available.
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

X0 = (-0.3, -0.3, 2.0)  # must match turbPipe.udf::feBox3d
X1 = (0.3, 0.3, 4.0)
DIMS = (11, 7, 21)

FIELD = sys.argv[1] if len(sys.argv) > 1 else "vz"


def slice_from(d, field):
    """Mid-y x-z slice (nz, nx) of a loaded box dict."""
    nx, ny, nz = d["dims"]
    X = d["X"].reshape((nz, ny, nx))[:, ny // 2, :]
    Z = d["Z"].reshape((nz, ny, nx))[:, ny // 2, :]
    F = np.asarray(d[field], dtype=np.float64).reshape((nz, ny, nx))[:, ny // 2, :]
    return X, Z, F


def synthetic():
    nx, ny, nz = DIMS
    x = np.linspace(X0[0], X1[0], nx)
    y = np.linspace(X0[1], X1[1], ny)
    z = np.linspace(X0[2], X1[2], nz)
    Z3, Y3, X3 = np.meshgrid(z, y, x, indexing="ij")
    F3 = fe_read.ref_vz(X3, Y3, Z3, 0.0) if FIELD == "vz" else fe_read.ref_T(X3, Y3, Z3)
    return X3[:, ny // 2, :], Z3[:, ny // 2, :], F3[:, ny // 2, :]


bp = sorted(glob.glob("box3d.bp")) and os.path.isdir("box3d.bp")
panels = []
if bp:
    db = fe_read.load("box3d.bp")
    panels.append((f"box3d.bp (t={db['time']:g})", *slice_from(db, FIELD)))
    vts = sorted(glob.glob("box3d[0-9]*.vts"))
    if vts:
        dv = fe_read.load(vts[-1])
        panels.append((f"{vts[-1]} (t={dv['time']:g})", *slice_from(dv, FIELD)))
else:
    panels.append(("(synthetic analytic)", *synthetic()))

fig, axes = plt.subplots(1, len(panels), figsize=(6.2 * len(panels), 4.2), squeeze=False)
vmin = min(F.min() for _, _, _, F in panels)
vmax = max(F.max() for _, _, _, F in panels)
for ax, (title, X, Z, F) in zip(axes[0], panels):
    pc = ax.pcolormesh(Z, X, F, shading="gouraud", cmap="jet", vmin=vmin, vmax=vmax)
    ax.scatter(Z, X, s=4, c="k", alpha=0.4)
    fig.colorbar(pc, ax=ax).set_label(FIELD)
    ax.set_xlabel("Z")
    ax.set_ylabel("X")
    ax.set_title(f"{title}: mid-y {FIELD}")

fig.tight_layout()
fig.savefig("test_bp_viz.png", dpi=150)
print("wrote test_bp_viz.png")
if matplotlib.get_backend().lower() != "agg":
    plt.show()
