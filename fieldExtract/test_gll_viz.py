"""Visualize the GLL grid via the box3dgll 1D-average output (Plan D).

Usage:  python test_gll_viz.py [field]     # field defaults to "vz" (also try "temp")

Loads the latest box3dgll_avgyPlane*.vts (feBox3dGll->doAvg("y"), gather mode:
y collapses -> x-z plane on GLL points) and draws a pcolormesh whose cell edges
sit at the node coordinates, with the nodes overlaid -- the edge clustering of
the Gauss-Lobatto-Legendre distribution is directly visible by eye. Falls back
to a synthetic GLL grid carrying the analytic y-average when the .vts does not
exist yet (i.e. before the NekRS case has been run).
Writes test_gll_viz.png; shows the figure if a display is available.
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

X0 = (-0.3, -0.3, 2.0)  # must match turbPipe.udf::feBox3dGll
X1 = (0.3, 0.3, 4.0)
DIMS = (11, 7, 21)

FIELD = sys.argv[1] if len(sys.argv) > 1 else "vz"


def load_plane(field):
    """(X, Z, F, time, pointDist, src) of the y-averaged x-z plane."""
    files = sorted(glob.glob("box3dgll_avgyPlane*.vts"))
    if files:
        d = fe_read.load(files[-1])
        nx, ny, nz = d["dims"]
        X = d["X"].reshape((nz, ny, nx))[:, 0, :]
        Z = d["Z"].reshape((nz, ny, nx))[:, 0, :]
        F = np.asarray(d[field], dtype=np.float64).reshape((nz, ny, nx))[:, 0, :]
        return X, Z, F, d["time"], d["pointDist"], files[-1]

    # synthetic fallback: GLL grid + y-average of the analytic field
    nx, ny, nz = DIMS
    t = 0.0
    x, _ = fe_read.axis_nodes_weights(nx, 1, X0[0], X1[0])
    y, _ = fe_read.axis_nodes_weights(ny, 1, X0[1], X1[1])
    z, _ = fe_read.axis_nodes_weights(nz, 1, X0[2], X1[2])
    Z3, Y3, X3 = np.meshgrid(z, y, x, indexing="ij")  # (nz, ny, nx)
    F3 = fe_read.ref_vz(X3, Y3, Z3, t) if field == "vz" else fe_read.ref_T(X3, Y3, Z3)
    F = fe_read.grid_avg(F3.ravel(), (nx, ny, nz), "y", (1, 1, 1))[:, 0, :]
    Z2, X2 = np.meshgrid(z, x, indexing="ij")
    return X2, Z2, F, t, (1, 1, 1), "(synthetic)"


X, Z, F, time, dist, src = load_plane(FIELD)

fig, ax = plt.subplots(figsize=(10, 4.2))

# flat shading with X/Z as cell corners: every mesh line IS a node line, so the
# GLL clustering toward the box edges shows up directly in the cell sizes
pc = ax.pcolormesh(Z, X, F[:-1, :-1], shading="flat", cmap="jet", edgecolors="k", linewidth=0.4)
ax.scatter(Z, X, s=6, c="k", zorder=3, label="sample points")

cb = fig.colorbar(pc, ax=ax)
cb.set_label(f"<{FIELD}>_y")
ax.set_xlabel("Z")
ax.set_ylabel("X")
ax.set_title(
    f"{src}: y-avg of {FIELD} on the GLL x-z grid "
    f"(pointDist={tuple(dist)}, 1=GLL), time={time:g}"
)
ax.legend(loc="upper right", fontsize=8)

fig.tight_layout()
fig.savefig("test_gll_viz.png", dpi=150)
print("wrote test_gll_viz.png")
if matplotlib.get_backend().lower() != "agg":
    plt.show()
