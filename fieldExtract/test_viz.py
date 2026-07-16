"""Visualize all six fieldExtract samplers in one figure.

Usage:  python test_viz.py [field]      # field defaults to "vz" (also try "temp")

Loads the latest .vts for each config (box1d/box2d/box3d/line/plane/cyl) via fe_read,
and picks a view by effective dimensionality:
  1D -> line plot,  2D -> pcolormesh,  3D -> 3D scatter.
Writes test_viz_<field>.png and shows it if a display is available.
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

CONFIGS = [
    ("box1d", "box1dBox*.vts"),
    ("box2d", "box2dBox*.vts"),
    ("box3d", "box3dBox*.vts"),
    ("line", "lineLine*.vts"),
    ("plane", "planePlane*.vts"),
    ("cyl", "cylBox*.vts"),
]

FIELD = sys.argv[1] if len(sys.argv) > 1 else "vz"


def latest(pattern):
    files = sorted(glob.glob(pattern))
    return files[-1] if files else None


def ref_err(d, field):
    if field == "vz":
        ref = fe_read.ref_vz(d["X"], d["Y"], d["Z"], d["time"])
    elif field == "temp":
        ref = fe_read.ref_T(d["X"], d["Y"], d["Z"])
    else:
        return float("nan")
    return fe_read.relerr(d[field], ref)


fig = plt.figure(figsize=(15, 9))
labels = ["x", "y", "z"]

for i, (name, pattern) in enumerate(CONFIGS, start=1):
    fname = latest(pattern)
    if fname is None:
        ax = fig.add_subplot(2, 3, i)
        ax.text(0.5, 0.5, f"{name}\n(no {pattern})", ha="center", va="center")
        ax.set_axis_off()
        continue

    d = fe_read.load(fname)
    if FIELD not in d:
        ax = fig.add_subplot(2, 3, i)
        ax.text(0.5, 0.5, f"{name}\n(no field '{FIELD}')", ha="center", va="center")
        ax.set_axis_off()
        continue

    nx, ny, nz = d["dims"]
    shape = (nz, ny, nx)
    coords = [np.squeeze(d[c].reshape(shape)) for c in ("X", "Y", "Z")]
    Fs = np.squeeze(d[FIELD].reshape(shape))
    err = ref_err(d, FIELD)

    if Fs.ndim <= 1:  # 1D -> line
        ax = fig.add_subplot(2, 3, i)
        k = int(np.argmax([np.ptp(c) for c in coords]))
        ax.plot(coords[k].ravel(), Fs.ravel(), "-o", ms=3)
        ax.set_xlabel(labels[k])
        ax.set_ylabel(FIELD)
    elif Fs.ndim == 2:  # 2D -> colormap
        ax = fig.add_subplot(2, 3, i)
        varying = [k for k in range(3) if np.ptp(coords[k]) > 1e-12]
        a, b = varying[0], varying[1]
        pc = ax.pcolormesh(coords[a], coords[b], Fs, shading="auto", cmap="jet")
        fig.colorbar(pc, ax=ax)
        ax.set_xlabel(labels[a])
        ax.set_ylabel(labels[b])
        ax.set_aspect("equal")
    else:  # 3D -> scatter
        ax = fig.add_subplot(2, 3, i, projection="3d")
        p = ax.scatter(d["X"], d["Y"], d["Z"], c=d[FIELD], cmap="jet", s=8)
        fig.colorbar(p, ax=ax, shrink=0.6)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

    ax.set_title(f"{name} {d['dims']}\nrelerr={err:.1e}")

fig.suptitle(f"fieldExtract samplers — {FIELD}")
fig.tight_layout()
out = f"test_viz_{FIELD}.png"
fig.savefig(out, dpi=120)
print(f"wrote {out}")
plt.show()
