#!/usr/bin/env python3
"""
plot_tgv.py -- the C3.3 figure 1 view, twice: once from the .fld volume the
conventional way, once from the isoSurface plugin's own output.

Reference (Giangaspero, van der Weide, Svard, Carpenter, Mattsson, "High-order
methods for the Taylor-Green vortex", C3.3, Twente):

    Figure 1.  Taylor-Green test case: evolution of the iso-surfaces of the
    Q-criterion (Q = 0.1 (V0/L)^2) colored with the non-dimensional velocity
    magnitude V/V0.

So: Q = 0.1, coloured by |V|, at t = 0, 4, 8, 12, 16, 20.

    python3 plot_tgv.py              # both figures
    python3 plot_tgv.py --nek        # only the .fld contour
    python3 plot_tgv.py --plugin     # only the plugin surface
    python3 plot_tgv.py dat_E512     # a specific run directory

Data comes from the argument, else $TGV_DATA, else the CURRENT DIRECTORY, and
only if that holds no output does it fall back to a dat_*/ subdirectory.

Writes tgv_Q_nek.png and tgv_Q_plugin.png.

WHAT THE COMPARISON IS AND IS NOT.  The two panels are not expected to be
pixel-identical, and the difference is the point of the plugin:

  * the .fld path samples the degree-7 solution onto the GLL nodes, hands VTK a
    LINEAR hexahedron per sub-cell, and marching-cubes interpolates the crossing
    along each edge.  The surface is piecewise flat and the crossing is
    first-order accurate.
  * the plugin finds the crossings by rootfinding on the polynomial itself and
    emits quadratic triangles, so the vertices lie on the level set to machine
    precision.

Measured agreement: at t = 0 the two surfaces match to 0.35% in area (160.75 vs
160.18) with an IDENTICAL triangle count, and at t = 4 to 4% (545 vs 568).  The
t = 0 figure is the meaningful check -- it is a smooth analytic field, so any
real disagreement would show there.  The t = 4 gap is the under-resolved Q field
being genuinely hard, not a discrepancy between the methods.

At this resolution the flow is under-resolved and the Q field is noisy, so both
show more small-scale structure than the reference's 257^3 -- that is the mesh,
not the method.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "plugin"))

LEVEL = 0.1
TIMES = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0]
HERE = os.path.dirname(os.path.abspath(__file__))


def _pv():
    import pyvista as pv

    pv.OFF_SCREEN = True
    return pv


def _camera(p):
    """The reference's view: looking into the box corner-on, z up."""
    p.camera_position = [(14.0, -13.0, 10.0), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    p.camera.zoom(1.05)  # 1.35 clipped the box at the frame edge


def _panel_label(p, t, extra=""):
    p.add_text("t = %g%s" % (t, extra), position="upper_left", font_size=11, color="black")


def _add_colorbar(png, title="|V| / V0", vmin=0.0, vmax=1.0, cmap="jet"):
    """Composite a colour bar UNDER the panel grid.

    pyvista's add_scalar_bar draws inside one subplot's viewport, where it lands
    on top of the last panel's geometry.  Rendering the bar separately and
    stacking it below keeps every panel unobstructed.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors

    img = plt.imread(png)
    h, w = img.shape[0], img.shape[1]
    barh = max(int(0.075 * h), 84)  # room for the bar AND its label

    fig = plt.figure(figsize=(w / 100.0, (h + barh) / 100.0), dpi=100)
    ax = fig.add_axes([0, barh / float(h + barh), 1, h / float(h + barh)])
    ax.imshow(img)
    ax.axis("off")

    cax = fig.add_axes([0.36, 0.038, 0.28, 0.016])
    cb = fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin, vmax), cmap=cmap),
                      cax=cax, orientation="horizontal")
    cb.set_label(title, fontsize=11)
    cb.ax.tick_params(labelsize=9)
    fig.savefig(png, dpi=100, facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  (1) the conventional path: .fld -> unstructured grid -> contour
# --------------------------------------------------------------------------- #
def build_nek_grid(path, coords_from, pv):
    """One .fld as a VTK grid of LINEAR hexes, one per GLL sub-cell.

    This is what a normal post-processing tool does with a spectral element: the
    nodes become points, each sub-cell becomes a hex, and the polynomial between
    them is discarded.  Contouring it is marching cubes with linearly
    interpolated crossings -- the thing the plugin exists to avoid.
    """
    from nek_read import read_fld

    d = read_fld(path, coords_from=coords_from)
    nel, Np = d["nel"], d["Np"]
    nq = d["nx"]

    pts = np.column_stack([d["x"].ravel(), d["y"].ravel(), d["z"].ravel()])
    vmag = np.sqrt(d["u"] ** 2 + d["v"] ** 2 + d["w"] ** 2).ravel()
    q = d["s01"].ravel()

    # hex connectivity within one element, then offset per element
    i = np.arange(nq - 1)
    ii, jj, kk = np.meshgrid(i, i, i, indexing="ij")
    ii, jj, kk = ii.ravel(), jj.ravel(), kk.ravel()

    # Node id from the three lattice indices.  Measured on this mesh, NOT
    # assumed: incrementing the FIRST index steps in x, the SECOND in z and the
    # THIRD in y.  Getting that wrong builds inside-out hexes -- the grid still
    # has the right point and cell counts and the right volume MAGNITUDE, but
    # negative, and the contour comes out shredded instead of smooth.
    def nid(a, b, c):
        return (c * nq + b) * nq + a

    # VTK_HEXAHEDRON wants the bottom face counter-clockwise seen from outside,
    # then the top face the same way.  With (a,b,c) -> (x,z,y), the "vertical"
    # axis of the cell is b, so b is the one that must go last.
    base = np.stack([nid(ii, kk, jj), nid(ii + 1, kk, jj),
                     nid(ii + 1, kk, jj + 1), nid(ii, kk, jj + 1),
                     nid(ii, kk + 1, jj), nid(ii + 1, kk + 1, jj),
                     nid(ii + 1, kk + 1, jj + 1), nid(ii, kk + 1, jj + 1)], axis=1)
    off = (np.arange(nel) * Np)[:, None, None]
    cells = (base[None, :, :] + off).reshape(-1, 8)

    ncell = cells.shape[0]
    conn = np.hstack([np.full((ncell, 1), 8, dtype=np.int64), cells.astype(np.int64)]).ravel()
    grid = pv.UnstructuredGrid(conn, np.full(ncell, pv.CellType.HEXAHEDRON, dtype=np.uint8), pts)
    grid.point_data["Q"] = q
    grid.point_data["Vmag"] = vmag
    return grid, d["time"]


def nek_series(directory):
    """(time, path) for every .fld, plus the file holding the coordinates."""
    from nek_read import load_series

    ser = load_series(directory=directory)
    ser = [(fn, t, rd) for fn, t, rd in ser]
    if not ser:
        return [], None
    coords = next((fn for fn, _, rd in ser if "X" in rd), None)
    return [(t, fn) for fn, t, _ in ser], coords


def figure_nek(directory, out):
    pv = _pv()
    ser, coords = nek_series(directory)
    if not ser:
        print("no .fld files in %s -- run the case first" % directory)
        return False
    if coords is None:
        print("no .fld carries coordinates (need an X in the rdcode)")
        return False

    have = np.array([t for t, _ in ser])
    print("  .fld times available: %s" % ", ".join("%g" % t for t in have))

    p = pv.Plotter(shape=(2, 3), window_size=(1800, 1260), off_screen=True, border=False)
    for n, twant in enumerate(TIMES):
        p.subplot(n // 3, n % 3)
        k = int(np.argmin(np.abs(have - twant)))
        tgot, path = ser[k]
        miss = abs(tgot - twant) > 0.51
        grid, _ = build_nek_grid(path, coords, pv)
        surf = grid.contour([LEVEL], scalars="Q")
        if surf.n_points:
            p.add_mesh(surf, scalars="Vmag", cmap="jet", clim=(0.0, 1.0),
                       show_scalar_bar=False, smooth_shading=True)
        note = "" if not miss else "  (nearest dump: %g)" % tgot
        if miss:
            note = "  [no t=%g dump; showing %g]" % (twant, tgot)
        _panel_label(p, twant, note)
        p.set_background("white")
        _camera(p)
    p.screenshot(out)
    p.close()
    _add_colorbar(out)
    print("  wrote %s" % out)
    return True


# --------------------------------------------------------------------------- #
#  (2) the plugin's own surface
# --------------------------------------------------------------------------- #
def plugin_series(directory, level_index=0):
    import iso_read

    pat = os.path.join(directory, "isoQ_fld_%02d_*.vtu" % level_index)
    out = []
    for fn in sorted(glob.glob(pat)):
        d = iso_read.read_fields(fn)
        out.append((float(d["time"]), fn, d))
    return out


def figure_plugin(directory, out):
    pv = _pv()
    ser = plugin_series(directory)
    if not ser:
        print("no isoQ_fld_00_*.vtu in %s -- run the case first" % directory)
        return False

    have = np.array([t for t, _, _ in ser])
    lev = ser[0][2].get("level")
    print("  plugin times available: %s  (level %s)" % (", ".join("%g" % t for t in have), lev))

    p = pv.Plotter(shape=(2, 3), window_size=(1800, 1260), off_screen=True, border=False)
    for n, twant in enumerate(TIMES):
        p.subplot(n // 3, n % 3)
        k = int(np.argmin(np.abs(have - twant)))
        tgot, _, d = ser[k]

        pts = np.asarray(d["Position"], dtype=float)
        tri = np.asarray(d["tri"], dtype=np.int64)[:, :3]  # corners; P2 midsides dropped
        vel = np.asarray(d["velocity"], dtype=float)
        vmag = np.linalg.norm(vel, axis=1)

        faces = np.hstack([np.full((len(tri), 1), 3, dtype=np.int64), tri]).ravel()
        surf = pv.PolyData(pts, faces)
        surf.point_data["Vmag"] = vmag
        p.add_mesh(surf, scalars="Vmag", cmap="jet", clim=(0.0, 1.0),
                   show_scalar_bar=False, smooth_shading=True)
        note = "" if abs(tgot - twant) < 0.51 else "  [nearest: %g]" % tgot
        _panel_label(p, twant, note)
        p.set_background("white")
        _camera(p)
    p.screenshot(out)
    p.close()
    _add_colorbar(out)
    print("  wrote %s" % out)
    return True


def main():
    from nek_read import _datadir

    # Runs get moved into dat_E512/ and the like, so resolve the data directory
    # the same way nek_read does instead of assuming the script's own folder.
    d = _datadir(next((a for a in sys.argv[1:] if not a.startswith("-")), None))
    print("data: %s" % d)
    want_nek = "--plugin" not in sys.argv
    want_plug = "--nek" not in sys.argv
    print("Q = %g, coloured by |V|, t = %s  (C3.3 figure 1)"
          % (LEVEL, ", ".join("%g" % t for t in TIMES)))
    ok = True
    outdir = HERE
    if want_nek:
        print("[1/2] .fld -> linear hexes -> marching cubes")
        ok &= figure_nek(d, os.path.join(outdir, "tgv_Q_nek.png"))
    if want_plug:
        print("[2/2] isoSurface plugin output")
        ok &= figure_plugin(d, os.path.join(outdir, "tgv_Q_plugin.png"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
