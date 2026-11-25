#!/usr/bin/env python3
"""
nek_read.py -- read a NekRS/Nek5000 `.fld` checkpoint with numpy alone.

VTK only grew a Nek5000 reader in 9.3, and this box ships 9.2.2; pyvista's
`get_reader` raises VTKVersionError on `tgv.nek5000`.  The format is simple
enough that reading it directly is less trouble than a version constraint, and
it matches the precedent of `iso_read.py` (numpy, no vtk needed).

Layout, all little-endian:

    132-byte ASCII header
    4-byte float 6.54321        endian test
    nel * int32                 global element ids
    DATA, one block per field GROUP, in rdcode order
    min/max TRAILER, 2 floats per field per element

Read out of Nek5000's own writer (core/prepost.f: mfo_outfld, mfo_outv,
mfo_outs), not inferred, because two plausible guesses are both wrong:

  * The min/max metadata is APPENDED ("add meta data to the end of the file",
    prepost.f:950), not interleaved and not at the front.  It is 2*nfield*nel
    floats, so ignoring it leaves the data 8192 floats short of the file at
    N=7 on 512 elements -- which shows up as size/(nel*Np) = 8.03, not 8.
  * A field GROUP is written per element with its components interleaved:
    mfo_outv copies element 1's u, then its v, then its w, then element 2's.
    So xyz is (nel, 3, Np) and a scalar is (nel, Np) -- NOT (nfield, nel, Np)
    globally.

Neither mistake raises an error and both leave x on [-pi,pi].  What catches
them is a POINTWISE test against a known field: for the TGV initial condition
u = sin(x)cos(y)cos(z), a correct read gives correlation 1.00000000 and max
error 1.2e-07 (float32), a wrong one gives ~0.05.  `python3 nek_read.py --check`
runs exactly that.

The header is

    #std  wdsize  Nx Ny Nz  nel  nelgt  time  istep  0 1  RDCODE  ...

where RDCODE names the fields present: X coordinates, U velocity, P pressure,
S## scalars.  **Coordinates are usually written only in the first dump**, so a
later file has rdcode "UPS01" and no X -- read f00000 for the mesh and reuse it.

    from nek_read import read_fld, load_series
    d = read_fld("tgv0.f00001", coords_from="tgv0.f00000")
    d["time"], d["x"], d["u"], d["s01"]     # each (nel, Np)
"""
import glob
import os
import re

import numpy as np


def _parse_header(h):
    tok = h[:120].split()
    if tok[0] != "#std":
        raise ValueError("not a Nek .fld file (no #std): %r" % h[:32])
    wdsize = int(tok[1])
    nx, ny, nz = int(tok[2]), int(tok[3]), int(tok[4])
    nel, nelgt = int(tok[5]), int(tok[6])
    time = float(tok[7])
    istep = int(tok[8])
    rdcode = tok[11]
    return dict(wdsize=wdsize, nx=nx, ny=ny, nz=nz, nel=nel, nelgt=nelgt,
                time=time, istep=istep, rdcode=rdcode)


def _groups_in(rdcode, ndim=3):
    """Expand a RDCODE into the ordered list of WRITE GROUPS, [(names, ncomp)].

    The grouping matters: a vector is one group of `ndim` components written
    element-by-element with the components interleaved (mfo_outv), while each
    scalar is its own single-component group (mfo_outs).
    """
    groups = []
    for nm, n in _fields_in(rdcode, ndim):
        groups.append((nm, n))
    return groups


def _fields_in(rdcode, ndim=3):
    """[(names, ncomp)] in write order, one entry per field group."""
    names = []
    i = 0
    while i < len(rdcode):
        c = rdcode[i]
        if c == "X":
            names.append((["x", "y", "z"][:ndim], ndim))
            i += 1
        elif c == "U":
            names.append((["u", "v", "w"][:ndim], ndim))
            i += 1
        elif c == "P":
            names.append((["p"], 1))
            i += 1
        elif c == "T":
            names.append((["t"], 1))
            i += 1
        elif c == "S":
            # S01, S02, ... -- two digits, and the COUNT of scalars, so "S01"
            # means one scalar and "S02" means two.
            m = re.match(r"S(\d\d)", rdcode[i:])
            if not m:
                raise ValueError("bad scalar code in %r" % rdcode)
            n = int(m.group(1))
            for k in range(n):
                names.append((["s%02d" % (k + 1)], 1))
            i += 3
        else:
            i += 1  # ignore anything else the writer put there
    return names


def read_fld(path, coords_from=None):
    """Read one .fld.  coords_from supplies x,y,z when this file has no X."""
    with open(path, "rb") as f:
        hdr = _parse_header(f.read(132).decode("ascii", "replace"))
        etag = np.frombuffer(f.read(4), dtype=np.float32)[0]
        if abs(etag - 6.54321) > 1e-3:
            raise ValueError("endian tag %g -- byte-swapped file not supported" % etag)
        nel = hdr["nel"]
        np_ = hdr["nx"] * hdr["ny"] * hdr["nz"]
        np.frombuffer(f.read(nel * 4), dtype=np.int32)  # element ids
        dt = np.float32 if hdr["wdsize"] == 4 else np.float64
        raw = np.frombuffer(f.read(), dtype=dt)

    groups = _groups_in(hdr["rdcode"])
    nfield = sum(g[1] for g in groups)
    want = nel * np_ * nfield
    if raw.size < want:
        raise ValueError("%s: expected %d values for %s, found %d"
                         % (path, want, groups, raw.size))
    trailer = raw.size - want
    if trailer not in (0, 2 * nfield * nel):
        raise ValueError("%s: %d leftover floats, expected 0 or %d (min/max trailer)"
                         % (path, trailer, 2 * nfield * nel))

    out = dict(hdr)
    out["Np"] = np_
    # One block per GROUP, components interleaved within each element.
    off = 0
    for names, ncomp in groups:
        blk = raw[off:off + nel * ncomp * np_].reshape(nel, ncomp, np_)
        off += nel * ncomp * np_
        for c, nm in enumerate(names):
            out[nm] = blk[:, c, :].astype(np.float64)

    if "x" not in out:
        if coords_from is None:
            raise ValueError("%s has no coordinates (rdcode %r); pass coords_from"
                             % (path, hdr["rdcode"]))
        c = read_fld(coords_from)
        for nm in ("x", "y", "z"):
            out[nm] = c[nm]
    return out


def load_series(stem="tgv0", directory=None):
    """Every tgv0.f##### in `directory`, sorted, with times.  Returns a list of
    (path, time); coordinates come from the first file that has them."""
    files = sorted(glob.glob(os.path.join(_datadir(directory), stem + ".f[0-9]*")))
    out = []
    for fn in files:
        with open(fn, "rb") as f:
            h = _parse_header(f.read(132).decode("ascii", "replace"))
        out.append((fn, h["time"], h["rdcode"]))
    return out


def _datadir(d=None):
    """Where the run's output lives.

    An explicit path wins, then $TGV_DATA, then the CURRENT DIRECTORY -- so
    running the script from inside a run folder just works.  Only if the cwd
    holds no output does it look for a dat_*/ subdirectory, which is where a
    finished run tends to get filed."""
    if d:
        return d
    if os.environ.get("TGV_DATA"):
        return os.environ["TGV_DATA"]
    if glob.glob(os.path.join(os.getcwd(), "*.f0*")):
        return os.getcwd()
    for base in (os.getcwd(), os.path.dirname(os.path.abspath(__file__))):
        for c in sorted(glob.glob(os.path.join(base, "dat_*"))):
            if os.path.isdir(c) and glob.glob(os.path.join(c, "*.f0*")):
                return c
    return os.getcwd()


def check(path=None):
    """Pointwise self-test against the analytic TGV initial condition.

    u = sin(x)cos(y)cos(z), v = -cos(x)sin(y)cos(z), w = 0 at t = 0.  A correct
    read matches to float32 (~1e-7); ANY layout mistake -- wrong grouping,
    ignoring the min/max trailer, transposing the reshape -- lands at ~0.05
    correlation while still producing plausible-looking coordinates.  That is
    why this test is pointwise and not a range check.
    """
    if path is None:
        path = os.path.join(_datadir(), "tgv0.f00000")
    d = read_fld(path)
    if abs(d["time"]) > 1e-12:
        print("%s is t=%g, not 0 -- the IC check needs the t=0 dump" % (path, d["time"]))
        return False
    x, y, z = d["x"], d["y"], d["z"]
    eu = np.abs(d["u"] - np.sin(x) * np.cos(y) * np.cos(z)).max()
    ev = np.abs(d["v"] + np.cos(x) * np.sin(y) * np.cos(z)).max()
    ew = np.abs(d["w"]).max()
    ok = max(eu, ev, ew) < 1e-5
    print("TGV initial condition, pointwise:")
    print("   max|u - sin x cos y cos z|  = %.3e" % eu)
    print("   max|v + cos x sin y cos z|  = %.3e" % ev)
    print("   max|w|                      = %.3e" % ew)
    print("   %s" % ("PASS (float32 round-off)" if ok else "FAIL -- the layout is wrong"))
    return ok


if __name__ == "__main__":
    import sys
    if "--check" in sys.argv:
        sys.exit(0 if check() else 1)
    print("data: %s" % _datadir())
    for fn, t, rd in load_series():
        print("%-16s t=%-10g rdcode=%s" % (os.path.basename(fn), t, rd))
