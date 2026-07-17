"""Validate the doAvg outputs (Plans C + D).

For each avg family, pair the latest _avg .vts with the same-step source box .vts
and compare against the numpy weighted average of the source data (exact
contract, tight tol; weights follow the source file's pointDist -- trapezoid for
uniform axes, GLL quadrature for GLL axes). Additionally cross-check against the
analytic reference fields averaged on the source grid (interp accuracy, loose tol).
"""
import glob
import sys

import numpy as np

import fe_read

# (avg glob, source glob, axes, mode)
CASES = [
    ("box3d_avgxy_g*.vts", "box3d[0-9]*.vts", "xy", "gather"),
    ("box3d_avgz_gs*.vts", "box3d[0-9]*.vts", "z", "gather-scatter"),
    ("box3d_avgxyz_g*.vts", "box3d[0-9]*.vts", "xyz", "gather"),
    # GLL box (Plan D): same reductions, GLL quadrature weights
    ("box3d_gll_avgxy_g*.vts", "box3d_gll[0-9]*.vts", "xy", "gather"),
    ("box3d_gll_avgz_gs*.vts", "box3d_gll[0-9]*.vts", "z", "gather-scatter"),
    ("box3d_gll_avgy_g*.vts", "box3d_gll[0-9]*.vts", "y", "gather"),
]


def check_analytic(avg_fname, src_fname, axes, mode, tol=1e-2):
    """Compare the avg file against the weighted average of the ANALYTIC fields
    evaluated on the source grid (validates interpolation + averaging end to end)."""
    da = fe_read.load(avg_fname)
    ds = fe_read.load(src_fname)
    nx, ny, nz = ds["dims"]
    X, Y, Z, t = ds["X"], ds["Y"], ds["Z"], da["time"]
    ok = True
    errs = []
    for key, ref_fld in (("vz", fe_read.ref_vz(X, Y, Z, t)), ("temp", fe_read.ref_T(X, Y, Z))):
        ref = fe_read.grid_avg(ref_fld, ds["dims"], axes, ds["pointDist"])
        if mode == "gather-scatter":
            ref = np.broadcast_to(ref, (nz, ny, nx))
        got = np.asarray(da[key], dtype=np.float64).reshape(ref.shape)
        # normalize by the UNaveraged field scale: averages can be ~0 (e.g. the
        # x-average of temp ~ sin(pi x) vanishes by symmetry)
        e = np.abs(got - ref).max() / np.abs(ref_fld).max()
        errs.append(e)
        ok = ok and (e < tol)
    print(
        f"{avg_fname}: vs analytic avg={axes} "
        f"relerr(vz)={errs[0]:.2e} relerr(temp)={errs[1]:.2e} -> {'PASS' if ok else 'FAIL'}"
    )
    return ok


ok = True
for avg_glob, src_glob, axes, mode in CASES:
    avgs = sorted(glob.glob(avg_glob))
    srcs = sorted(glob.glob(src_glob))
    if not avgs or not srcs:
        print(f"{avg_glob}: no files found -> FAIL")
        ok = False
        continue
    # process() and doAvg() run on the same cadence, so counters are aligned
    ok = fe_read.check_avg(avgs[-1], srcs[-1], axes, mode) and ok
    ok = check_analytic(avgs[-1], srcs[-1], axes, mode) and ok

sys.exit(0 if ok else 1)
