"""Validate the doAvg outputs (Plan C).

For each avg family, pair the latest _avg .vts with the same-step source box .vts
and compare against the numpy trapezoidal average of the source data (exact
contract, tight tol). Additionally cross-check against the analytic reference
fields trapezoid-averaged on the source grid (interp accuracy, loose tol).
"""
import glob
import sys

import numpy as np

import fe_read

# (avg glob, source glob, axes, mode)
CASES = [
    ("box3d_avgxyLine*.vts", "box3dBox*.vts", "xy", "gather"),
    ("box3d_avgzBox*.vts", "box3dBox*.vts", "z", "gather-scatter"),
    ("box3d_avgxyzPoint*.vts", "box3dBox*.vts", "xyz", "gather"),
    ("plane_avgyLine*.vts", "planePlane*.vts", "y", "gather"),
]


def check_analytic(avg_fname, src_fname, axes, mode, tol=1e-2):
    """Compare the avg file against the trapezoid-average of the ANALYTIC fields
    evaluated on the source grid (validates interpolation + averaging end to end)."""
    da = fe_read.load(avg_fname)
    ds = fe_read.load(src_fname)
    nx, ny, nz = ds["dims"]
    X, Y, Z, t = ds["X"], ds["Y"], ds["Z"], da["time"]
    ok = True
    errs = []
    for key, ref_fld in (("vz", fe_read.ref_vz(X, Y, Z, t)), ("temp", fe_read.ref_T(X, Y, Z))):
        ref = fe_read.trap_avg(ref_fld, ds["dims"], axes)
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
