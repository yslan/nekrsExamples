"""Validate the ADIOS2 .bp backend (Plan F).

turbPipe.udf emits the box3d sampler in BOTH formats on the same cadence:
  process()/doAvg(..., "adios")  -> box3d.bp, box3d_avg*_*.bp
  process()/doAvg(..., "vts")    -> box3d0000N.vts, box3d_avg*_*0000N.vts

Three checks:
  1. analytic  -- fe_read.check on box3d.bp (field vs analytic ref, tol 1e-2);
  2. avg       -- each _avg*.bp vs the weighted average of box3d.bp (exact
                  contract of doAvg, tight tol) + analytic cross-check;
  3. A/B       -- box3d.bp fields/coords match the paired last-step box3d*.vts to
                  float32 (both backends see the same interpolated data).

Run after `nrsmpi turbPipe <np>`; requires the adios2 python bindings.
"""
import glob
import sys

import numpy as np

import fe_read

# (avg .bp, source .bp, axes, mode) -- BP avg families vs the BP source box
AVG_CASES = [
    ("box3d_avgxy_g.bp", "box3d.bp", "xy", "gather"),
    ("box3d_avgz_gs.bp", "box3d.bp", "z", "gather-scatter"),
    ("box3d_avgxyz_g.bp", "box3d.bp", "xyz", "gather"),
]


def check_ab(bp_fname, vts_glob, tol=1e-6):
    """Assert the .bp fields/coords equal the paired last-step .vts to float32."""
    vts = sorted(glob.glob(vts_glob))
    if not vts:
        print(f"A/B {bp_fname}: no {vts_glob} to pair -> FAIL")
        return False
    db = fe_read.load(bp_fname)
    dv = fe_read.load(vts[-1])
    ok = True
    worst = []
    for key in ("X", "Y", "Z", "vz", "temp"):
        a = np.asarray(db[key], dtype=np.float64).ravel()
        b = np.asarray(dv[key], dtype=np.float64).ravel()
        e = np.abs(a - b).max() if a.shape == b.shape else np.inf
        worst.append(f"{key}={e:.1e}")
        ok = ok and (e < tol)
    print(f"A/B {bp_fname} vs {vts[-1]}: max|bp-vts| {' '.join(worst)} -> {'PASS' if ok else 'FAIL'}")
    return ok


ok = True

# 1. analytic field check on the BP box
ok = fe_read.check("box3d.bp") and ok

# 2. BP averaging families (exact contract + analytic), against the BP source
for avg_bp, src_bp, axes, mode in AVG_CASES:
    if not glob.glob(avg_bp):
        print(f"{avg_bp}: not found -> FAIL")
        ok = False
        continue
    ok = fe_read.check_avg(avg_bp, src_bp, axes, mode) and ok

# 3. cross-format A/B: BP must match the paired VTS to float32
ok = check_ab("box3d.bp", "box3d[0-9]*.vts") and ok

sys.exit(0 if ok else 1)
