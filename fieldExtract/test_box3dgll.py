"""Validate the box3dgll sampler (box mode, GLL points on all axes, {11,7,21}).

Checks 1) the .vts node positions match the zwgll-mapped GLL grid, 2) the
interpolated fields match the analytic references (fe_read.check).
"""
import glob
import sys

import numpy as np

import fe_read

X0 = (-0.3, -0.3, 2.0)  # must match turbPipe.udf::feBox3dGll
X1 = (0.3, 0.3, 4.0)

pattern = "box3dgllBox*.vts"
fname = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob(pattern))[-1]

d = fe_read.load(fname)
nx, ny, nz = d["dims"]

ok = d["pointDist"] == (1, 1, 1)
if not ok:
    print(f"{fname}: pointDist={d['pointDist']} (expected (1, 1, 1)) -> FAIL")

# node positions: each axis line must match the zwgll nodes mapped to [X0, X1]
extract = {
    "X": lambda A: A[0, 0, :],
    "Y": lambda A: A[0, :, 0],
    "Z": lambda A: A[:, 0, 0],
}
for ia, (key, n) in enumerate((("X", nx), ("Y", ny), ("Z", nz))):
    ref, _ = fe_read.axis_nodes_weights(n, d["pointDist"][ia], X0[ia], X1[ia])
    got = extract[key](d[key].reshape((nz, ny, nx)))
    e = np.abs(got - ref).max() / max(abs(X0[ia]), abs(X1[ia]))
    axis_ok = e < 1e-6  # float32 storage of the coordinates
    ok = ok and axis_ok
    print(f"{fname}: {key} nodes vs zwgll relerr={e:.2e} -> {'PASS' if axis_ok else 'FAIL'}")

# field values vs analytic references
ok = fe_read.check(fname) and ok

sys.exit(0 if ok else 1)
