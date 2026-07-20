"""Validate azimuthal (theta) averaging of the cylinder demo (Plan G).

turbPipe_cyl.udf samples a full-circle r-theta-z grid (points mode + setQuadrature)
and averages in theta through the plugin:
  cyldemo*.vts           raw Cartesian point cloud
  cyldemo_avgy_g*.vts    theta-collapsed (r,z) profile  (x=radius, y=0)
  cyldemo_avgy_gs*.vts   theta-average broadcast back onto the cloud (rings)

Two checks:
  1. exact contract -- the plugin's gather output equals the numpy uniform-periodic
     theta-mean of the same-step source cloud (float32 tol);
  2. analytic -- <vz>_theta = vz(r) (axisymmetric) and <temp>_theta = 0 (odd in x
     over the full circle), each normalized by the unaveraged field scale.

Run after `nrsmpi turbPipe_cyl <np>`.
"""
import glob
import sys

import numpy as np

import fe_read


def latest(pattern):
    fs = sorted(glob.glob(pattern))
    return fs[-1] if fs else None


ok = True

src = latest("cyldemo[0-9]*.vts")
avg = latest("cyldemo_avgy_g[0-9]*.vts")
if src is None or avg is None:
    print(f"missing inputs (src={src}, avg={avg}) -> FAIL")
    sys.exit(1)

ds = fe_read.load(src)
da = fe_read.load(avg)
nr, nt, nz = ds["dims"]

# 1. exact contract: plugin gather == numpy azimuthal mean of the source
errs = []
for k in ("vz", "temp"):
    ref = fe_read.azimuthal_avg(ds[k], ds["dims"]).reshape(nz, 1, nr)
    got = np.asarray(da[k], dtype=np.float64).reshape(nz, 1, nr)
    e = np.abs(got - ref).max() / np.abs(ds[k]).max()
    errs.append(e)
    ok = ok and (e < 1e-5)
print(f"{avg}: exact theta-mean vs source  relerr(vz)={errs[0]:.2e} relerr(temp)={errs[1]:.2e} "
      f"-> {'PASS' if errs[0] < 1e-5 and errs[1] < 1e-5 else 'FAIL'}")

# reduced-grid coordinate convention: x=radius, y=0
if np.abs(da["Y"]).max() > 1e-6:
    print(f"  WARN: expected y=0 on the averaged plane, got max|Y|={np.abs(da['Y']).max():.2e}")

# 2. analytic: <vz>_theta = vz(r,z), <temp>_theta = 0
X, Z, t = da["X"], da["Z"], da["time"]
vz_ref = fe_read.ref_vz(X, np.zeros_like(X), Z, t)
e_vz = np.abs(np.asarray(da["vz"]).ravel() - vz_ref).max() / np.abs(vz_ref).max()
e_temp = np.abs(np.asarray(da["temp"]).ravel()).max() / np.abs(ds["temp"]).max()
ok = ok and (e_vz < 1e-2) and (e_temp < 1e-2)
print(f"{avg}: analytic <vz>_theta=vz(r) relerr={e_vz:.2e}  <temp>_theta=0 relerr={e_temp:.2e} "
      f"-> {'PASS' if e_vz < 1e-2 and e_temp < 1e-2 else 'FAIL'}")

sys.exit(0 if ok else 1)
