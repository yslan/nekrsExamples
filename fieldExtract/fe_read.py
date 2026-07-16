"""Shared loader/validator for fieldExtract .vts output.

Analytic reference fields mirror turbPipe.udf::test_sol:
  vz   = 6/5 (1 - (r/R)^6) cos(t)          (r = sqrt(x^2+y^2), R = 0.5)
  temp = sin(pi x) cos(2.7 pi z + 0.2)

Each per-config test script (test_box1d.py, ...) imports check() from here.
"""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

R = 0.5  # pipe radius


def load(fname):
    """Return dict with dims, time, X/Y/Z arrays, and one entry per PointData field."""
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()
    grid = reader.GetOutput()

    pts = vtk_to_numpy(grid.GetPoints().GetData())
    fd = grid.GetFieldData()
    pd = grid.GetPointData()

    nx, ny, nz = vtk_to_numpy(fd.GetArray("boxDims")).ravel().astype(int)
    out = {
        "dims": (int(nx), int(ny), int(nz)),
        "time": vtk_to_numpy(fd.GetArray("TimeValue")).ravel().item(),
        "X": pts[:, 0],
        "Y": pts[:, 1],
        "Z": pts[:, 2],
    }
    for i in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(i)
        out[arr.GetName()] = vtk_to_numpy(arr)
    return out


def ref_vz(X, Y, Z, time):
    rr = np.sqrt(X * X + Y * Y) / R
    W = 6 / 5.0 * (1 - rr**6) * np.cos(time)
    W[rr > 1.0] = 0.0
    return W


def ref_T(X, Y, Z):
    rr = np.sqrt(X * X + Y * Y) / R
    T = np.sin(np.pi * X) * np.cos(2.7 * np.pi * Z + 0.2)
    T[rr > 1.0] = 0.0
    return T


def relerr(a, b):
    denom = np.abs(b).max()
    return np.abs(a - b).max() / denom if denom > 0 else np.abs(a - b).max()


def check(fname, tol=1e-2):
    """Load fname, compare vz/temp to the analytic refs, print + return PASS/FAIL."""
    d = load(fname)
    X, Y, Z, t = d["X"], d["Y"], d["Z"], d["time"]
    ev = relerr(d["vz"], ref_vz(X, Y, Z, t))
    eT = relerr(d["temp"], ref_T(X, Y, Z))
    ok = (ev < tol) and (eT < tol)
    print(
        f"{fname}: dims={d['dims']} time={t:g} "
        f"relerr(vz)={ev:.2e} relerr(temp)={eT:.2e} -> {'PASS' if ok else 'FAIL'}"
    )
    return ok
