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


def trap_avg(F, dims, axes):
    """Trapezoidal average of flat field F (dims=(nx,ny,nz)) over axes in 'xyz'.

    Mirrors fieldExtract::doAvg: endpoints half weight, normalized by (n-1) per
    averaged axis. Returns shape (nz', ny', nx') with averaged axes kept as 1.
    """
    nx, ny, nz = dims
    A = np.asarray(F, dtype=np.float64).reshape((nz, ny, nx))
    for ax in axes.lower():
        axis = {"x": 2, "y": 1, "z": 0}[ax]
        n = A.shape[axis]
        w = np.ones(n)
        w[0] = w[-1] = 0.5
        shape = [1, 1, 1]
        shape[axis] = n
        A = (A * w.reshape(shape)).sum(axis=axis, keepdims=True) / (n - 1)
    return A


def check_avg(avg_fname, src_fname, axes, mode, tol=1e-5):
    """Compare an _avg .vts against the numpy trapezoid-average of its source box .vts.

    Exercises the exact contract of fieldExtract::doAvg (same grid, same weights),
    so tol is tight -- limited only by float32 storage / summation order.
    """
    da = load(avg_fname)
    ds = load(src_fname)
    nx, ny, nz = ds["dims"]
    ok = True
    errs = []
    for key in ("vz", "temp"):
        ref = trap_avg(ds[key], ds["dims"], axes)
        if mode == "gather-scatter":
            ref = np.broadcast_to(ref, (nz, ny, nx))
        got = np.asarray(da[key], dtype=np.float64).reshape(ref.shape)
        # normalize by the UNaveraged source-field scale: averages can be ~0
        # (odd symmetry), while the float32 storage noise is absolute
        e = np.abs(got - ref).max() / np.abs(ds[key]).max()
        errs.append(e)
        ok = ok and (e < tol)
    print(
        f"{avg_fname}: dims={da['dims']} avg={axes} mode={mode} "
        f"relerr(vz)={errs[0]:.2e} relerr(temp)={errs[1]:.2e} -> {'PASS' if ok else 'FAIL'}"
    )
    return ok
