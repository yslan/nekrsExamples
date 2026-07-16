"""Reference reader for fieldExtract .vts output.

Loads a StructuredGrid .vts, reshapes fields using the boxDims stored in
FieldData, and validates a 2D (ny=1) slice against the analytic fields used in
turbPipe.udf::test_sol.
"""
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib.pyplot as plt

fname = "tttBox00003.vts"

reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName(fname)
reader.Update()

grid = reader.GetOutput()
points = vtk_to_numpy(grid.GetPoints().GetData())
print("dims:", grid.GetDimensions())
print("npts:", grid.GetNumberOfPoints())
print("points:", points.shape)

print("FieldData")
fd = grid.GetFieldData()
for i in range(fd.GetNumberOfArrays()):
    arr = fd.GetArray(i)
    print(" " * 4, arr.GetName(), vtk_to_numpy(arr))

print("PointData")
pd = grid.GetPointData()
for i in range(pd.GetNumberOfArrays()):
    arr = pd.GetArray(i)
    print(" " * 4, arr.GetName(), vtk_to_numpy(arr).shape)

time = vtk_to_numpy(fd.GetArray("TimeValue")).ravel().item()
nx, ny, nz = vtk_to_numpy(fd.GetArray("boxDims")).ravel().astype(int)
print("time:", time)
print("boxDims:", nx, ny, nz)

# reshape following the (nz, ny, nx) lexicographic ordering
X = points[:, 0].reshape((nz, ny, nx))
Z = points[:, 2].reshape((nz, ny, nx))
W = vtk_to_numpy(pd.GetArray("vz")).reshape((nz, ny, nx))
T = vtk_to_numpy(pd.GetArray("temp")).reshape((nz, ny, nx))

# ny = 1: take the only y-plane
X = X[:, 0, :]
Z = Z[:, 0, :]
W = W[:, 0, :]
T = T[:, 0, :]
print("slice shape:", X.shape)


# analytic references (see turbPipe.udf::test_sol)
def chk_W(X, W, time):
    R = 0.5
    rr = np.sqrt((X / R) ** 2)
    Wref = 6 / 5.0 * (1 - np.power(rr, 6)) * np.cos(time)
    Wref[rr > 1.0] = 0
    return abs(W - Wref).max() / abs(Wref).max()


def chk_T(X, Z, T):
    R = 0.5
    rr = np.sqrt((X / R) ** 2)
    Tref = np.sin(np.pi * X) * np.cos(2.7 * np.pi * Z + 0.2)
    Tref[rr > 1.0] = 0
    return abs(T - Tref).max() / abs(Tref).max()


# plot
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
err = chk_W(X, W, time)
ax.plot_surface(Z, X, W, edgecolor="k", linewidth=0.3, antialiased=True, cmap="jet")
ax.set_xlabel("X")
ax.set_ylabel("Z")
ax.set_zlabel("W")
ax.set_title(f"Time={time:g}, relerr={err:.2e}")

fig, ax = plt.subplots()
err = chk_T(X, Z, T)
pc = ax.pcolormesh(Z, X, T, cmap="jet", edgecolor="k", linewidth=0.0001)
cb = fig.colorbar(pc, ax=ax)
ax.set_aspect("equal")
ax.set_xlabel("Z")
ax.set_ylabel("X")
cb.set_label("T")
ax.set_title(f"Time={time:g}, relerr={err:.2e}")

plt.show()
