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

print('FieldData')
fd = grid.GetFieldData()
for i in range(fd.GetNumberOfArrays()):
    arr = fd.GetArray(i)
    data = vtk_to_numpy(arr)
    print(' '*4, arr.GetName(), data)

print('PointData')
pd = grid.GetPointData()
for i in range(pd.GetNumberOfArrays()):
    arr = pd.GetArray(i)
    data = vtk_to_numpy(arr)
    print(' '*4, arr.GetName(), data.shape)

time = vtk_to_numpy(fd.GetArray('TimeValue')).ravel().item()
print(time)

X = points[:,0]
Y = points[:,1]
Z = points[:,2]

W = vtk_to_numpy(pd.GetArray("vz"))
T = vtk_to_numpy(pd.GetArray("temp"))


nx, ny, nz = vtk_to_numpy(fd.GetArray("boxDims")).ravel().astype(int)
print(nx,ny,nz)

X = X.reshape((nz,ny,nx))
Y = Y.reshape((nz,ny,nx))
Z = Z.reshape((nz,ny,nx))
W = W.reshape((nz,ny,nx))
T = T.reshape((nz,ny,nx))

# Since ny = 1, take the only y-plane
X = X[:, 0, :]
Y = Y[:, 0, :]
Z = Z[:, 0, :]
W = W[:, 0, :]
T = T[:, 0, :]
print(X.shape)


# chk
def chk_W(X,Y,Z,W,time):
    zLength = 6
    R = 0.5
    xr = X / R
    yr = Y / R
    rr = np.sqrt(xr*xr + yr*yr)
    zo = 2 * np.pi * Z / zLength

    Wref = np.zeros_like(X)
    Wref = 6/5. * (1 - np.power(rr, 6)) * np.cos(time)
    Wref[rr>1.0] = 0
    relerr = abs(W-Wref).max() / abs(Wref).max()
    return relerr

def chk_T(X,Y,Z,T,time):
    zLength = 6
    R = 0.5
    xr = X / R
    yr = Y / R
    rr = np.sqrt(xr*xr + yr*yr)
    zo = 2 * np.pi * Z / zLength

    Tref = np.sin(np.pi* X) * np.cos(2.7 * np.pi * Z + 0.2)
    Tref[rr>1.0] = 0
    relerr = abs(T-Tref).max() / abs(Tref).max()
    return relerr


# plot
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
err = chk_W(X,Y,Z,W,time)
pc = ax.plot_surface(
        Z, X, W,
        edgecolor="k",
        linewidth=0.3,
        antialiased=True,
        cmap='jet',
    )
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_zlabel('W')
ax.set_title(f'Time={time:g}, relerr={err:.2e}')

fig, ax = plt.subplots()
err = chk_T(X,Y,Z,T,time)
pc = ax.pcolormesh(
        Z, X, T,
        cmap="jet",
        edgecolor='k',
        linewidth=0.0001,
    )
cb = fig.colorbar(pc, ax=ax)
ax.set_aspect("equal")
ax.set_xlabel('Z')
ax.set_ylabel('X')
cb.set_label('T')
ax.set_title(f'Time={time:g}, relerr={err:.2e}')

plt.show()
