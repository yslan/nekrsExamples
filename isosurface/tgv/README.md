# TGV — a Taylor–Green vortex with a live Q-criterion isosurface

The smallest case that shows the isosurface plugin doing its job: as the vortex
sheets roll up and break down, the Q-criterion isosurface is extracted on the
GPU every checkpoint step and written as a `.vtu` carrying velocity and pressure
sampled onto it.

```
nrsmpi tgv 1                 # 1 rank, 1 GPU
```

The example writes **two** levels, Q = 0.1 and 0.01, which share one solve pass:
the fibre roots are found once per element and every level reads them.

## What you get

```
isoQ_fld_00_00000.vtu        level 0 (Q = 0.1)    surface + fields  <- the deliverable
isoQ_fld_01_00000.vtu        level 1 (Q = 0.01)
isoQ_fld_00000.vtm           groups both levels of this dump: open this one
tgv0.f00001, ...             the usual NekRS field files (Q is in them as scalar00)
```

One `.vtu` per level per dump plus one `.vtm` grouping them, the first written at tstep 0
from the initial condition. Open the `.vtm` in ParaView to see both Q levels as one object.

Open the `.vtu` directly in ParaView and colour by `velocity`. At `cellOrder = 2`
the cells are `VTK_QUADRATIC_TRIANGLE`, which ParaView tessellates natively —
about 25× more accurate inside each triangle than flat cells, for ~2.5 extra
points per corner.

Or read it with numpy alone, no `vtk` required:

```python
from iso_read import read_fields
d = read_fields("isoQ_fld_00_00000.vtu")
d["Position"], d["fields"]["velocity"], d["tri"]
```

## The one setting to know about

```ini
checkpointControl  = simulationTime
checkpointInterval = 2.0     # in tgv.par
```

The `.udf` writes an isosurface on every checkpoint step, so this drives **both**
the `.fld` series and the isosurface series. With `endTime = 20` and `dt = 2e-3`
that is **10 isosurfaces**, at t = 2, 4, ... 20 — the breakdown at t ≈ 9 and the
decaying turbulence after it are both covered. The upstream `examples/tgv` case
has `checkpointInterval = 0` — checkpoints off — and with that value this case
runs to completion and writes *nothing*.

## The mesh

`tgv.re2` ships with the case (compressed as `tgv.re2.tar.gz` for git), so
nothing needs generating — just `tar xzf tgv.re2.tar.gz` first. It is 8×8×8
elements on `[-π,π]³`, fully periodic: 512 elements, 0.26M points at
`polynomialOrder = 7`.

To change the resolution, edit `input.box` and re-run `genbox` (a Nek5000 tool,
**not** part of NekRS and not in `$NEKRS_HOME/bin`):

```bash
echo input.box | genbox            # -> box.re2
mv box.re2 tgv.re2
```

`base.rea` ships with the case: genbox reads it for the parameter template only,
not for geometry, so the same file serves any box.

On this machine the working binary is `/Disk_D/home/sean/Nek5000_v19/bin/genbox`;
the v17 one fails to load `libgfortran.so.4`.

## The figure

```bash
python3 plot_tgv.py            # both; ~2 min for the .fld one
python3 plot_tgv.py --plugin   # just the plugin surface
```

Reproduces figure 1 of the C3.3 reference — Q = 0.1 coloured by |V|/V0, at
t = 0, 4, 8, 12, 16, 20 — twice:

| file | source | how the surface is found |
|---|---|---|
| `tgv_Q_nek.png` | `tgv0.f*` | GLL nodes → linear hexes → VTK marching cubes |
| `tgv_Q_plugin.png` | `isoQ_fld_00_*.vtu` | rootfinding on the polynomial, P2 triangles |

The first is what a normal post-processing tool does with a spectral element:
the nodes become points, each sub-cell a linear hex, the polynomial between them
discarded, and each crossing interpolated along an edge. The second is the
plugin's own output, whose vertices are on the level set to machine precision.

The two agree to **0.35% in area at t = 0**, with an identical triangle count
(37184), and to 4% at t = 4. The t = 0 panel is the meaningful check: it is a
smooth analytic field, so a real disagreement would show up there. The t = 4 gap
is the Q field being genuinely hard at this resolution, not a difference between
the methods.

This mesh is 8³ where the reference is 257³, so expect far more small-scale
structure than the paper's figure.

`nek_read.py` reads the `.fld` with numpy alone — VTK only gained a Nek5000
reader in 9.3, and this box has 9.2.2. It has a self-test:

```bash
python3 nek_read.py --check    # pointwise against the analytic TGV IC
python3 nek_read.py            # list the dumps and their times
```

That test earned its place. The `.fld` layout has two traps — the min/max block
is a **trailer**, not a header, and a vector field is written per element with
its components interleaved — and getting either wrong still yields coordinates
on [-π,π] and a plausible-looking picture. Only a pointwise comparison against
`u = sin(x)cos(y)cos(z)` catches it: correct reads match to 1.1e-07 (float32),
a wrong one correlates at 0.05.

Both scripts read the **current directory** by default, so running them from
inside a run folder just works. An explicit argument
(`python3 plot_tgv.py dat_E512`) or `$TGV_DATA` overrides that, and if the cwd
holds no output they fall back to a `dat_*/` subdirectory.

## What this case is, and is not

It is a **plugin usage guide**. It is *not* a validated Taylor–Green solution:
8³ elements is far too coarse for Re = 1600, and the `hpfrt` regularization, the
energy/enstrophy diagnostics and the coordinate shear of the upstream case have
all been removed to keep the `.udf` readable. For the physics, see
`examples/tgv` and the C3.3 reference (Twente).

## Using the plugin in your own case

Three things:

```cpp
// 1. the two includes
#ifdef __okl__
#include "isoSurface.okl"
#endif
#include "isoSurface.hpp"

// 2. UDF_Setup -- build the extractor on whatever scalar you want a surface of
IsoSurfaceParams prm;
prm.levels    = {yourLevel};   // several levels share one solve pass
prm.fileName  = "isoQ";
prm.cellOrder = 2;
iso = std::make_unique<isoSurface>(nrs->meshV,
                                   isoSurface::field{"Q", {o_Q}}, prm);
iso->setOutputFields({{"velocity", {o_ux, o_uy, o_uz}}});

// 3. UDF_ExecuteStep -- when the scalar has changed
iso->process(time, tstep, true);   // true = Q moved, re-solve the surface
```

`process()` writes the surface-plus-fields file. `writeRoots()` and
`writeSurface()` are debug products and each writes only itself — there are no
dump flags, the caller picks a verb.

**While developing**, note that editing `isoSurface.hpp` or `isoSurface.okl` does
**not** trigger a NekRS rebuild (only the `.udf` and `<case>.okl` are checked).
`touch tgv.udf` or run with `--build-only`, or you will be testing the previously
built kernels.
