# isoSurface — a developer's map

For modifying `isoSurface.hpp` and `isoSurface.okl`. To *use* the plugin, read
`README.md`; for why the design is what it is, `doc/main.pdf` §8.

Everything below is in namespace scope in `isoSurface.hpp` unless marked OKL.

## The pipeline

```
   scalar on the device
        |
   isoCullOKL            drop elements whose nodal range cannot reach the level
        |
   isoSolveFibresOKL     3(N+1)^2 fibre solves/element -> roots in (r,s,t)
        |
   isoCompactPointsOKL   gather the sparse per-element slabs
        |
   (host)                physical x,y,z and both gradients; the MultiEdge flag
        |
   isoTopologyOKL        marching cubes on the GLL lattice + the ring walk
        |
   (host)                P2 mid-side node creation: the unique-edge dedup and
                         the face-lock mask
        |
   isoProjectOKL         Newton-project each mid-side node onto the level set
        |
   isoInterpOKL          sample the registered fields at every surface point
        |
   .vtu                  collective MPI-IO, appended binary
```

The kernels work **entirely in the reference element** — that is what makes the
rootfinding exact, and it is why the host still has work to do.

## Data structures

| struct | holds | filled by | read by |
|---|---|---|---|
| `IsoSurfaceParams` | every tunable; see `README.md` | the caller | everything |
| `RootSet` | one level's crossings: `(r,s,t)`, `(x,y,z)`, both gradients, the owning element, fibre `(dir,j,k)`, flags | `solveRoots` | `buildTopology`, `writeRoots` |
| `SurfaceSet` | one level's triangulation: points, normals, flat connectivity, `cellType`/`nodesPerCell` | `buildTopology` | `interpFields`, both surface writers |
| `FieldSet` | the sampled field values, one array per component | `interpFields` | `writeFields` |
| `IsoStats` | counters: roots, cells, ambiguous faces, tangencies, per-element occupancy | every stage | the report lines |
| `IsoCaps` | the three per-element capacities | `setupCommon` | every kernel launch |
| `isoDetail::AberthOpts`, `CubeEdge` | solver options; the generated cube-edge table | — | the numerics |

`RootSet` and `SurfaceSet` are **structure-of-arrays**: the writers pack SoA
straight into the appended VTK blocks.

Points in a `SurfaceSet` come in two groups, in order: `[0, nCorner)` are
CORNERS and *are* roots (`cornerRoot[i]` indexes the `RootSet`); `[nCorner,
nPoint)` are the P2 mid-side nodes.

## Public API

```cpp
isoSurface(mesh_t*, const field& scalar, const IsoSurfaceParams&);
isoSurface(mesh_t*, const std::vector<dfloat>& scalar, const IsoSurfaceParams&);
void setField(const field&);            // or a host vector
void setOutputFields(const std::vector<field>&);
void setOutputFieldHost(const std::string&, int nComp, const std::vector<dfloat>&);
```

`field` is `std::tuple<std::string, std::vector<deviceMemory<dfloat>>>`, the
NekRS plugin convention. The extraction scalar must be ONE component; an output
field may be several (`{"velocity", {o_ux, o_uy, o_uz}}`).

**The scalar is held by HANDLE.** `deviceMemory`'s copy constructor copies the
reference-counted `occa::memory`, so the plugin aliases the caller's buffer and
re-reads it at the top of every solve — a field overwritten in place needs no
re-registration. The plugin cannot *detect* that the values changed, which is
what `process(t, step, true)` is for.

**Three verbs, three products.** There are no dump flags; the caller chooses:

| verb | writes | |
|---|---|---|
| `process(t, step[, solveGeom])` | `_fld_` | the deliverable: surface + fields |
| `writeSurface(t, step)` | `_srf_` | debug: geometry, normals, cell flags |
| `writeRoots(t, step)` | `_roots_` | debug: the raw crossings as a point cloud |

With more than one level each verb ends with `writeVtm(product, counter)`: rank 0 writes
`<fileName><tag><SSSSS>.vtm`, a `vtkMultiBlockDataSet` whose block `i` is level `i` (name
`level_<LL> <field>=<value>`, file = the `.vtu` basename). It runs after the per-level loop, i.e.
after every collective `MPI_File_close`, so it never names a missing file. The XML is the pure
function `isoVtmXml`, which `plugin/test_vtm.sh` extracts from the header and tests.

The stages are also callable directly — `solveRoots(force)`,
`buildTopology(force)`, `interpFields()` — and each skips its work if its input
has not changed. Accessors `roots(l)`, `surface(l)`, `fields(l)`, `stats()`,
`params()`, `caps()` expose the state without copying.

## The eight kernels (OKL)

**`isoCullOKL`** — one element per block, `(N+1)^2` threads reduce the nodal
min/max. An element is kept when the level lies in
`[min - θ·K_N·R, max + θ·K_N·R]`. `K_N = (Λ_N−1)/2` is the exact worst-case GLL
overshoot as a fraction of the nodal range, so θ = 1 cannot miss a crossing;
θ = 0 (the naive test) can, because an interpolant really does leave its nodal
range.

**`isoSolveFibresOKL`** — the hot kernel. **One fibre per thread**, `@inner(0)`
spanning all `3(N+1)^2`. Register pressure decides that shape: the Aberth state
is `6·Nq` doubles, and unrolling the three directions risks 216–288 registers
against an architectural 255. Two passes over the fibres, because OKL has no
atomic capture — count into `@shared`, prefix-sum, write at an exclusive offset.
That also makes the output order deterministic. `@shared`: `s_fld[Np]`,
`s_nRoot[3(N+1)^2]`.

**`isoTopologyOKL`** — one element per block. The marching lattice IS the GLL
grid, so crossings are already exact; no edge interpolation. Cut points are
connected by the **ring walk**, not a triangle table: each face contributes
segments, and chaining them closes one or more rings, each fanned separately.
A face is shared by exactly two cells and both see its own four edges, so the
pairing cannot disagree across a cell boundary. Where four crossings make the
pairing ambiguous, the interpolant is **evaluated at the face centre** to decide.
`@shared`: `s_fld`, `s_edgePt`, `s_edgePt2` (the second crossing on an edge),
`s_triCount` — about 18 kB at N=7 against a typical 48 kB budget.

**`isoProjectOKL`** — the P2 mid-side nodes, fixed `niterProj` Newton steps in
the reference frame, no early exit so a warp never diverges. `midLock` is a
3-bit mask: a node whose parents lie on a shared element face has that
coordinate pinned, and the locked gradient components are zeroed **before** the
step — which keeps it an exact Newton step on the restricted problem, hence
quadratic. Zeroing afterwards would be only linear.

**`isoInterpOKL`** — samples every registered field at every surface point. The
basis weights are field-independent, so they are computed once per point and
reused across all fields and components.

**`isoCompactPointsOKL` / `isoCompactCellsOKL` / `isoCompactFieldOKL`** — gather
the sparse per-element slabs into dense arrays before the host copy.

## Host / device split

The interface is 18 `deviceMemory` buffers, allocated once in `allocDevice()`
and reused. Four things stay on the **host** because no kernel produces them:

1. **Physical coordinates and both gradients.** The kernels emit only reference
   `(r,s,t)`; the isoparametric map and the `J⁻ᵀ` rotation run here.
2. **P2 mid-side node creation** — the unique-edge dedup is a hash;
   `isoProjectOKL` takes `midCount`/`midLock` as given.
3. **Lattice refinement**, where v1 inserts a plane.
4. **The `MultiEdge` flag.**

A split that dropped (1) would write a well-formed `.vtu` full of zeros.

## Capacities

OKL has no dynamic allocation, so every per-element output needs a fixed bound.
These are **derived from the algorithm, not measured** — that distinction cost
two production aborts, and lowering them is a buffer overrun, not a tuning
choice.

| | value | why |
|---|---|---|
| `capRoot` | `capScale · 3N(N+1)²` | `3(N+1)²` fibres, each a degree-N restriction with at most N real roots |
| `capCell` | `p_isoMaxTri · N³` | each sub-cell owns a FIXED slice, so a smaller value writes past the buffer |
| `capPoint` | `capRoot + 3·capCell` | corners, plus up to three mid-side nodes per triangle |

`p_isoMaxCut = 24`, **not 12**: an edge whose endpoints share a sign can carry
*two* crossings, and since the ring walk both are kept.
`p_isoMaxTri = p_isoMaxCut − 2 = 22`, **not** the textbook marching-cubes 5 —
this scheme has no triangle table, and a ring of `np` points emits `np − 2`.

Overflow is counted, not silent: the run aborts and names the stage.

## Traps

- **Editing this file does not trigger a rebuild.** NekRS mtime-checks only
  `<case>.okl` and the `.udf`. `touch` the `.udf` or run `--build-only`.
- **Four things OKL rejects that g++ accepts**: a cast on the right of a binary
  operator needs parentheses; there is no atomic *capture*; CUDA/HIP need
  `@atomic x += 1`, not `x = 1`; a `@shared` scalar must be an array. A
  post-increment also cannot index (`walk[nw++]`).
- **`deviceMemory` carries its dtype** and throws on a mismatch, so an index
  buffer must be `dlong`, not `int` — they are the same type in some builds and
  not others.
- **Derived tables must be derived.** The cube-edge, face-edge and marching
  tables are generated and self-verified. Restating one by hand is how three
  separate bugs were introduced.

## Testing without a GPU

`plugin/` holds the development tree: a second host implementation kept as the
oracle, five analytic cases, `okl_extract.py` (compiles the OKL as ordinary C++,
proving arithmetic and indexing but *not* the absence of races), and two checks
against the real toolchain — `host_compile_check.sh` (real NekRS headers plus
udfMake's shim) and `okl_translate_check.sh` (the real OCCA parser, 20
backend × order × precision combinations, then `occa compile` of all eight
kernels). A hand-written stub proves the shape of an API and nothing more; both
of those checks exist because a stub let a real failure through.
