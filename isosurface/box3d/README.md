# box3d — the plugin against five cases with known answers

The companion to `tgv/`, and its complement. `tgv/` shows the production call on
a field with no known answer; this case does the opposite — every field here has
a closed-form level set, so the extracted surface can be **checked** — and it
exercises the two debug verbs `tgv/` never calls.

```bash
nrsmpi box3d 1
```

`numSteps = 0`, so nothing is solved: the five fields are constructed on the
mesh, extracted, and compared against their analytic surfaces at `tstep 0`.

## The five cases

| case | what it is for |
|---|---|
| `plane` | degree 1. The interpolant is exact, so **any** error is a bug, not a floor. |
| `sphere` | curved. The error *is* the element's interpolation floor — the right answer is "as small as the discretisation allows", not zero. |
| `poly5` | degree-deficient: the top Legendre coefficient is tiny. This is what exercises **deflation** (`aFloor`); a solver that divides by the leading coefficient fails here. |
| `tangency` | \|∇f\| → 0 at the level. A double root is only √eps-accurate, so the requirement is that the roots are **flagged**, not that they are exact. |
| `tgv` | a real flow field, with an analytic velocity to interpolate onto the surface. |

Levels and centres are derived from the mesh, so a case cannot silently become
vacuous if the mesh changes.

Each runs at both cell orders, which is the only way to see what `cellOrder = 2`
buys: the difference is *inside* the triangles, and at the nodes both orders are
exact by construction.

## What it writes

Per case, per cell order:

```
iso<name>_roots_00_00000.vtu     the crossings, a point cloud   writeRoots()
iso<name>_srf_00_00000.vtu       the surface alone              writeSurface()
iso<name>_fld_00_00000.vtu       surface + fields               process()
iso<name>P2_*.vtu                the same at cellOrder = 2
iso<name>0.f00000                the scalar itself, as a NekRS checkpoint
isok_fld_*.vtu                   the live field (k), two levels, the production shape
isok_fld_<SSSSS>.vtm             groups the two isok levels of each dump (>1 level only)
```

Each verb writes exactly its own file, so a debug stage costs nothing unless
asked for. Only `process()` advances the dump counter, so a debug write shares
the index of the deliverable it belongs to.

Read any of them with `../py_demo/` conventions, or with `iso_read.py` from the
development tree — numpy only, no `vtk` needed.

## When to reach for the debug verbs

- **`writeRoots()`** when the *surface* looks wrong and you need to know whether
  the rootfinding or the triangulation is at fault. It writes the crossings with
  their reference coordinates, both gradients, the owning element and the
  tangency/on-face flags — everything the topology stage consumed.
- **`writeSurface()`** when the geometry is suspect but the fields are not. It
  writes the triangulation with normals and cell flags and skips interpolation
  entirely, so it is also the cheap way to look at a surface.

## Files

`box3d.{udf,par,box,usr}`, plus `isoTestCases.hpp` — the five named
constructors and their checks. Drop that last one in production; it exists to
have something to compare against. The mesh ships compressed as
`box3d.re2.tar.gz` (git); `tar xzf box3d.re2.tar.gz` before running.
