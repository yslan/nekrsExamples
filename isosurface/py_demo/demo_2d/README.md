# demo_2d — the 3D isosurface strategy, one dimension down

Contour extraction on a **single (possibly curved) GLL element**: solve `2(N+1)` fibres with
route C, insert lines where an edge holds two roots, march the `N²` sub-cells reading **exact**
crossings, resolve the ambiguous saddle case by *evaluating* the interpolant, and map out through
the element's own geometry nodes.

This is section 6 of the report in miniature — the same pipeline as 3D, in the one dimension
where you can draw it and check it by eye.

> **Self-contained by design.** Every numerical routine is written from its defining recurrence in
> `kernel2d.py`: GLL nodes by Newton on `P'_N`, barycentric weights, the Legendre three-term
> recurrence, Ehrlich–Aberth, bracketed bisection. **No scipy, no `numpy.polynomial`, no LAPACK,
> no third-party contouring**, and no imports from the 1D demos. numpy is used only as an array
> container. This ports directly to C/CUDA, and it means the comparison against the reference
> methods shares no library code.

## Usage

```bash
python3 run_demo.py                      # all studies, writes figures
python3 run_demo.py --study selftest     # verify the hand-written kernel
python3 run_demo.py --study cases        # nine contour cases   (+ figure)
python3 run_demo.py --study accuracy     # vs R1 and R2         (+ figure)
python3 run_demo.py --study geometry     # three geometries     (+ figure)
python3 run_demo.py --study cost         # 2D cost + 3D extrapolation

python3 run_demo.py --case saddle -N 8   # one case, per-fibre detail
python3 run_demo.py --geom annulus       # straight | curved | annulus
python3 run_demo.py --no-fig
```

Cases: `circle` `ellipse` `wavycircle` `twoblobs` `saddle` `tangency` `wave` `band` `corner`.

As a library:

```python
import sys; sys.path.insert(0, 'demo_2d')
import contour2d as C
r, R, S = C.gll_grid(N)
res = C.extract_contour(N, field(R, S), level)
phys = C.map_to_physical(res.segments, Xnodes, Ynodes)
#  res.crossings[d][i]  exact roots on fibre i, direction d
#  res.slopes[d][i]     p'(root): sign = crossing direction, magnitude flags tangency
#  res.inserted         (direction, coordinate) of inserted lines
#  res.ambiguous        saddle cells resolved by evaluation
```

## The pipeline

| step | what happens | cost |
|---|---|---|
| 1. fibre solve | `2(N+1)` rootfinds by route C (Aberth, fixed 12 iterations) | dominant |
| 2. insertion | where an edge holds ≥2 roots, insert a line; sign known *a priori* | ~2% |
| 3. marching | marching squares on `N²` cells, **exact** crossings looked up | table lookups |
| 4. ambiguity | saddle case decided by evaluating the interpolant at the cell centre | 1 eval |
| 5. map out | `(r,s) → (x,y)` by barycentric interpolation of the geometry nodes | `O(N)`/vertex |

**What the fibre solve gives you for free:** exact crossings, the *count* per edge (the thing
corner signs cannot encode), `p'` at each root (Aberth computes `p/p'` anyway — its sign is the
crossing direction, its magnitude flags tangency), and the complete sign structure with **no**
evaluation between roots (between two consecutive roots the sign is constant by definition).

## Results

Figures: `../figs/E_contour_cases.png`, `E_contour_accuracy.png`, `E_contour_geometry.png`
Console: `output.txt`

### Geometry

The fibre solve happens entirely in the **reference** element and is unaffected by deformation —
the roots in `(r,s)` are *identical* on straight, curved, and annulus geometries. Only the
map-out differs. This is not an approximation: the discrete solution **is** the tensor polynomial
in `(r,s)`, so `f∘x` is already what we're solving.

### Accuracy vs the two conventional approaches

Both references use the **same** marching-squares topology written here; only the vertex source
differs (linear interpolation vs exact rootfinding):

- **R1** — treat the GLL nodal values as an `N×N` bilinear grid
- **R2** — resample to a uniform `M×M` grid (`M ≈ 1.5N`), then contour

| case | N | this work | R1 bilinear | R2 resampled | R1/ours |
|---|---|---|---|---|---|
| `circle` | 8 | **1.0e-16** | 1.8e-02 | 5.7e-03 | 1.8e+14 |
| `ellipse` | 8 | **1.0e-16** | 5.0e-02 | 9.2e-03 | 5.0e+14 |
| `ellipse` | 12 | **1.3e-16** | 2.4e-02 | 3.8e-03 | 1.9e+14 |
| `twoblobs` | 8 | 4.0e-02 | 6.1e-02 | 3.6e-02 | 1.5 |
| `wavycircle` | 8 | 5.6e-02 | 3.0e-02 | 4.7e-02 | **0.5** |

**Read this in two halves.** `circle` and `ellipse` are *polynomial* in `(r,s)`, so the
interpolant is exact and the only error left is the contouring — there we are at machine
precision and the references are 10¹³–10¹⁵× worse, limited by linear interpolation along a cell
edge, which resampling shrinks but never removes.

`twoblobs` and `wavycircle` are **not** polynomials. The interpolation error swamps everything
and all three methods tie — ours is even slightly *worse* on `wavycircle`, because placing a
vertex exactly on the wrong curve buys nothing. This is the conditioning floor of report §2 in
2D: **a rootfinder cannot beat the interpolant it is given.** Exact crossings pay off when the
element resolves the field — the regime a converged simulation should be in, and the one where
the reference methods discard the accuracy the discretization paid for.

### Cost (N=8, one element)

```
2D:  fibre solves 2(N+1) = 18      marching cells N² = 64
3D:  fibre solves 3(N+1)² = 243    marching cells N³ = 512
     roots worst case 1944         realistic ~87  (0.36/fibre)
     storage realistic ~1.4 kB     worst case 30.4 KiB  (r+xyz, fp32)
```

## Cross-refs
report §6 (`tex/06_threed.tex`) · `../demo_C` (the fibre solver) · `../notes/05` (conditioning
floor) · `../notes/14` (Aberth initialization)
