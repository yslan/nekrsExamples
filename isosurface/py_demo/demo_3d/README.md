# demo_3d — three hard isosurface cases on one deformed element

Each case is shown at the **four stages** of the §6 pipeline, so the method can be inspected
rather than taken on trust:

| stage | what it shows |
|---|---|
| **(1) roots** | the raw output: exact fibre crossings on the GLL lattice, coloured by direction |
| **(2) triangles** | marching cubes on the `N³` sub-cells — flat facets, ambiguous cells counted |
| **(3) reconstructed** | the curved surface rebuilt from the **stored payload only** — vertex positions + unit normals, **no further access to the element data** |
| **(4) exact** | the analytic level set, finely sampled: the reference |

Figures: `../figs/F_case_tunnel.png`, `F_case_gyroid.png`, `F_case_neck.png`. Console: `output.txt`.

## The three cases, and why these

Chosen from what the visualization literature treats as the hard configurations
(Chernyaev's MC33; Nielson–Hamann's asymptotic decider; Lewiner et al.).

**`tunnel` — `rs + st + tr = −0.18`.** The Marching Cubes **case 13** configuration: the
element's corner signs alternate `(+,−,−,−,−,−,−,+)`, which is exactly where face ambiguity and
interior ambiguity interact. This is the case MC33's ~730-entry subcase tables exist for. Here
the interpolant is available, so the ambiguity is settled by **evaluating it** — 18 ambiguous
cells at N=8, no subcase table. The tunnel is captured.

**`gyroid` — `sin r cos s + sin s cos t + sin t cos r = 0`.** A triply-periodic minimal surface:
several disconnected sheets threading one element, high genus, and the most ambiguous cells of
any case here (**81**). The standard stress test for reproducing topology rather than merely
producing triangles.

**`neck` — `((r−0.46)² + s² + t²)((r+0.46)² + s² + t²) = 0.045`.** Two lobes at the instant of
merging. `min |∇f| = 2.0e-02` against O(1) elsewhere — the √eps tangency of report §2 in 3D, and
the one case where gradient projection must be guarded.

## Results

```
  case      roots   tris  ambig  min|grad|   payload |   flat    PN recon  projected    floor
  tunnel      180    294     18   6.29e-01      7.0K | 4.31e-02  2.65e-02   3.79e-03  7.92e-08
  gyroid      333    528     81   3.25e+00     12.7K | 1.67e-01  8.10e-02   1.18e-02  4.05e-04
  neck         14     24      0   1.96e-02      0.6K | 3.05e-02  1.45e-02   7.30e-03  2.91e-05
```

All four errors are `max |f − c|` against the **analytic** field, sampled **inside** every
triangle — not at the nodes, where the answer is trivially exact.

- **flat** — flat marching-cubes triangles
- **PN recon** — curved from the stored payload alone (**the operational constraint**)
- **projected** — subdivide + Newton on the interpolant; needs element access, kept as reference
- **floor** — distance from the *interpolant's* level set to the analytic one: the discretization
  error, which no surface method can beat

**Reading it honestly.** PN reconstruction roughly **halves** the flat error at ~5 floats/vertex
and zero element access. Projection does 7–11× better but needs the element data. On `tunnel`
both remain far above the interpolation floor (7.9e-08), so there the limit is the
**reconstruction**, not the discretization — a richer payload would pay. On `gyroid` the floor is
4.1e-04 and the field is barely resolved at N=8, so raising N matters more than refining the
surface.

## Usage

```bash
python3 run_demo.py                  # all cases, 4-stage figures
python3 run_demo.py --case tunnel    # one case
python3 run_demo.py --study report   # numbers only
python3 run_demo.py -N 8 --nsub 2    # element order, subdivision depth
```

## Reconstruction (the constraint)

Stored per vertex: position (3 floats) + unit normal (2 after normalization). Curvature is
recovered by the PN-triangle rule — purely geometric, no field access:

```
m = (p0+p1)/2 ;   m ← m − ((m−p0)·n0)n0/2 − ((m−p1)·n1)n1/2
```

The normal is free at extraction: the along-fibre derivative is `p'(r*)` from Aberth, the
transverse ones come from a Lagrange contraction of nodal values already held
(`../demo_2d/refine2d.py::tangent_from_stored`).

## Cross-refs
report §6 · `../demo_2d` (the 2D pipeline and the `surface3d.py` extractor this uses) ·
`../notes/05` (the √eps tangency) · Chernyaev MC33 · Lewiner et al.
