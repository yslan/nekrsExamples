# py_demo — the method, one piece at a time

Five self-contained Python demos. They are how the algorithm in `isoSurface.okl`
was designed and checked: each one isolates a single decision, measures it, and
writes the figure that justifies it. Nothing here is needed to *use* the plugin —
read these to understand or modify it.

```bash
python3 demo_C/run_demo.py            # start here
python3 demo_3d/run_demo.py --case turbulent
```

Each demo has its own `README.md` with the detail; this page says which to read.

## Where to start

| | what it isolates |
|---|---|
| **`demo_C`** | **the shipped algorithm.** Matrix-free Ehrlich–Aberth on the Legendre series: no matrix, ~3N doubles per fibre, fixed 12 iterations, no branches. |
| **`demo_3d`** | **the 3D strategy.** Five hard cases on one element, each at four stages: roots → marching cubes → payload-only reconstruction → analytic reference. |
| `demo_2d` | the same strategy one dimension down, where you can see it. Fibre solves → marching squares with *exact* crossings → isoparametric map-out. |
| `demo_A` | route A: the modal Legendre comrade matrix + degree deflation. The GPU candidate before route C won. |
| `demo_B` | route B: the nodal barycentric pencil. Slower, never divides by a leading coefficient — kept as the CPU oracle everything else is checked against. |

`demo_A`/`demo_B`/`demo_C` take `-N <order>` (3–15), `--field`, `--level` and
`--study`; they share `src/gllbase.py` and `src/studies.py`.

## Layout

```
src/         gllbase.py (GLL nodes, barycentric weights, nodal<->modal),
             studies.py (the convergence and cost sweeps)
demo_*/      one directory per demo, each with its own README
figs/        where a run writes its figures (--figdir, default ../figs)
```

`figs/` here is deliberately separate from `release/figs/`, which holds the four
curated figures the top-level README shows — a demo run would otherwise
overwrite them.

## Two honest notes

- **`demo_A` and `demo_B` need `scipy`** (`scipy.linalg` for the eigensolvers
  they exist to demonstrate). `demo_C`, `demo_2d` and `demo_3d` are numpy-only,
  which is deliberate: those three implement GLL, Lagrange, Legendre and Aberth
  from their definitions, so nothing in the measured path comes from a library.
- **`demo_C/aberth_bench.c` is a separate C build**, not run by the Python. It
  exists because a Python timing of a per-fibre solve measures the interpreter,
  not the algorithm — see `notes/12`.

  ```bash
  gcc -O3 -march=native -o aberth_bench demo_C/aberth_bench.c -lm && ./aberth_bench 8 200000
  ```

## What the demos are evidence for

The write-up in `../doc/main.pdf` cites them throughout; its figures are the ones
these scripts produce. If a claim in the README or in `DEVELOPER.md` has a number
attached, a demo here is usually where it came from.
