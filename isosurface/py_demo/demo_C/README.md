# demo_C — Route C: matrix-free Ehrlich–Aberth on the Legendre series

Find **all** real roots of `f(x) − c` on a **single GLL reference element** `[-1,1]` of order `N`,
by iterating on all `n` root estimates simultaneously — **without ever forming a matrix**.

This is the **GPU candidate** for the 3D target (`../notes/08`): one fiber per thread,
`3(N+1)²` fibers per element, a **fixed** iteration count, and no branches in the hot loop.

> Route C shares route A's front end exactly — same nodal→modal transform, same degree
> deflation (`../demo_A`). Only the root-extraction step differs: route A builds an `n×n`
> comrade matrix and calls an eigensolver; route C iterates on `n` complex numbers.

## Usage

```bash
python3 run_demo.py                      # every study, orders 3..15, writes figures
python3 run_demo.py --study iterations   # iterations|convergence|cost|gpu|reference
python3 run_demo.py --niter 10           # change the fixed iteration count
python3 run_demo.py --orders 5,7,9,11    # custom order list
python3 run_demo.py --no-fig             # skip figure generation

python3 run_demo.py -N 8                             # single-order detail dump
python3 run_demo.py -N 8 --field sine3               # pick a test field
python3 run_demo.py -N 6 --field gauss --level 0.7   # override the isolevel
```

The `-N` dump prints the deflated coefficients, the **per-fiber state in doubles**, the initial
estimates, a bar-chart iteration trace, the final complex spectrum, and found-vs-exact.

As a library:

```python
from route_c import roots_route_c
r = roots_route_c(N, f)                             # fixed 12 iterations
r, info = roots_route_c(N, f, return_info=True)     # info['final_corr'] = safety net
r, hist = roots_route_c(N, f, track=True)           # per-iteration correction history
```

## The math

Steps 1–2 are route A's: project onto Legendre (`a = M f`, with `γ_N = 2/N`), then deflate
trailing negligible coefficients so `n` is the polynomial's **true** degree.

**Ehrlich–Aberth.** Given estimates `z_1..z_n` of *all* roots at once, each is updated by a Newton
step corrected for the pull of the others:

```
                  w_i                         p(z_i)                      1
  z_i  ←  z_i  −  ──────────  ,     w_i  =  ─────────  ,    S_i  =   Σ   ───────
                  1 − w_i S_i                p′(z_i)               j≠i   z_i−z_j
```

This is Newton applied to `p(x)/∏_{j≠i}(x−z_j)` — the polynomial with the *other* roots divided
out. `S_i` is what stops two estimates collapsing onto the same root, and is why this finds **all**
roots rather than the same one `n` times. Convergence is **cubic** for simple roots.

**Evaluating `p` and `p′` with no matrix.** The same Sturm–Liouville recurrence that *builds* the
comrade matrix, run forwards as an evaluation:

```
P_{k+1}  = ((2k+1)·x·P_k − k·P_{k−1}) / (k+1)
P′_{k+1} = ((2k+1)·(P_k + x·P′_k) − k·P′_{k−1}) / (k+1)
```

accumulating `p = Σ a_k P_k` and `p′ = Σ a_k P′_k`. `O(n)` work, `O(1)` storage, no allocation.
Verified against `legval`/`legder`: **8.9e-16** value, **3.6e-15** derivative.

**Initialization — the one real trap.** The textbook monomial start (`n` points on a circle of
radius ≈ the Cauchy bound) **diverges** here: corrections grew from 1e+00 to 8e+01 over 12
iterations. The roots of interest lie in `[-1,1]`, so start on Chebyshev points nudged off the
real axis:

```
z_i = cos(π(i+½)/n) + 0.1i·sin(π(i+½)/n)
```

The imaginary part matters: with all estimates exactly real and collinear, `S_i` is ill-behaved,
and genuinely complex roots (ATAP p.139–140, the Bernstein-ellipse clustering) are unreachable
from a purely real start.

## "Eigensolver" characterisation

| property | value |
|---|---|
| problem type | **none** — no matrix, no eigenproblem; a simultaneous fixed-point iteration on `n` complex numbers |
| per-fiber state | `a[n+1]` real + `z[n]` complex ≈ **3n doubles** |
| nonzero pattern | n/a. The only all-pairs term is `S_i`, `O(n²)` per iteration, `n ≤ 15` |
| work | `niter · [n evals × O(n) + n repulsions × O(n)]` = **O(niter·n²)** |
| branches in hot loop | **none** — no pivoting, no shifts, no deflation-during-iteration, no convergence test |
| iteration count | **fixed**, default 12 (see below) |
| accuracy | matches route B to **2.0e-13** on every non-tangency case tested |

(For the record: Ehrlich–Aberth *is* mathematically Newton on the comrade matrix's characteristic
polynomial — it solves the same problem route A's eigensolver does, by a different route.)

## Results

Convergence: `../figs/C_convergence.png` · **Iterations: `../figs/C_iterations.png`** ·
Cost: `../figs/C_cost.png` · Baseline: `../figs/00_testfields.png` · Console: `output.txt`

### The central claim: a fixed iteration count is safe

Correction magnitude per iteration, N=8 — cubic, then **stalls** at machine epsilon:

| field | it1 | it2 | it3 | it4 | it5 | it6 | it7 | it8 |
|---|---|---|---|---|---|---|---|---|
| poly5 | 2e-01 | 5e-02 | 4e-04 | 2e-10 | **5e-17** | 3e-17 | 3e-17 | 3e-17 |
| gauss | 2e+00 | 2e+00 | 9e-01 | 1e+00 | 1e-01 | 5e-03 | 3e-07 | **2e-16** |
| sine2 | 5e-01 | 3e+00 | 2e+00 | 4e-03 | 1e-06 | **4e-16** | 4e-16 | 4e-16 |
| sine3 | 2e-01 | 8e-02 | 7e-02 | 5e-02 | 2e-03 | 3e-08 | **3e-16** | 4e-16 |
| tanh | 3e+00 | 2e+00 | 5e-01 | 2e-01 | 9e-03 | 1e-06 | **1e-16** | 9e-17 |

Machine precision in **6–8 iterations**, then flat — so overshooting the count is free, which is
exactly what a fixed-count kernel needs.

Over **20,000 random fibers**, against route B:

| niter | root-count mismatches | max \|C − B\| |
|---|---|---|
| 6 | 57 | 5.7e-08 |
| 8 | 2 | 4.0e-08 |
| **10** | **0** | **2.0e-13** |
| **12 (default)** | **0** | **2.0e-13** |
| 24 | 0 | 2.0e-13 |

On the 8 analytic fields × orders 3–15, the only disagreements are **11/104, all on `tangent`** —
a double root, √eps-conditioned for *every* method (`../notes/05`). Zero non-tangency failures.

### Why not just write a QR? (measured, not assumed)

A hand-written complex single-shift Givens QR on route A's last-column comrade form — which is
**already upper Hessenberg** for free (ATAP Ex. 18.2), so no reduction is needed — is **39
non-comment lines** and does work. But at equal fixed cost it loses. 4000 random fibers vs route B:

| method | root-count mismatches | max \|·−B\| |
|---|---|---|
| QR, 4 sweeps/root | 169 | 1.4e+00 |
| QR, 6 sweeps/root | 25 | 2.3e-07 |
| QR, 8 sweeps/root | 1 | 2.2e-08 |
| **Aberth, 10 iters** | **0** | **7.1e-14** |
| **Aberth, 12 iters** | **0** | **7.1e-14** |

QR needs ~8 sweeps/root to approach reliability and *still* doesn't reach zero, while costing
2.5–7.5× more flops — and it carries a Wilkinson shift formula, a deflation criterion (a
data-dependent **branch**), and per-fiber-varying work. The conclusion is not "a QR is too hard to
write"; it is **"a QR is the wrong tool for n ≤ 15 on a fixed budget."**

### Cost: read the compiled number, not the python one

**The `--study cost` table in `output.txt` is misleading and is kept only for continuity.** It
reports ~500–1200 µs/solve for route C against ~28–52 µs for routes A and B — *18–25× slower for
the route with the fewest flops*. That is an artifact of packaging, not algorithms: route C is
~600–1600 interpreter-level numpy calls on size-8 arrays, while routes A and B each make **one**
compiled LAPACK call. The implied flop rate gives it away — route C at N=5 reads **1.0 Mflop/s**,
when a core does 1000–10,000. See `../notes/12`.

Three measurements of the same kernel, N=5 → N=15:

| N | python scalar (the demo) | python batched over fibers | **compiled C, 1 core** |
|---|---|---|---|
| 5 | 485 µs | 9.2 µs | **3.10 µs** |
| 8 | 673 µs | 16.3 µs | **7.71 µs** |
| 11 | 913 µs | 30.6 µs | **14.80 µs** |
| 15 | 1180 µs | 58.8 µs | **27.63 µs** |

`aberth_bench.c` (in this directory) is the compiled version — it mirrors `route_c.py` exactly.
Build and run:

```bash
gcc -O3 -march=native -o aberth_bench aberth_bench.c -lm
./aberth_bench 8 200000        # degree, repetitions
```

At 3.1 µs/fiber (N=5) it sustains **4.9 Gflop/s** single-core — a realistic rate, confirming the
compiled kernel is arithmetic-bound, which is what a GPU port needs. 322,583 fibers/s on one core;
this mesh needs 93,312 fibers at N=5.

> **Flop-model caveat.** The `O(niter·n²)` model used in the tables below counts *complex*
> operations as single flops. Real flops (`c*c`=6, `c/c`≈11, `c±c`=2) are **25× higher**: 14,640
> at n=5, not 600. The *scaling* is right, the constant is not — so the "2.1–6.3× fewer flops than
> QR" ratio below is computed in the wrong currency and needs recounting on both sides before it
> means anything. `../notes/12`.

### GPU budget (one fiber per thread)

A CUDA thread holds ~255 32-bit registers ≈ 128 doubles before spilling:

| N | **route C (~3N)** | route A (N²) | route B (2(N+2)²) | fibers/elem | total fibers (864 el) |
|---|---|---|---|---|---|
| 5 | **16** | 25 | 98 — spill | 108 | 93,312 |
| 8 | **25** | 64 | 200 — spill | 243 | 209,952 |
| 11 | **34** | 121 | 338 — spill | 432 | 373,248 |
| 15 | **46** | 225 — spill | 578 — spill | 768 | 663,552 |

Work per fiber:

| N | route C (Aberth, 12) | dense QR (~10N³) | ratio |
|---|---|---|---|
| 5 | 600 | 1,250 | 2.08× |
| 8 | 1,536 | 5,120 | 3.33× |
| 11 | 2,904 | 13,310 | 4.58× |
| 15 | 5,400 | 33,750 | 6.25× |

Both columns are idealized models in the same (undercounted) currency — see the flop-model caveat
above. The robust statement is the **scaling**: `O(niter·n²)` vs `O(n³)`, i.e. the `O(N²)` reached
by abandoning the matrix rather than by a structured eigensolver (`../notes/06`). Treat the ratio
column as indicative only until both sides are recounted in real flops.

## Kernel sketch (~40 lines, branch-free)

```
a = M·f                       # tensor GEMM, once for the whole field — not per thread
n = deflate(a, 1e-12)         # O(N), one pass
z_i = cos(π(i+½)/n) + 0.1i·sin(π(i+½)/n)
repeat 12×:                   # fixed; no convergence test
    (v,d) = legendre_val_der(a, z)    # O(n) recurrence per root
    w = v/d ;  S_i = Σ_{j≠i} 1/(z_i − z_j)
    z -= w/(1 − w·S)
keep real z in [−1,1]
```

## Limitations

- **No theorem** guarantees all roots at a *fixed* count. Empirically 0/20,000 at niter ≥ 10, but a
  production kernel should keep the safety net: `return_info=True` reports `final_corr`, and a
  large value means that fiber did not settle.
- **Double roots** are √eps-accurate, as for every method. On an isosurface these are *tangencies* —
  detect them (a returned pair separated by ~1e-8) rather than trusting either root.
- **Complex arithmetic** doubles the register count versus a real-only scheme. A real-only variant
  (Bairstow-style quadratic factors) avoids that but reintroduces branches.
- **The demo's µs/solve does not rank the routes** — it compares python dispatch counts. Use
  `aberth_bench.c` (`../notes/12`).
- Validation so far is against **synthetic** data. Synthetic data has hidden structure from this
  project twice (`../notes/03` parity, `../notes/06` `fᵀb`); real fibers from `../iddes.fld` are
  the outstanding check.

## Cross-refs
`../notes/05` (accuracy floor, double roots) · `../notes/07` (deflation) · `../notes/08` (3D/GPU
framing) · `../notes/10` (in-house eigensolver analysis) · `../notes/11` (the QR measurement that
corrected it) · `../notes/12` (why the demo µs/solve is an artifact; compiled numbers) · `../demo_A` (shared front end) · `../demo_B` (the CPU oracle) · ATAP ch. 18,
Ex. 18.2
