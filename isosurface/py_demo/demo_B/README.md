# demo_B — Route B: nodal barycentric arrowhead pencil

Find **all** real roots of `f(x) − c` on a **single GLL reference element** `[-1,1]` of order `N`,
by building a generalized eigenvalue **pencil directly from the nodal values** — no basis
transform. This is the "eigenvalue problems constructed from function values themselves" idea
Trefethen flags as under-exploited at *ATAP* p.139 (Amiraslani–Corless–Lancaster 2004/2006).

Source mesh is always one element; `N` is the user's knob (tested 3–15, SEM-typical 5–10).

## Usage

```bash
python3 run_demo.py                      # every study, orders 3..15, writes figures
python3 run_demo.py --study identity     # identity|convergence|cost|degeneracy|reference
python3 run_demo.py --orders 5,7,9,11    # custom order list
python3 run_demo.py --no-fig             # skip figure generation

python3 run_demo.py -N 8                             # single-order detail dump
python3 run_demo.py -N 8 --field sine3               # pick a test field
python3 run_demo.py -N 6 --field gauss --level 0.7   # override the isolevel
```

Fields: `poly5 tangent gauss tanh runge sine2 sine3 sine5` (see `../figs/00_testfields.png`).
The `-N` dump prints the nodes, the barycentric weights **and their spread**, the assembled
arrowhead `A` and `diag(B)` (for N ≤ 8), the finite/infinite eigenvalue split, the full spectrum,
and found-vs-exact.

As a library:

```python
from route_b import roots_route_b
r, info = roots_route_b(N, f, return_info=True)   # info['n_finite'], info['n_infinite']
```

`../src/gllbase.py` holds the GLL basis, test fields and error metric; `../src/studies.py` holds
the studies and plotting, shared with `demo_A`.

## The math

Expand in the **Lagrange cardinal basis** on the GLL nodes, `ℓ_j(x_i) = δ_ij`:

```
p(x) = Σ_j f_j ℓ_j(x)
```

The expansion coefficients **are** the nodal values — no projection, no `γ_N` subtlety, nothing to
get wrong. That is the appeal.

**The obstruction.** `{ℓ_j}` has *no* three-term recurrence: `x·ℓ_j(x)` has degree `N+1` and is not
a combination of the `ℓ_i`. So "multiply by x" is not a matrix on `span{ℓ_j}` — step 2 of the
colleague recipe fails outright.

**The fix.** Allow a *generalized* eigenproblem `A v = λ B v` with `B` singular. With barycentric
weights `b_j = 1/∏_{i≠j}(x_j − x_i)` and node polynomial `ℓ(x) = ∏_i (x − x_i)`,

```
p(x) = ℓ(x) · Σ_j  b_j f_j / (x − x_j)
```

so roots of `p` are zeros of the **secular function** `s(x) = Σ_j b_j f_j/(x − x_j)`. Clearing
denominators of a secular equation is the classical **arrowhead** eigenproblem:

```
      ⎡ 0    f_0  f_1  …  f_N ⎤            ⎡ 0              ⎤
      ⎢ b_0  x_0              ⎥            ⎢    1           ⎥
  A = ⎢ b_1       x_1         ⎥      B  =  ⎢       1        ⎥
      ⎢  ⋮            ⋱       ⎥            ⎢          ⋱     ⎥
      ⎣ b_N               x_N ⎦            ⎣             1  ⎦
```

Nodes on the diagonal, nodal values across the top row, barycentric weights down the first column.
`B` is the identity with a **single zero at (0,0)**.

**Why it is exact.** Expanding the arrowhead determinant along the first row/column,

```
det(A − λB) = − Σ_j b_j f_j ∏_{i≠j} (x_i − λ)  =  ∓ p(λ)
```

— the pencil's characteristic polynomial *is* the interpolant. An identity, same status as
Thm 18.1. Verified by `--study identity`: `det/p` is **constant in λ** to
`max|normalized − 1| = 7.9e-14` over N ∈ {4,6,8,11,15} × 6 values of λ.

> **Scale caveat.** `gll_bary_weights()` returns the closed form `(−1)^j √w_j`, which is
> proportional to the true barycentric weights, not equal (κ = 2.77 at N=4, −3241 at N=15). The
> pencil is homogeneous of degree 1 in `b`, so **every root is unaffected**; only the determinant
> carries κ, which is why the identity study normalizes. See `../notes/02`.

The isolevel goes in the **top row** as `f_j − c`. Do *not* shift the diagonal — it holds nodes.

**Root exactly at a node** looks fatal (the secular function has an apparent pole) but is fine:
when `f_j = 0` the residue `b_j f_j` vanishes with it, and the pencil is assembled from `f_j` and
`b_j` directly, never from a quotient. Verified to 1e-15 (`../notes/02`).

## Eigenproblem characterisation

| property | value |
|---|---|
| problem type | **generalized**, non-symmetric, **singular B**: `A v = λ B v` |
| size | `(N+2) × (N+2)` — larger than route A's `N × N` |
| structure of A | **arrowhead** (bordered diagonal): diagonal `x_j`, first row `f_j`, first column `b_j` |
| structure of B | diagonal `(0,1,…,1)`, exactly singular |
| nonzeros | `nnz(A) = 3N + 3` (measured): 18/49 at N=5, 27/100 at N=8, 48/289 at N=15; `nnz(B) = N+1` |
| symmetry | **none** (row holds `f`, column holds `b`); symmetrizable only if all `f_j b_j > 0`, which the data does not guarantee |
| eigenvalues | `N+2` total; **≥1 infinite always**, plus one more per degree deficiency; finite ones = roots |
| LAPACK path | `scipy.linalg.eig(A,B)` → `dggev` = `dggbal` + `dgghrd` + `dhgeqz` (QZ) |
| infinity filter | **relative** test on the homogeneous pair: `homogeneous_eigvals=True`, keep `|β| > tol·|α|`. Never threshold `|λ|` by magnitude |
| flop model | `dgghrd` ~(20/3)n³ + `dhgeqz` ~30n³ ≈ **37n³**, `n = N+2` |

## Complexity

`O(N³)` asymptotically — but as with route A, **not what is measured at SEM orders**
(`../notes/04`). Fitted exponent over N=3…15 is **p = 0.44**. The flop model says route B should
cost ~7.8× route A at N=5; measured is **1.2–1.5×**, because both are dominated by call overhead,
not arithmetic. Breakdown at N=8: build pencil 2.50 µs, `eig(A,B)` 20.24 µs, total 32.9 µs.

Arrowhead pencils admit `O(n²)` (or `O(n log n)` with fast multipole) structured solvers via
secular deflation — *more* exploitable than route A's comrade matrix. Not implemented; pointless
at these sizes.

## Results

Baseline: `../figs/00_testfields.png`, `../figs/00_testfield_resolution.png`
Convergence: `../figs/B_convergence.png`  ·  Cost: `../figs/B_cost.png`
Full console output: `output.txt`

### Max root error vs N

| field | roots | N=5 | N=7 | N=9 | N=11 | N=13 | N=15 |
|---|---|---|---|---|---|---|---|
| poly5 | 5 | 6.7e-16 | 1.3e-15 | 1.0e-15 | 1.4e-15 | 1.3e-15 | 2.1e-15 |
| tangent | 3 | 1.1e-08 | 2.6e-09 | 0.0 | 2.4e-08 | 7.8e-16 | 2.0e-08 |
| gauss | 2 | 5.1e-02 | 2.4e-03 | 1.6e-03 | 2.1e-04 | 7.2e-06 | 2.6e-06 |
| tanh | 1 | 1.7e-02 | 2.8e-02 | 2.6e-02 | 1.3e-02 | 4.5e-03 | 7.7e-04 |
| runge | 2 | 1.7e-01 | 1.2e-01 | 7.4e-02 | 3.7e-02 | 1.3e-02 | 2.8e-04 |
| sine2 | 4 | 1.4e-01 | 1.5e-02 | 1.9e-03 | 2.4e-04 | 7.2e-06 | 2.8e-07 |
| sine3 | 6 | 1.0e+00 | 7.2e-01 | 5.5e-02 | 9.3e-03 | 2.9e-04 | 1.8e-04 |
| sine5 | 10 | 1.3e+00 | 6.0e-01 | 8.2e-01 | 8.4e-01 | 1.3e+00 | 3.8e-02 |

Identical to route A on every field **except `poly5`**, where route B holds ~1e-15 at every order
while route A degrades to 8.7e-11 at N=7 (and would fail outright at exactly `a_N = 0`). A
degree-5 polynomial sampled at N > 5 *is* degree-deficient — the case route B is built for.
Everything else sits exactly on the conditioning floor (see `../demo_A/README.md` for the floor
table and `../notes/05` for what it means).

### Degeneracy — the reason route B exists

`--study degeneracy`, N=7, driving the top Legendre mode to zero:

| `a_N` | route A error | route B error |
|---|---|---|
| 1e-2 | 2.7e-03 | 2.7e-03 |
| 1e-6 | 2.7e-07 | 2.7e-07 |
| 1e-10 | 4.7e-09 | 1.4e-09 |
| 1e-12 | 7.4e-07 | **1.4e-09** |
| 1e-14 | 1.1e-04 | **1.4e-09** |
| 0 | 3.4e-04 | **1.4e-09** |

Route B is flat across 12 orders of magnitude of `a_N`. Degree deficiency surfaces as an **extra
infinite eigenvalue** (`n_infinite` rises above its floor of 1) — information, not a failure. A
well-resolved SEM element is *precisely* the small-`a_N` regime, so this is the normal case.

### Why GLL specifically

The barycentric form is only as good as the conditioning of `b_j`. Spread `max|b|/min|b|`:

| N | GLL | equispaced |
|---|---|---|
| 8 | 3.66 | 7.0e+01 |
| 16 | 5.09 | 1.3e+04 |
| 32 | 7.15 | 6.0e+08 |
| 64 | 10.07 | 1.8e+18 |

`O(N)` on GLL vs `O(4^N/N)` equispaced. **Do not reuse this construction on equispaced points.**

## Cross-refs
`../notes/02` (full derivation, scale caveat, node-root edge case) · `../notes/04` (complexity) ·
`../notes/05` (the conditioning floor) · `../plans/A_gll_rootfinder.md` · ATAP p.139
