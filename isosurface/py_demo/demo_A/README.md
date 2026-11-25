# demo_A — Route A: modal Legendre comrade matrix **+ degree deflation**

Find **all** real roots of `f(x) − c` on a **single GLL reference element** `[-1,1]` of order `N`,
by transforming the nodal values to a Legendre (modal) expansion and computing the eigenvalues of
its **comrade matrix**. This is Trefethen *ATAP* Thm 18.1 (colleague, Chebyshev) transplanted to
Legendre, i.e. ATAP Exercise 18.4.

> **Degree deflation is ON by default** (`a_floor = 1e-12`). Plain route A divides by the top
> Legendre coefficient `a_N` and is *silently wrong* when the polynomial's true degree is below
> `N`; deflation trims trailing negligible coefficients first, which removes that failure entirely
> — **0 disagreements with route B over 20,000 randomized draws** — while keeping a **standard**
> eigenproblem. Run `--plain` to see the original behaviour. This is the form recommended as the
> GPU kernel in `../notes/08`.

Source mesh is always one element; `N` is the user's knob (tested 3–15, SEM-typical 5–10).

## Usage

```bash
python3 run_demo.py                      # every study, orders 3..15, writes figures
python3 run_demo.py --study convergence  # one study: convergence|cost|degeneracy|reference
python3 run_demo.py --orders 5,7,9,11    # custom order list
python3 run_demo.py --no-fig             # skip figure generation
python3 run_demo.py --plain              # DISABLE deflation (original, fragile)
python3 run_demo.py --a-floor 1e-10      # custom deflation threshold
python3 run_demo.py --study deflation    # threshold study: plain vs deflated vs B

python3 run_demo.py -N 8                             # single-order detail dump
python3 run_demo.py -N 8 --field sine3               # pick a test field
python3 run_demo.py -N 8 --field gauss --level 0.7   # override the isolevel
```

Fields: `poly5 tangent gauss tanh runge sine2 sine3 sine5` (see `../figs/00_testfields.png`).
The `-N` dump prints the GLL nodes/weights, every Legendre coefficient `a_k` with a magnitude bar,
the flagged ratio `|a_N|/max|a|`, the comrade matrix itself, `cond(C)`, the full complex spectrum,
and the found-vs-exact comparison.

As a library:

```python
from route_a import roots_route_a
r = roots_route_a(N, f)                          # deflation on (a_floor = 1e-12)
r, info = roots_route_a(N, f, return_info=True)  # info: degree, deflated, a_N_ratio
r = roots_route_a(N, f, a_floor=0.0)             # plain route A, for comparison
```

`../src/gllbase.py` holds the GLL basis, the test fields, and the error metric;
`../src/studies.py` holds the studies and plotting, shared with `demo_B`.

## The math

**Step 1 — nodal → modal.** Expand the interpolant in Legendre polynomials,
`p(x) = Σ_{k=0}^{N} a_k P_k(x)`, via the discrete projection

```
a_k = (1/γ_k) Σ_j w_j f_j P_k(x_j),        γ_k = <P_k,P_k>_GLL
```

GLL quadrature is exact to degree `2N-1`, so `γ_k = 2/(2k+1)` for `k < N` but

```
γ_N = 2/N        NOT 2/(2N+1)
```

because `P_N²` has degree `2N`. This top-mode inflation is the single easiest thing to get wrong,
and it matters here more than usual: `a_N` is the number step 3 divides by.

**Step 2 — multiply by x.** The Sturm–Liouville three-term recurrence
`(k+1)P_{k+1} = (2k+1) x P_k − k P_{k−1}` rearranges to

```
x P_k = k/(2k+1) · P_{k−1} + (k+1)/(2k+1) · P_{k+1},        x P_0 = P_1
```

so on `span{P_0..P_{N−1}}` the operator "multiply by x" is **tridiagonal with zero diagonal**.

**Step 3 — close the system.** Row `N−1` needs `P_N`, outside the span. At a root,
`P_N = −(1/a_N) Σ_{k<N} a_k P_k`, giving the spill-over row

```
C[N−1, :] −= (N/(2N−1)) · a[0:N] / a_N
```

Eigenvalues of `C` are exactly the roots of `p`. The isolevel enters once, as `f_j − c`.

**Step 2.5 — degree deflation (the fix).** Before building `C`, drop trailing coefficients that
are negligible:

```
keep a[0:d+1]   where   d = max{ k : |a_k| > a_floor · max_j |a_j| }
```

and build the `d × d` comrade matrix. The division is then by a *genuine* leading coefficient
rather than by roundoff.

**Why it is needed — and what the real failure is.** Not simply "`a_N` is small": a well-resolved
element has a small `a_N` and plain route A handles it fine (0 discrepancies in 4000 randomized
trials with realistic spectral decay, `../notes/07`). The failure is `a_N` small **while the low
modes stay O(1)** — a polynomial whose true degree is below `N` forced through a degree-`N`
comrade matrix. That happens for real (a degree-5 field sampled at N=6) and it fails *silently*:
at N=6 on `poly5`, plain route A returns five plausible real roots in `[-1,1]` with max error
**4.2e-01**, while `cond(C) = 2.4e17` and the last row carries entries ~1e15.

**Threshold choice** (`scratch/16`, 20,000 draws with 25% forced degree-deficient; "mismatch" =
disagreeing with route B on the *number* of roots):

| `a_floor` | mismatches | max \|A−B\| | deflation fired |
|---|---|---|---|
| 0 (off) | 24 | 2.8e-01 | 0.0% |
| 1e-16 | 17 | 2.8e-01 | 8.3% |
| 1e-14 | 0 | 3.1e-04 | 54.1% |
| **1e-12 (default)** | **0** | **2.8e-13** | **61.1%** |
| 1e-10 | 0 | 2.6e-11 | 68.0% |
| 1e-08 | 9 | 7.5e-09 | 75.0% |
| 1e-06 | 792 | 1.1e-06 | 82.1% |
| 1e-02 | 6810 | 1.8e-02 | 96.1% |

`1e-12` and `1e-10` are both in the safe plateau; `1e-12` is the tighter. Above `1e-8` the
threshold starts deleting genuine leading coefficients.

**What deflation does not fix.** A double root is √eps-conditioned for every method
(`../notes/05`). And on GPU the deflated size is data-dependent, so batches must be bucketed by
degree — measured distributions say this is rare (at N=8 with realistic decay the degree is full
~100% of the time). See `../notes/08`.

## Eigenproblem characterisation

| property | value |
|---|---|
| problem type | **standard**, non-symmetric: `C v = λ v` |
| size | `N × N` |
| structure | tridiagonal (zero diagonal) + **full last row**; a rank-one update of a Jacobi-like matrix |
| nonzeros | `nnz = 3N − 3` exactly (measured): 12/25 at N=5, 21/64 at N=8, 42/225 at N=15 |
| symmetry | **none** — the tridiagonal part is symmetrizable, the last row destroys it |
| LAPACK path | `scipy.linalg.eig` → `dgeev` = `dgebal` + `dgehrd` + `dhseqr` |
| balancing | **essential** (ATAP p.139, Parlett–Reinsch); never pass `nobalance` |
| flop model | `dgehrd` (10/3)N³ + `dhseqr` ~10N³ ≈ **13N³**, plus `2(N+1)²` for the transform |
| with deflation | same, at the deflated size `d ≤ N`; deflation itself is one `O(N)` pass |

**Note:** `Cᵀ` (coefficients in the last *column*) has the same eigenvalues and is **already upper
Hessenberg**, which would delete the (10/3)N³ reduction — ATAP Ex. 18.2. scipy exposes no
`dhseqr` binding to exploit it; measured gain is only 1.05–1.34×. See `../notes/04`.

## Complexity

`O(N³)` *asymptotically*, but **that is not the measured cost at SEM orders** — see `../notes/04`.
Fitting `t = c·N^p` over N=3…15 gives **p = 0.39**: the solve is a few thousand flops buried under
LAPACK and python call overhead. Breakdown at N=8: nodal→modal 0.64 µs, build `C` 3.15 µs,
`eig(C)` 13.33 µs, total 22.6 µs. True cubic behaviour only appears near N≈80.

## Results

Baseline — what is being solved: `../figs/00_testfields.png`, `../figs/00_testfield_resolution.png`
Convergence: `../figs/A_convergence.png`  ·  Cost: `../figs/A_cost.png`
Full console output: `output.txt` (deflation on) and `output_plain.txt` (`--plain`, for contrast)

### Max root error vs N (deflation ON — the default)

| field | roots | N=5 | N=7 | N=9 | N=11 | N=13 | N=15 |
|---|---|---|---|---|---|---|---|
| poly5 | 5 | 1.0e-15 | 1.4e-15 | 2.1e-15 | 2.2e-15 | 1.4e-15 | 2.1e-15 |
| tangent | 3 | 1.7e-08 | 4.4e-16 | 4.5e-08 | 1.5e-08 | 3.9e-08 | 2.2e-15 |
| gauss | 2 | 5.1e-02 | 2.4e-03 | 1.6e-03 | 2.1e-04 | 7.2e-06 | 2.6e-06 |
| tanh | 1 | 1.7e-02 | 2.8e-02 | 2.6e-02 | 1.3e-02 | 4.5e-03 | 7.7e-04 |
| runge | 2 | 1.7e-01 | 1.2e-01 | 7.4e-02 | 3.7e-02 | 1.3e-02 | 2.8e-04 |
| sine2 | 4 | 1.4e-01 | 1.5e-02 | 1.9e-03 | 2.4e-04 | 7.2e-06 | 2.8e-07 |
| sine3 | 6 | 1.0e+00 | 7.2e-01 | 5.5e-02 | 9.3e-03 | 2.9e-04 | 1.8e-04 |
| sine5 | 10 | 1.3e+00 | 6.0e-01 | 8.2e-01 | 8.4e-01 | 1.3e+00 | 3.8e-02 |

Compare against the **root-displacement floor** `|δp(x*)|/|p′(x*)|` — the conditioning of the
problem, i.e. the error a *perfect* rootfinder would still make given the same interpolant:

| field | N=5 | N=7 | N=9 | N=11 | N=13 | N=15 |
|---|---|---|---|---|---|---|
| poly5 | 6.6e-16 | 9.9e-16 | 1.7e-15 | 1.6e-15 | 2.1e-15 | 3.2e-15 |
| gauss | 5.2e-02 | 2.4e-03 | 1.6e-03 | 2.1e-04 | 7.2e-06 | 2.6e-06 |
| tanh | 1.8e-02 | 2.8e-02 | 2.6e-02 | 1.3e-02 | 4.5e-03 | 7.7e-04 |
| runge | 1.5e-01 | 1.6e-01 | 8.2e-02 | 3.8e-02 | 1.3e-02 | 2.8e-04 |
| sine2 | 1.6e-01 | 1.4e-02 | 1.9e-03 | 2.4e-04 | 7.2e-06 | 2.8e-07 |
| sine3 | 5.5e-01 | 6.0e-01 | 7.3e-02 | 9.6e-03 | 2.9e-04 | 1.8e-04 |

**The two tables agree digit-for-digit** on every non-degenerate field: route A is exact, and all
remaining error is interpolation. (Confirmed independently: roots match 40-digit mpmath roots of
the *same interpolant* to 4e-16 — `../notes/05`.) `tangent` is a double root and is intrinsically
√eps-accurate for every method (`../notes/05`).

`poly5` is the row to watch: with deflation it holds **~1e-15 at every order**. Without it
(`--plain`, saved in `output_plain.txt`) it reads 8.7e-11 at N=7 and **4.2e-01 at N=6**.

### Deflation: before and after

`--study deflation` (excerpt; `*` = deflation materially improved the answer):

| field | N | `\|a_N\|/max` | deg | fired | plain A | deflated | route B |
|---|---|---|---|---|---|---|---|
| poly5 | 5 | 1.00e+00 | 5 | – | 9.99e-16 | 9.99e-16 | 6.66e-16 |
| poly5 | 6 | 8.20e-17 | 5 | yes | **4.20e-01** | \* **1.11e-15** | 1.11e-15 |
| poly5 | 8 | 1.42e-15 | 5 | yes | 1.85e-12 | \* 5.55e-16 | 8.88e-16 |
| poly5 | 11 | 2.21e-15 | 5 | yes | 6.99e-15 | 2.22e-15 | 1.44e-15 |
| gauss | 8 | 3.05e-02 | 8 | – | 4.60e-03 | 4.60e-03 | 4.60e-03 |
| gauss | 15 | 1.18e-04 | 15 | – | 2.65e-06 | 2.65e-06 | 2.65e-06 |

Deflation fires only where the degree genuinely is deficient, and is a no-op elsewhere.

Across **91 non-tangency (field, N) cases**: plain A disagrees with route B on **2**, deflated A
on **0**. Across 20,000 randomized draws: **0 root-count mismatches, max diff 3.1e-11**.

### Cost

| N | A µs | B µs | B/A | A nnz/dense | B nnz/dense |
|---|---|---|---|---|---|
| 3 | 17.2 | 26.4 | 1.54× | 6/9 | 12/25 |
| 5 | 20.7 | 28.8 | 1.39× | 12/25 | 18/49 |
| 7 | 22.4 | 31.6 | 1.41× | 18/49 | 24/81 |
| 9 | 25.4 | 34.4 | 1.36× | 24/81 | 30/121 |
| 11 | 28.0 | 37.4 | 1.34× | 30/121 | 36/169 |
| 13 | 32.3 | 43.1 | 1.33× | 36/169 | 42/225 |
| 15 | 38.7 | 47.0 | 1.22× | 42/225 | 48/289 |

Route A is the cheaper route, but by **1.2–1.5×**, not the ~4–8× the flop model predicts — both
are overhead-bound.

### The degeneracy sweep, with and without deflation

`--study degeneracy`, N=7, driving the top mode to zero; error vs the exact roots of the
degree-6 limit:

| `a_N` | plain A (`--plain`) | **A + deflation** | route B |
|---|---|---|---|
| 1e-2 | 2.7e-03 | 2.7e-03 | 2.7e-03 |
| 1e-6 | 2.7e-07 | 2.7e-07 | 2.7e-07 |
| 1e-10 | 4.7e-09 | 4.7e-09 | 1.4e-09 |
| 1e-12 | 7.4e-07 | 7.4e-07 | 1.4e-09 |
| 1e-14 | 1.1e-04 | **1.4e-09** | 1.4e-09 |
| 0 | **3.4e-04** | **1.4e-09** | 1.4e-09 |

Deflation flattens the tail: once `|a_N|/max|a|` drops below `1e-12` the degree is trimmed and the
answer matches route B exactly. Between `1e-12` and `1e-10` the sweep still shows the plain
behaviour (the threshold has not yet fired) — tighten `--a-floor` if that band matters, at the
cost of firing more often.

## Cross-refs
`../notes/01` (structure) · `../notes/03` (parity trap) · `../notes/04` (why N³ is not the measured
cost) · `../notes/05` (the conditioning floor) · `../notes/07` (A vs B, and when A really fails) ·
`../notes/08` (3D/GPU: this deflated form is the recommended kernel) ·
`../plans/A_gll_rootfinder.md` · ATAP ch. 18 Thm 18.1, Ex. 18.2/18.3/18.4
