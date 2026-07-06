# CONTRACT: fs-ivl

## Purpose and layer
Certified arithmetic (plan §6.4): outward-rounded interval arithmetic and
affine (noise-symbol) forms. This crate defines what "certified" means
system-wide: every operation's postcondition is ENCLOSURE — the result
contains the true image of the inputs. Layer: **L1**. Depends only on
fs-math (strict functions, declared budgets, nudging primitives; the
high-precision oracle rungs `eft`/`dd` live THERE at L0, single
implementation shared with fs-la — recorded relocation, beads
6ys.8/6ys.12).

## Public types and semantics
- `Interval` — closed [lo, hi], ±∞ endpoints allowed, NaN rejected at
  construction, INVARIANT lo ≤ hi. Std operator traits (+, −, ×, ÷, unary −)
  plus `sqrt/exp/ln/tanh/sin/cos`, `abs`, set ops (`hull`, `intersect`,
  `encloses`, `contains`), `width` (rounded up — never understates),
  `midpoint` (representative only, NOT certified central).
- Directed rounding: computed endpoints are nudged outward via
  fs-math `next_up`/`next_down` — basic ops 1 ULP (rigorous: IEEE
  correctly-rounded ops err ≤ ½ ULP), elementary functions by fs-math's
  DECLARED budgets (exp/expm1/ln/sin/cos 3, tanh 5, sqrt 1 after the exact
  hardware sqrt). NO global rounding-mode state anywhere (thread-safe,
  SIMD-mixing-safe, grep-lintable).
- `Interval::WHOLE` = [−∞, +∞]: the "no useful enclosure" answer. Division
  by a zero-containing interval returns WHOLE (documented; not a panic —
  certified callers branch on it). 0·∞ ambiguities in mul/div also return
  WHOLE (conservative, never wrong).
- sin/cos: non-monotone handling by conservative critical-point counting
  against an interval enclosure of π; "possibly contains a critical point"
  counts as containment. Width ≥ 2π (or beyond `det::TRIG_DOMAIN`, or
  unbounded) → [−1, 1] (honest, always correct).
- `AffineCtx`/`Affine` — affine forms `c₀ + Σ cᵢ·εᵢ + [−err, +err]`.
  Symbol identity IS correlation (same-context symbols cancel; fresh
  symbols don't). All rounding and the first-order mul residue
  (rad·rad) are absorbed into `err` with upward rounding. Reference
  operator impls (`&x + &y`), `scale`, `to_interval` (rigorous collapse),
  `radius`. Symbol ids are deterministic (context counter — replay-stable).

- `expansion` — Shewchuk floating-point expansion arithmetic (exact sums,
  scaled products, expansion×expansion products, exact-diff two-component
  constructors, sign extraction). Every op is error-free: the output's
  exact sum equals the exact real result (residual-law tested).
- `predicates` — adaptive-precision EXACT `orient2d`/`orient3d`/`incircle`/
  `insphere` (+ `*_with_stage` telemetry variants): fast float evaluation
  guarded by Shewchuk's proven stage-A error bounds (arrangement-specific
  constants with ε = 2⁻⁵³), escalating to exact expansion evaluation of
  the standard difference determinants (orient2d runs the full faithful
  B/C/D adaptive ladder). `orient2d_sos` implements the Edelsbrunner–Mücke
  tie-breaking ladder (total, antisymmetric, deterministic);
  `orient3d_sos` is a deterministic projection cascade.

## Invariants
1. **Containment law (G0)**: for every op and every point tuple inside the
   inputs, the true result lies inside the output. Tested: 100k arithmetic
   checks vs the dd oracle, 20k rewrite pairs (x·(y+z) vs x·y+x·z both
   contain truth), 50k intervals × 6 elementary functions vs an independent
   platform-libm point oracle, 20k point checks over 4000 random expression
   DAGs, Rump's polynomial containing the exact −54767/66192 (wide-and-
   right beats narrow-and-wrong — asserted width > 1 documents the honest
   outcome).
2. Affine linear terms cancel EXACTLY on shared symbols: x − x collapses to
   width < 1e−13 (tightness ratio vs plain intervals > 1e10, tested);
   first-order identity-zero DAGs stay O(radius²) wide, independent of
   center magnitude.
3. Trig enclosures capture extrema: intervals straddling π/2 + kπ (sin) or
   kπ (cos) include ±1 exactly.
4. Enclosures never claim outside function ranges (sin/cos ⊆ [−1,1],
   tanh ⊆ [−1,1], exp ≥ 0, sqrt ≥ 0).

## Error model
Domain violations that admit NO enclosure panic with structured messages:
NaN/inverted constructor endpoints, sqrt of entirely-negative interval,
ln with hi ≤ 0. Partial-domain overlaps degrade gracefully (sqrt clips at
0; ln returns lower bound −∞). Division by zero-containing intervals
returns `WHOLE`, never panics.

## Determinism class
Bit-deterministic CROSS-ISA by construction (straight-line IEEE arithmetic
+ fs-math strict functions). Evidence: FNV-64 golden hash over 500 random
DAG enclosure endpoints + affine collapses = `0x3712_a4c1_2d5e_5864`,
recorded on aarch64-apple (M4 Pro), required to match on x86-64
(Threadripper) in tests/conformance.rs. Golden-evidence policy applies.

## Cancellation behavior
Straight-line arithmetic; no poll points needed (certified BULK kernels
belong to future consumers, which poll at their own tile boundaries).

## Unsafe boundary
None. `unsafe_code` denied; no capsules.

## Feature flags
None.

## Conformance tests
`tests/conformance.rs`: random-DAG containment battery + golden hash.
`tests/predicates.rs`: adversarial degeneracy batteries (cocircular /
cospherical / collinear lattice configurations), dyadic and 1-ulp
perturbation classification, SoS tie determinism, measured stage-A filter
rates. In-crate predicate batteries: i128 lattice oracle agreement (2D and
3D, thousands of rounds with degeneracy hits asserted), the Kettner-class
half-ulp grid where the naive float determinant demonstrably misjudges
while the exact ladder matches integer ground truth, expansion residual
laws, and SoS totality/antisymmetry.
In-crate: dd-oracle arithmetic containment, rewrite pairs, elementary
containment vs libm oracle, trig extrema capture, zero-divisor semantics,
Rump's polynomial, set operations, affine cancellation/tightness/
independence, structured rejections. Reimplementations must pass the
golden-hash case bit-for-bit.

## No-claim boundaries
- **Rigor is conditional on fs-math's declared ULP budgets** (empirically
  enforced there): if a budget were violated, elementary enclosures could
  under-cover by the violation amount. Basic-op enclosures are
  unconditionally rigorous.
- No Taylor models yet (future bead 6ys.13).
- Predicates are certified for inputs whose difference monomials (degree
  ≤ 5) stay inside the normal f64 range — no intermediate overflow or
  underflow (Shewchuk's standard caveat, inherited; out-of-range inputs
  are outside the certificate).
- `incircle`/`insphere` symbolic perturbation is NOT provided (hooks and
  2D/3D orientation SoS only); the full 3D Edelsbrunner–Mücke ladder for
  `orient3d_sos` is replaced by a documented deterministic projection
  cascade — revisit when fs-mesh's Delaunay needs E-M-exact semantics.
- Predicate throughput: only stage-A filter rates are measured; the exact
  tail is Vec-based and deliberately unoptimized (cold path).
- No quad-double (dd suffices for current oracles; recorded on 6ys.12).
- No SIMD lanes: scalar everywhere v1 (correctness-over-lanes, per plan).
- Affine elementary functions (exp/sin on affine forms via Chebyshev
  linearization) are not implemented — interval fallback applies.
- `Interval::midpoint` of half-infinite intervals returns 0.0 (documented
  representative-point semantics, not an enclosure claim).
