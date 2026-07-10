# CONTRACT: fs-math

## Purpose and layer
Deterministic elementary functions (strict mode) and the workspace
floating-point POLICY: FMA contraction, subnormals, NaN, ULP budgets
(patch Rev O; plan §5.4/§6.4). Layer: L0.

## Public types and semantics
- `det::{exp, expm1, ln, sin, cos, tanh, sqrt}` — strict-mode functions
  built EXCLUSIVELY from IEEE arithmetic (+,−,×,÷, mul_add, sqrt): bit-
  identical cross-ISA BY CONSTRUCTION, empirically PROVEN (golden hash
  0xeb79cab7a01643e5 identical on aarch64-apple M4 Pro and x86-64 TR 5995WX).
- Declared ULP budgets (measured maxima in parentheses, vs platform-libm
  oracle, 200k samples + edges): exp 3 (1), expm1 3 (2), ln 3 (1),
  sin 3 (2), cos 3 (2), tanh 5 (3). sqrt is 0 ULP (IEEE-correctly-rounded
  hardware).
- EXTENSION FAMILY (bead wf9.14, additive): `det::{tan, atan, atan2,
  erf, erfc, pow}`. Declared budgets (measured, 200k+ samples + edges):
  tan 8 (shared Cody–Waite reduction; tan = sin/cos BITWISE on even
  quadrants — an identity, not an approximation), atan/atan2 4 (+1 for
  atan2 composition), erf 6, erfc 10 (deep tail included: exact-dd x²
  before exp — plain f64 there would cost ~x²/2 ULP), pow = the HONEST
  formula 3·(|y·ln x|+1)+5 (the exp∘ln magnification is intrinsic;
  dd-ln refinement recorded). erf/erfc run their cancellation-prone sums
  in double-double (why dd lives at L0). Oddness of tan/atan/erf and
  atan2's y-sign symmetry hold BITWISE by construction. In the erf
  band x ∈ (1.5, 3.5) the external test oracle is weaker than the
  implementation; budget-grade evidence there is DISJOINT-PATH
  cross-validation (Taylor-dd vs CF-dd agree ≤ 3 ULP, tested).
- PINNED-ORDER INTEGER POWERS (bead 4xnt, additive): `det::powi(x, n)` —
  LSB-first binary exponentiation with a fixed source-level operation tree
  (negative powers use one final reciprocal, with a reciprocal-base
  range-recovery pass only when an overflowed intermediate would erase a
  representable subnormal; n = 0 → 1.0 for every x; i32::MIN handled).
  Exists because `f64::powi`'s rounding
  is optimization-level-dependent (llvm.powi has no pinned order;
  observed 1-ULP debug/release divergence from n = 4 up), which is a
  build-mode determinism hazard in any golden-feeding path. Positive
  n ≤ 512 agrees BITWISE with `pow`'s integer fast path (same order,
  every exponent tested). NOT correctly
  rounded: one rounding per executed multiply; measured ≤ 2 ULP vs
  platform powi for |n| ≤ 64. Own golden hash, identical in both build
  modes by construction.
- INVERSE-TRIG COMPLETION (bead t88x, additive): `det::{asin, acos}` via
  atan2 on the FACTORED complement √((1−x)(1+x)) (endpoint-conditioned;
  1 − x² cancels catastrophically at |x| → 1). Declared budget 6
  (measured worst 3, 200k samples). asin odd BITWISE (atan2's sign fold
  + commuting factors); acos(±1) = {0, π} and asin(±1) = ±π/2 exact
  through atan2's special table; |x| > 1 → NaN. Identity checks measure
  at the IDENTITY's scale (π − acos re-measures at the small result's
  scale and inflates conditioning ~16× — measured, documented in the
  battery).
- `det::TRIG_DOMAIN` = 2²⁰: the Cody–Waite/Payne–Hanek dispatch boundary
  (4-part Cody–Waite at and below; the `payne` module's 1280-bit reduction
  above). Trig budgets now hold for ALL finite arguments: declared
  `SIN_LARGE_ULP_BUDGET` = 4 beyond the boundary, measured max 1 ULP over a
  4000-sample exponent sweep 2²¹..2¹⁰⁰⁰ against the platform-libm oracle,
  0 ULP on the published worst-case double 6381956970095103·2⁷⁹⁷ (reduced
  |r| = 4.7e-19). Odd/even symmetry BITWISE at every landmark.
- `payne`: SELF-VERIFYING constants — the 2/π limbs are hardcoded AND
  regenerated at test time by an all-integer Machin bignum (π = 16·atan(1/5)
  − 4·atan(1/239) in u64-limb fixed point, binary long division for 2/π);
  the regeneration test compares every limb and the π hex expansion against
  published digits (G5: integer arithmetic, bit-identical on every ISA).
  The `Fx` bignum doubles as reference machinery for hard-case tests.
- Policy vocabulary: `canonical_nan`, `next_up/next_down`, `nudge_out`
  (fs-ivl's directed-rounding primitive), `ulp_distance`.
- `c64::C64` — complex f64 (bead urvw): operator traits with strict
  arithmetic, overflow-safe magnitude (max-scaled, no libm hypot —
  tested at 1e±300), Smith division/reciprocal (scaling-robust),
  principal sqrt via stable half-angle formulas (both half-planes
  tested). The shared complex home going forward; fs-fft's private
  mini-type migration is recorded cleanup.
- `eft::{two_sum, quick_two_sum, two_prod}` — error-free transformations:
  the returned (result, error) pair reconstructs the EXACT real value
  (bitwise-testable identities; `quick_two_sum` requires |a| ≥ |b|,
  debug-asserted). Relocated here from fs-la's mixed-precision scope so
  fs-ivl and fs-la share one implementation (beads 6ys.8/6ys.12).
- `dd::Dd` — double-double (~106-bit significand) via std operator traits
  (+, −, ×, ÷) plus `abs/sqrt/lt`. Documented error bounds: add/sub/mul
  ≤ 2⁻¹⁰⁴ relative, div/sqrt ≤ 2⁻¹⁰³, finite non-over/underflowing
  operands. Normalization invariant `hi = fl(hi+lo)` property-tested.
  Quad-double is recorded follow-up scope (not needed by current oracles).

## Invariants
- No platform libm on any strict path (sqrt excepted: IEEE-exact).
- Reduction constants are EXACT bit patterns with trailing-zero mantissas
  (j·part products exact) — decimal literals are forbidden there (a 184-ULP
  bug class, regression-tested).
- tanh/sin odd and cos even BITWISE (symmetry by construction).
- exp(0)=1, ln(1)=+0, sin(0)=0, cos(0)=1 exactly; NaN in → NaN out;
  subnormals never flushed.
- Golden hash changes require a schema-bump-style justification.

## Error model
Total functions; domain violations return NaN/±inf per IEEE conventions.

## Determinism class
Deterministic CROSS-ISA (the strongest class in the workspace) — proven.

## Cancellation behavior
Straight-line arithmetic; no poll points needed.

## Unsafe boundary
None.

## Feature flags
None (fast-mode platform-libm variants are recorded follow-up scope).

## Conformance tests
Per-function ULP batteries (budget-gated, measured maxima printed as JSONL),
tiny-x expm1 cancellation battery, near-1 ln battery, bitwise symmetry
sweeps, special-value policy table, nudge bracketing, cross-ISA golden hash,
core-only + worst-case-point + constant-integrity regressions
(tests/core_regression.rs). All verified on BOTH reference ISAs.

## No-claim boundaries
- cbrt/log1p: not yet implemented (follow-up bead). (tan/atan/atan2/pow/
  erf/erfc/asin/acos/powi landed via wf9.14, t88x, and 4xnt — see above;
  this line previously understated the implemented surface.)
- `det::powi`: f32 variant not provided; no claim that its bits match
  platform `f64::powi` (the source-level tree is intentionally pinned).
- Trig beyond |x| > 2²⁰: RESOLVED (bead r6r5, `payne` module — see above).
- Fast mode (lower-accuracy feature-flagged variants): BLOCKED on consumers
  declaring tolerable budgets (fs-material/LUMEN) — deliberately not built
  speculatively, per the bead's own instruction.
- The nightly ULP-ledger re-measurement lane: the budget-vs-measured tests
  ship here and run in every suite; wiring them into a dedicated nightly
  regression lane belongs to the CI/CD bead family (huq.4 closed; the perf-CI
  bead fz2.4 owns nightly gates).
- Correctly-rounded (0.5 ULP) results: NOT claimed — budgets above.
- dd-oracle billions-scale nightly battery arrives with fs-ivl.
