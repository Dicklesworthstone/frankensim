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
- `det::TRIG_DOMAIN` = 2²⁰: trig budgets valid within |x| ≤ 2²⁰ (4-part
  Cody–Waite reduction); beyond → deterministic but budget-void (no-claim).
- Policy vocabulary: `canonical_nan`, `next_up/next_down`, `nudge_out`
  (fs-ivl's directed-rounding primitive), `ulp_distance`.
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
- tan/atan2/pow/cbrt/log1p/erf: not yet implemented (follow-up bead).
- Trig beyond |x| > 2²⁰ (Payne–Hanek is recorded follow-up scope).
- Correctly-rounded (0.5 ULP) results: NOT claimed — budgets above.
- dd-oracle billions-scale nightly battery arrives with fs-ivl.
