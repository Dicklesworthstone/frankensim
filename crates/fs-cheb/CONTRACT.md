# CONTRACT: fs-cheb

> Status: PARTIAL — the 1D core, collocation, ORR–SOMMERFELD,
> 2D low-rank, Fourier-periodic, and colleague-root sections are in
> force; 3D low-rank and Qty integration remain follow-up scope.

## Purpose and layer
Chebfun-style function objects (plan §6.5): smooth 1D functions as
adaptively truncated Chebyshev expansions, plus spectral collocation
differentiation matrices. Layer: **L1**. Deps: fs-fft (DCT/FFT pair),
fs-la (LU/eigen paths), fs-math (strict elementary functions),
fs-ivl (interval root certification), and fs-exec (`Cx` for the
budgeted, cancellable entry points).

## Budgets and admission (`budget` module, bead sj31i.55 slice 1)
- `ChebBudget` (schema `CHEB_BUDGET_SCHEMA_VERSION = 2`, non-exhaustive,
  explicit caps: retained coefficients, total adaptive samples,
  collocation dimension, abstract work ops, peak temporary bytes) plus
  `admit_adaptive_build` / `admit_dirichlet_eigs` / `admit_root_scan`:
  conservative worst-case samples/coefficients/work/temporary-byte
  formulas run with CHECKED `u128` arithmetic BEFORE allocation or
  function evaluation. The adaptive envelope includes sampled values,
  DCT-II complex data/scratch, twiddles, output, and six-step headroom;
  root admission includes normalization, derivative recurrence/output,
  every-cell refinement, and retained candidates; eigensolve admission
  includes persistent matrices, cyclic-Jacobi copies, blocked-LU update
  and GEMM-pack workspace, iteration vectors, and result storage.
  huge requests refuse as typed `ChebError::CapExceeded` or
  `ChebError::Overflow` — never a saturated size that still iterates,
  panics, or allocates. `ChebAdmission` (sealed) is the evidence the
  preflight ran.
- Budgeted, cancellable entry points thread an explicit `Cx` and poll
  at bounded boundaries: `try_build_budgeted` (per adaptive round;
  bitwise-identical to `Cheb1::build` on the happy path; `Cancelled`
  carries a `resume_from` grid whose resumption is deterministic and
  bitwise-equivalent to an uncancelled run),
  `dirichlet_laplace_eigs_budgeted` (before/after opaque matrix kernels,
  per shift, and per 10 inverse-power sweeps; cancellation retains only
  a prefix of completed fixed-sweep estimates, with NO convergence or
  residual certificate), and `Cheb1::roots_budgeted`
  (per 64 scan cells; cancellation returns NO partial — an incomplete
  scan is not a root-set claim). Terminal states are explicit:
  `Complete`/`Cancelled` run enums with deterministic `WorkReceipt`s,
  typed `ChebError` refusals (`Domain`/`Shape`/`CapExceeded`/`Overflow`/
  `Unresolved`/`NonFinite`/`Numerical`/`Cancelled`) where the classic
  APIs panic. The budgeted eigensolve also refuses `k > 64` (the fixed
  FD surrogate supplies at most 64 shifts; the classic API silently
  shorts).

## Public types and semantics
- `Cheb1` — coefficients over FIRST-KIND Chebyshev points (roots grid):
  values ↔ coefficients is exactly fs-fft's DCT-II/III pair. `build`
  doubles the grid until the trailing quarter of coefficients sits at
  the configured numerical plateau (2.2e-15 relative), then truncates;
  unresolvable functions panic at `max_degree` with a structured
  message. `eval` (Clenshaw), `differentiate` (coefficient recurrence,
  domain chain rule), `integral` (even-coefficient formula), `add`,
  `mul` (resample + rebuild), `roots` (fixed reference-grid sign scan +
  safeguarded reference-coordinate bisection/Newton polish for isolated,
  well-conditioned sign-changing roots).
- `lobatto_points`, `diff_matrix` — Chebyshev–Lobatto collocation:
  Trefethen construction with the negative-sum-trick diagonal (rows sum
  to EXACT zero, tested bitwise).
- `orr_sommerfeld::{growth_rates, max_growth, critical_reynolds}` —
  plane-Poiseuille temporal stability via clamped Chebyshev collocation
  (Trefethen D4c construction from the φ = (1−y²)u substitution),
  generalized problem reduced through fs-la's complex LU, spectrum via
  the complex QR eigensolver. `growth_rates` is the "modal growth rates
  σ₁..σ_k" first-class query (descending real part, deterministic
  tie-breaks). ACCEPTANCE EVIDENCE: the neutral crossing at α = 1.02056
  reproduces the published Re_c = 5772.22 (displayed digits exact at
  N = 48); stability verdicts correct on both sides; golden hash
  `0x7b3b_e74e_d5a6_faad` cross-ISA.
- `dirichlet_laplace_eigs` — deflated inverse-power-iteration demo of
  the collocation eigen path (validates against analytic (kπ/2)²).

## Invariants
1. Machine-precision recovery on analytic fixtures with expected degree
   growth (exp ≤ 20, Runge in (exp, 300], sin(20x) on [0,3] in
   (40, 200]) — tested.
2. Calculus identities: d/dx exp = exp to 1e-11; definite integrals to
   1e-12 with domain scaling — tested.
   Affine domain maps preserve the established ordinary-domain rounding path
   while avoiding either an overflowing `b-a` or a finite-width doubled
   numerator overflow on extreme one-sided and same-sign domains.
   Center/radius evaluation retains representable offsets around the center;
   calculus scaling is combined before intermediate overflow/underflow;
   integral accumulation uses an error-free partial expansion so representable
   cancellation residuals survive even when every naive prefix is finite, then
   falls back to exact common power-of-two normalization when an expansion
   prefix would overflow. Tests cover constants,
   physical/reference linears, subnormal and scale-separated integrals, and
   roots across the finite endpoint range. Polynomial values or derivatives whose
   final f64 result is itself unrepresentable remain an explicit no-claim.
3. Plateau detection does NOT chase noise floors (tested with a
   deterministic ~1e-18 jitter fixture).
4. `diff_matrix` rows sum to exact zero (differentiation annihilates
   constants bitwise).
5. Deterministic per ISA: all state is built on strict fs-math cos/sin
   and fixed-order arithmetic. The last admitted radix-8/4/2 downstream
   golden was identical in debug and release on M4 Pro. Current root/integral
   semantics intentionally invalidate it; replacement and x86-64 equality are
   pending an admitted two-profile, two-ISA replay. The upstream FFT stage-path
   golden remains independently verified in all four ISA/profile quadrants.
6. Policy-filtered colleague candidates agree with the sign scanner on the
   admitted isolated fixtures, recover the retained even-multiplicity fixtures
   the scanner cannot see, and
   certification boxes are reported in physical-domain coordinates.
   Approximate colleague coefficients use exact power-of-two normalization and
   refuse any coefficient exponent range, matrix half-ratio, or recurrence-row
   addition that would lose a non-zero term before eigenanalysis. Certified roots likewise refuse
   coefficient information loss, perform the c0/2
   convention in interval arithmetic, and enclose the exact-real derivative
   through interval automatic differentiation. Affine
   images are outward-rounded with interval arithmetic and clamped to the
   finite physical domain; a widened physical box is revalidated over its
   full inverse image before retaining Certified existence-and-uniqueness
   authority. `certified_roots` currently interprets `min_width` in the
   dimensionless reference coordinate; a typed physical-width API is tracked
   separately.
7. `Cheb2` captures separable rank exactly on fixture functions, keeps
   deterministic pivot tie-breaking, and converges spectrally on the
   smooth non-separable fixture. All public components have one common x-domain,
   one common y-domain, and finite non-zero inverse pivots. Three-factor component products try the
   established order first and then safe pairings so a representable result is
   not silently lost to an overflowing or underflowing intermediate. Component
   sums use the same error-free expansion as 1D integration, and fixed-slice DCT
   terms apply `2/n` before accumulation so a representable coefficient is not
   rejected merely because its unscaled prefix would overflow.
8. `FourierSeries` exactly recovers trigonometric fixture modes,
   differentiates `sin` to `cos`, and uses c₀ for the periodic integral.

## Error model
Structured panics for programmer/modeling errors: non-finite or
inverted domains, non-finite samples/coefficients, unresolved functions
at `max_degree`, domain mismatches in algebra, invalid colleague
policies, an unrepresentable colleague normalization, matrix half-ratio, or
recurrence-row addition,
non-positive certification widths, non-identical algebra domains,
an unrepresentable algebra/transform coefficient or Cheb2 inverse pivot,
an exact root normalization (scanner or certifier) that would lose coefficient information,
a root query on the identically-zero polynomial (whose root set is a
continuum), a detected root candidate whose local slope the fixed-grid
fallback cannot resolve honestly, a physical
derivative whose coefficient representation is not finite `f64`, malformed
public `Cheb2` or `FourierSeries` fields, and non-power-of-two Fourier sample
counts.

## Determinism class
Bit-deterministic per ISA by construction. The golden hashes coefficients
+ integral + derivative sample + roots + collocation eigenvalues. The
last admitted radix-8/4/2 value was recorded in both build profiles on
aarch64-apple, but the current reference-coordinate root-refinement and
exceptional integral semantics intentionally invalidate that downstream
value. The test remains armed and known-red pending a current admitted
debug/release replay on both reference ISAs; no replacement value or G5 claim
is accepted from a stale binary. The coupling remains registered in
`golden-couplings.json` against `fs-fft:transform-bits=1` so later transform
changes cannot strand it.

## Cancellation behavior
The classic entry points are bounded (max_degree cap) with no poll
points — retained unchanged this slice for their existing callers. The
budgeted twins poll before allocation/evaluation, at adaptive sample
boundaries, before/after opaque transforms and dense kernels, per eigen
shift and 10-sweep tile, throughout root normalization/refinement, and
per 64 scan cells. A final poll precedes every complete result. They
drain to explicit `Cancelled` states: resumable for construction,
diagnostic-prefix-retaining for the eigensolve, and refusing without
partials for the root scan.

## Unsafe boundary
One registered capsule: `src/fma/mod.rs` (+ SAFETY.md beside it, entry
in `unsafe-capsules.json`) — the bead-nabk x86 FMA-codegen pattern.
The `unsafe` is ONLY the `target_feature(enable = "avx2,fma")` calling
contract around `#[inline(always)]` safe bodies (Clenshaw `eval`, the
Dirichlet D·D product and Rayleigh matvec, the Orr–Sommerfeld matmul);
runtime-detected, and the portable body IS the twin. Bit-identical by
construction (native fused instruction vs correctly-rounded libm
`fma()`; chain shapes untouched). NO performance claims are made for
these paths — the crate has no perf lane; the capsule exists to close
the baseline-x86 per-element-libm-call hazard class, not to certify a
number. The census sites deliberately NOT capsuled: `FourierSeries::
eval` (trig-call-bound), `Cheb1::differentiate` (once-per-call, alloc-
dominated), the `os_matrices` D4-clamp assembly (one O(n²) pass beside
three O(n³) matmuls).

## Feature flags
None.

## Conformance tests
tests/cheb_battery.rs (recovery, calculus, plateau robustness, roots,
collocation accuracy, eigen demo, golden hash).

tests/budget_battery.rs (bead sj31i.55, cases cb-001..cb-006): G0
admission boundary tables incl. huge requests refusing before
allocation and at-cap/one-over work and temporary envelopes; typed
domain/shape/root-normalization refusals;
bitwise parity between budgeted and classic construction/eigensolve/
root paths; real cancellation with deterministic resume equivalence
and empty-prefix refusal; typed `Unresolved`/`NonFinite` refusals where
the classic API panics; same-profile receipt determinism. Cross-ISA G5
evidence remains pending as stated above. Slice-2/3 cases
cb-007..cb-009: colleague admission boundaries + bitwise budgeted/
classic candidate parity + pre-eigen cancellation drain; admission
tables for Cheb2/Fourier/Orr-Sommerfeld incl. typed shape refusals and
huge requests; budgeted builder twins (Fourier bitwise mirror with
typed non-finite/shape refusals, Cheb2 and OS admission-gated wrappers
with boundary drains, solver failures mapped to `Numerical`).

## Variants (bead kw89)

- `colleague::colleague_roots` — the Chebyshev companion matrix
  (three-term-recurrence rows, coefficient-loaded last row scaled by
  −1/(2aₙ)), eigenvalues via the fs-la complex nonsymmetric stack,
  filtered by a DOCUMENTED [`ColleaguePolicy`] (trailing-coefficient
  trim, imaginary tolerance, domain slack, √ε-scale cluster dedupe —
  a double root's eigenvalue pair splits at ~5e-9, measured). This COVERS the
  retained even-multiplicity fixture that motivated the path: an
  (x−r)²(x−s) case the sign scanner misses is found under the retained policy
  (cheb-102). Close candidates can be clustered or filtered, so this API does
  not claim complete enumeration. Cheb1 stores
  the Σ′ convention (c₀ un-halved) — the colleague and interval
  paths halve it on entry (a measured 2.2e-1 root error before).
- `colleague::certified_roots` — fs-ivl interval Newton on Clenshaw
  evaluated in interval arithmetic: eligible isolated interior roots can come
  back CERTIFIED (unique-root proofs, widths ~6e-15 measured). Multiple or
  endpoint roots, `min_width` termination, and the finite subdivision budget
  can return honest `Possible` boxes. Returned boxes are mapped back to the
  physical Cheb1 domain.
- `cheb2::Cheb2` — Chebfun2-style adaptive cross approximation:
  deterministic max-residual pivots, rank-1 slice updates at FIXED
  resolution (ACA residual slices carry absolute cancellation noise,
  so the adaptive builder's machine-plateau test cannot pass on them
  — measured panic, documented in-module), exact rank on separable
  fixtures, spectral rank convergence on smooth ones, separable
  integration.
- `fourier::FourierSeries` — trigonometric interpolants on [0, 2π)
  via fs-fft's real transform (power-of-two samples): eval,
  ik-multiply differentiation (Nyquist zeroed, the real-signal
  convention), integral off c₀, tail-magnitude spectral-decay probe.

`tests/variants.rs`: cheb-101 colleague vs subdivision vs analytic;
cheb-102 even-multiplicity recovery; cheb-103 1e-3 clustered roots;
cheb-104 physical-domain interval certification (+ honest
non-certification of the double root); cheb-105 ACA
ranks/accuracy/integral; cheb-106 Fourier recovery/derivative/decay/
Bessel integral; cheb-107 bitwise replays; cheb-108 fail-fast guards
for invalid policies, domains, samples, widths, and public spectral
structs.

## No-claim boundaries
- Even-multiplicity roots: the colleague path recovers the retained admitted
  fixtures under its declared filtering policy, including a case missed by the
  sign-grid path. It does not establish generic recovery, multiplicity, root
  count, or completeness. The v1 sign-grid `roots` keeps documented no-claims
  for even-multiplicity,
  clustered, and multiple/ill-conditioned roots, and does not certify that its
  returned vector is complete. It remains the zero-dependency fallback.
- No 3D low-rank (2D ships; tensor-train is the successor), no
  complex-root REPORTING policy (real-only surfaced, documented), no
  Fourier rootfinding-on-the-circle, no Qty-dimensioned functions,
  no FrankenScipy cross-checks yet.
- `mul` may overshoot the minimal degree (resample-based); fine for
  correctness, recorded for the perf lane.
- Budget coverage after slice 2: EVERY module now has an exact-u128
  admission preflight (`admit_adaptive_build`, `admit_dirichlet_eigs`,
  `admit_root_scan`, `admit_colleague_roots`, `admit_cheb2_build`,
  `admit_fourier_build`, `admit_growth_rates` — the last four typed-
  refuse the shape violations the classic constructors panic on), and
  the heaviest paths have Cx-threaded budgeted twins (adaptive build,
  Dirichlet eigensolve, sign-grid root scan, colleague candidates).
  `colleague_roots_budgeted`, `cheb2_build_budgeted`, and
  `growth_rates_budgeted` poll at the boundaries AROUND their
  admission-bounded tiles — each tile is one non-preemptible unit whose
  classic path runs unchanged inside, numeric-evidence asserts
  (exponent-span normalization, non-finite samples, zero pivots,
  fixture-scale eigensolver convergence) retained and documented;
  `fourier_build_budgeted` is a faithful fallible MIRROR (bitwise
  spectrum parity, typed non-finite/shape refusals, 4096-sample polls).
  Algebra (`add`/`mul`) and calculus budgeting remain the bead's
  outstanding scope; the classic APIs keep their panicking contracts
  unchanged for existing callers.
- The abstract op counts in receipts are ADMITTED worst-case bounds,
  not measured cycle counts; no performance claim is attached.
- `EigsRun` values are fixed-sweep estimates. `Complete` means all
  requested shifts finished; neither complete nor cancelled output is a
  convergence, residual, ordering, or uniqueness certificate.
