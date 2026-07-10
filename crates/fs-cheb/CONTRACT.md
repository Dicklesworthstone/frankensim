# CONTRACT: fs-cheb

> Status: PARTIAL — the 1D core, collocation, ORR–SOMMERFELD,
> 2D low-rank, Fourier-periodic, and colleague-root sections are in
> force; 3D low-rank and Qty integration remain follow-up scope.

## Purpose and layer
Chebfun-style function objects (plan §6.5): smooth 1D functions as
adaptively truncated Chebyshev expansions, plus spectral collocation
differentiation matrices. Layer: **L1**. Deps: fs-fft (DCT/FFT pair),
fs-la (LU/eigen paths), fs-math (strict elementary functions), and
fs-ivl (interval root certification).

## Public types and semantics
- `Cheb1` — coefficients over FIRST-KIND Chebyshev points (roots grid):
  values ↔ coefficients is exactly fs-fft's DCT-II/III pair. `build`
  doubles the grid until the trailing quarter of coefficients sits at
  the machine-precision plateau (8.9e-16 relative), then truncates;
  unresolvable functions panic at `max_degree` with a structured
  message. `eval` (Clenshaw), `differentiate` (coefficient recurrence,
  domain chain rule), `integral` (even-coefficient formula), `add`,
  `mul` (resample + rebuild), `roots` (subdivision + bisection +
  Newton polish).
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
3. Plateau detection does NOT chase noise floors (tested with a
   deterministic ~1e-18 jitter fixture).
4. `diff_matrix` rows sum to exact zero (differentiation annihilates
   constants bitwise).
5. Deterministic per ISA: all state is built on strict fs-math cos/sin
   and fixed-order arithmetic. The radix-4 downstream golden is recorded
   on M4 Pro; x86-64 equality is armed pending RCH admission (the previous
   radix-2 golden was trj-verified).
6. Colleague roots agree with the subdivision scanner on simple roots,
   recover even-multiplicity roots the scanner cannot see, and
   certification boxes are reported in physical-domain coordinates.
7. `Cheb2` captures separable rank exactly on fixture functions, keeps
   deterministic pivot tie-breaking, and converges spectrally on the
   smooth non-separable fixture.
8. `FourierSeries` exactly recovers trigonometric fixture modes,
   differentiates `sin` to `cos`, and uses c₀ for the periodic integral.

## Error model
Structured panics for programmer/modeling errors: non-finite or
inverted domains, non-finite samples/coefficients, unresolved functions
at `max_degree`, domain mismatches in algebra, invalid colleague
policies, non-positive certification widths, malformed public `Cheb2`
or `FourierSeries` fields, and non-power-of-two Fourier sample counts.

## Determinism class
Bit-deterministic per ISA by construction. The golden hashes coefficients
+ integral + derivative sample + roots + collocation eigenvalues. The
current radix-4 value is recorded on aarch64-apple and must match on
x86-64 before the cross-ISA row is restored to verified.

## Cancellation behavior
Construction is bounded (max_degree cap); no poll points needed at v1
scales.

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
tests/cheb_battery.rs (recovery, calculus, plateau robustness, roots,
collocation accuracy, eigen demo, golden hash).

## Variants (bead kw89)

- `colleague::colleague_roots` — the Chebyshev companion matrix
  (three-term-recurrence rows, coefficient-loaded last row scaled by
  −1/(2aₙ)), eigenvalues via the fs-la complex nonsymmetric stack,
  filtered by a DOCUMENTED [`ColleaguePolicy`] (trailing-coefficient
  trim, imaginary tolerance, domain slack, √ε-scale cluster dedupe —
  a double root's eigenvalue pair splits at ~5e-9, measured). This
  RESOLVES the v1 even-multiplicity no-claim: (x−r)²(x−s) fixtures
  the sign scanner provably misses are found (cheb-102). Cheb1 stores
  the Σ′ convention (c₀ un-halved) — the colleague and interval
  paths halve it on entry (a measured 2.2e-1 root error before).
- `colleague::certified_roots` — fs-ivl interval Newton on Clenshaw
  evaluated in interval arithmetic: simple roots come back CERTIFIED
  (unique-root proofs, widths ~6e-15 measured); multiple roots come
  back honestly `Possible` (their derivative encloses zero, as it
  must). Returned boxes are mapped back to the physical Cheb1 domain.
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
- Even-multiplicity roots: RESOLVED for the colleague path (above);
  the v1 subdivision `roots` keeps its documented limitation and
  remains the zero-dependency fallback.
- No 3D low-rank (2D ships; tensor-train is the successor), no
  complex-root REPORTING policy (real-only surfaced, documented), no
  Fourier rootfinding-on-the-circle, no Qty-dimensioned functions,
  no FrankenScipy cross-checks yet.
- `mul` may overshoot the minimal degree (resample-based); fine for
  correctness, recorded for the perf lane.
