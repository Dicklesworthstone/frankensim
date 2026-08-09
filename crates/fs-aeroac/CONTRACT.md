# fs-aeroac — Contract

## Purpose and layer

Layer L3. Low-Mach aeroacoustic source models (bead
frankensim-fsim-aeroacoustic-sources-9ok02, IN PROGRESS — this crate
is the first slice): the acoustic-analogy half of the honest hybrid
over fs-lbm's incompressible base flows. Everything here is a
SHAPE/SCALING authority; ABSOLUTE SPL IS NEVER CLAIMED (the crate's
scope law, pinned as data and by test).

## Public types and semantics

- `bessel::{j0, j1, y0, y1, hankel0_outgoing, hankel1_outgoing}` —
  deterministic Bessel/Hankel functions built EXCLUSIVELY from
  fs_math primitives (det transcendentals, Dd double-double series
  accumulation below the x = 18 crossover, recurrence-GENERATED
  asymptotic sums above it — deliberately no transcribed coefficient
  tables). `e^{-i omega t}` convention: the outgoing cylindrical wave
  is `H^(1) = J + iY`. Y at non-positive arguments is NaN (real
  domain).
- `curle2d::dipole_pressure` — 2D frequency-domain Curle compact
  dipole over the outgoing Hankel Green's function
  `G = (i/4) H0^(1)(kr)`: `p = (ik/4) H1^(1)(kr) (rhat . F)` per unit
  span. The 2D-vs-3D Green's-function trap from the bead's polish
  round is structurally excluded (frequency domain, Hankel kernel);
  every output embeds [`SCOPE_STATEMENT`].
- `bickley::bickley_rayleigh_mode` — Rayleigh-equation shooting
  solver (RK4 + complex secant, half-line with symmetry boundary
  conditions and the exact `phi' = -alpha phi` far-field decay) for
  the Bickley jet `U = sech^2(y)`: the inviscid instability oracle
  the fs-lbm jet validation chain runs against. Non-convergence
  REFUSES with the residual disclosed; no partial eigenvalue.
- `bickley::rayleigh_residual_closed_form` — the self-verification
  surface for the analytic pins.
- `SCOPE_STATEMENT` — the no-absolute-SPL / 2D-to-3D-span-correction
  law as data (the marketing-mutation guard asserts it).

## Invariants

1. EXACT Wronskian: `J1 Y0 - J0 Y1 = 2/(pi x)` holds within 5e-13
   over a log grid 1e-3..3e3 including the series/asymptotic
   crossover from both sides. Honest limits (review-measured against
   a 120-digit reference): the identity is blind to COMMON phase
   errors, and pointwise accuracy beyond the certified grid degrades
   as ~x*eps from f64 phase representation (~2e-11 at x = 1e6) —
   stated in the module doc; intended kr ranges sit inside the
   certified band.
2. Independent-path derivatives: central-difference `J0'` matches
   `-J1` and `Y0'` matches `-Y1` at FD accuracy (1e-8) across both
   regimes.
3. Cross-implementation oracle: fsci-special's cephes-heritage
   j0/j1/y0/y1 agree within 5e-12 of envelope (fsci's plain-f64
   Y-series carry ~1e-12 cancellation near their own crossover —
   measured; the Wronskian is the precision arbiter).
4. Small-argument limits exact (J0(0) = 1, J1 ~ x/2, bounded Y0 log
   constant); far-field `|H0| -> sqrt(2/(pi x))` within the O(x^-2)
   envelope correction.
5. Curle dipole physics: far-field decay exponent is -1/2
   (CYLINDRICAL spreading; the 3D-Green's-swap mutation reads -1 and
   fails the gate), the perpendicular directivity null is EXACT, the
   45-degree amplitude matches cos(theta) within 1e-3, and the scope
   statement rides on every output.
6. Bickley analytic pins are SELF-VERIFIED per run: `phi = sech^2` at
   `(alpha, c) = (2, 2/3)` and `phi = sech tanh` at `(1, 2/3)` drive
   the closed-form Rayleigh residual to machine zero (< 1e-14) at
   every probe point, and the residual is proven LIVE (perturbing
   alpha or c by 0.01-0.05 moves it above 1e-3). No literature value
   is transcribed on trust.
7. The shooting solver reproduces the pins from the unstable side
   (alpha = 1.95 sinuous: c within 0.02 of 2/3 with small positive
   Im c; alpha = 0.95 varicose: within 0.03), shows positive
   mid-band growth exceeding the near-neutral growth, and the
   eigenvalue is grid-converged (2048 vs 4096 RK4 steps within
   1e-8).
8. VERIFICATION CHAIN link 1 (dev-dep composition with fs-lbm): the
   D2Q9 Bickley jet's measured linear growth rate CONVERGES toward
   the inviscid Rayleigh oracle under simultaneous grid + Reynolds
   refinement — measured -21.4% (b = 8 lu, Re 240) closing to -12.0%
   (b = 12 lu, Re 432), both biased LOW as finite-Re physics demands,
   with clean-exponential window consistency, low-Mach and linearity
   diagnostics green per run, and gates authored just outside the
   measured values (a fixed loose envelope would hide a broken
   convergence direction). Executed lesson in the test doc: the free
   jet's own viscous diffusion (~nu/b^2 per step) decays the base
   flow and drags the measured rate down over long runs — fit early,
   at high Re.
9. Refusals typed by name (non-positive wavenumber, coincident
   observer/source, non-finite inputs, bad solver parameters,
   non-convergence with residual); Hankel functions outside their
   domain are FULLY NaN (both components — a half-valid complex once
   leaked, review-caught); the Rayleigh solver's result is
   documented GUESS-DEPENDENT (any eigenvalue of the symmetry class
   zeroes the mismatch); determinism bitwise.

## Error model

Typed `AeroacError` (non-finite, invalid-parameter, not-converged
with disclosed residual). No silent degradation, no fabricated
eigenvalues.

## Determinism class

Deterministic: fs_math::det transcendentals + Dd arithmetic only, no
platform libm, no RNG.

## Cancellation behavior

Synchronous, milliseconds-class. No `Cx` integration (workspace
`frankensim-ccmn`).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/bickley_lbm.rs` (1): the LBM-vs-Rayleigh growth-rate
convergence fixture (invariant 8).

`tests/aeroac.rs` (11): Wronskian identity; derivative cross-checks;
fsci-special oracle; small-argument limits; Hankel far-field
amplitude; Curle spreading + directivity + scope guard; Curle typed
refusals; Bickley pin self-verification + falsifiers; shooting solver
vs pins + grid convergence; Bickley typed refusals; bitwise
determinism.

## No-claim boundaries (the bead's remaining scope — OPEN)

- fs-lbm source EXTRACTION (surface-pressure spectra from jet-labium
  runs) and the absorbing-layer treatment with a measured reflection
  coefficient: not implemented (the bead's pilot showed 6% flux
  imbalance at Re 200 from outlet reflections — spectra are NOT
  trustworthy before the sponge lands).
- Edge-tone Strouhal staging vs published data: not implemented.
- Fitted flute-noise tables + demo-synth consumption: not
  implemented.
- Quadrupole (volume) sources: dipoles dominate at low Mach over
  rigid surfaces; the Lighthill volume term is out of scope v1.
- Only orders 0 and 1 of the Bessel functions are provided (all the
  2D monopole/dipole machinery needs); higher orders are not
  implemented.
- The Rayleigh oracle is INVISCID (no critical-layer viscous
  correction); exactly-neutral real-c modes are ill-posed for
  shooting and are validated by approach from the unstable side.
