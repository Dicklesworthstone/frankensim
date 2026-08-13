# fs-vfit — Contract

## Purpose and layer

Layer L3. Passive rational approximation of tabulated frequency
responses (bead frankensim-fsim-vector-fitting-pvv40): the bridge that
turns offline frequency-domain results — BEM radiation loads, TMM bore
impedances, coupled-body mobilities — into runtime-realizable stable
filters. Identification (two independent front ends), passivity
certification with convex repair, and prewarped bilinear discretization
to parallel biquad sections and discrete state-space.

## Conventions (load-bearing)

- Laplace axis `s = i*omega` (`e^{+i*omega*t}`); real impulse responses
  by conjugate-closed storage (`PoleTerm::Pair` stores only the
  `Im > 0` member — conjugacy cannot drift).
- FrankenSim's acoustics stack (fs-duct, fs-bem) uses `e^{-i*omega*t}`,
  the COMPLEX CONJUGATE convention: conjugate such data before fitting.
  Executed failure pinned in the clarinet casebook: unconjugated
  acoustic data rotates phase the wrong way through every resonance and
  NO stable rational model fits it (~50% error floor for every method).
  `|H|` and `Re H` are conjugation-invariant.
- Passivity is impedance-form positive realness `Re H(i*omega) >= 0`
  with `d >= 0`, `e >= 0`.

## Public types and semantics

- `model::RationalModel` — pole-residue terms + direct `d` + improper
  `s*e`; `eval`/`eval_iw`; real block-diagonal `state_space()` whose
  independent LU-solve `StateSpace::eval` route is the
  realization-parity oracle.
- `vf::vector_fit(omega, h, opts)` — relaxed vector fitting
  (Gustavsen-Semlyen): homogeneous sigma system with the relaxed
  non-triviality row, pole relocation via the sigma-zero eigensolve on
  the real block realization, stable-pole flipping, conjugate closure.
  `d >= 0` / `e >= 0` enforced IN the solve — exact by enumeration of
  all four bound combinations, keeping the feasible minimum. Data normalized to O(1) by geometric-mean magnitude and
  every LS column-equilibrated (physical scales like |Z| ~ 1e7
  otherwise poison the QR). Weight presets: uniform, inverse-magnitude
  (honest at antiresonances), log-band — recorded in `FitReport`.
- `vf::fit_auto_order` — ascending-order selection; a plateau
  terminates the ascent only AFTER a >= 10x cliff has been seen
  (executed observation: under-resolved orders sit on a
  pre-convergence plateau before the cliff at the true order); the
  noise floor terminates unconditionally (overfit refusal). Returns
  the full order-vs-error curve.
- `loewner::loewner_fit` — the INDEPENDENT second front end: Loewner /
  shifted-Loewner pencil on a subsampled grid, conjugate-augmented and
  real-transformed, SVD rank reveal, projected-pencil eigenvalues;
  direct terms handled by ITERATED STRIPPING (crude d/e estimate,
  strip, pencil, refit d/e via the shared residue pass, repeat —
  executed: an unstripped d biases poles ~6e-4 relative, Q-amplified
  to percent-level response error). Residues always from the same
  final LS pass as vector fitting; nothing upstream of the pole
  estimates is shared.
- `loewner::cross_check` — pole and response agreement diagnostics
  between the two front ends.
- `passivity::check_passivity` — grid arm (log sweep + per-resonance
  refinement, violations compressed to per-band representatives) plus
  the Hamiltonian eigenvalue crossing test (exact) when `d > tol`;
  otherwise the certificate is the named weaker `GridOnly` class.
  Descriptor-form statement: on the imaginary axis `Re(i*omega*e) = 0`
  — the improper term is lossless and the Hamiltonian test applies to
  the proper part unchanged; general improper descriptor passivity
  (even matrix pencils) is a no-claim.
- `passivity::repair_passivity` — convex residue perturbation
  (poles, d, e FIXED): active-set QP (KKT equality solves via fs-la
  LU) over accumulated violation frequencies, iterated with
  re-certification; reports rounds, relative perturbation, and the
  final KKT stationarity residual.
- `discretize::bilinear` — prewarped bilinear transform to PARALLEL
  first-order/biquad sections (`K = omega_pw / tan(omega_pw T/2)`);
  exact at the prewarp frequency; the improper term maps to the exact
  lossless differentiator section (pole at `z = -1`, admitted by the
  stability check as marginal BY DESIGN). `eval_f32_quantized` is the
  coefficient-quantization probe. `bilinear_state_space` — the same
  map in Tustin state-space form; the improper coefficient is returned
  as `e_leftover` for the caller's extra section (a named seam, not a
  silent drop).
- `discretize::DigitalFilter::step` — transposed DF-II runtime of the
  parallel bank. `reflectance` is the scattering map
  `(Z−Zc)/(Z+Zc)`. `modulate_delay` peels or applies a known
  bulk delay. `realize_tabulated` is the foundry-to-runtime door:
  vector-fit + bilinear. `realize_tabulated_impedance` is the same
  door after `repair_passivity` (`Re Z ≥ 0`); repair exhaustion
  keeps the raw stable fit. `DigitalFilter::enforce_abs_bound` is
  the scattering projection `|H| ≤ bound`. `DelayedFilter` is a
  characteristic port (delay ⊕ residual filter) with
  `enforce_scattering_passivity` (`|R| ≤ 1` on a caller grid). A
  TMM bore, a muffler, a pulse tube, and a BEM radiation load are
  fillings of that port. Music is not a special case.

## Invariants

1. Realization parity: partial-fraction evaluation and the state-space
   LU route agree to 1e-12 relative (independent arithmetic).
2. Known-answer identification: a 6-pole system is recovered from
   clean samples to 1e-9 relative response and 1e-6 relative poles.
3. Fitted models are stable and conjugate-closed by construction.
4. The two front ends agree on clean data and disagree DIAGNOSTICALLY
   on aliased data (the cross-check is the artifact detector).
5. A repaired model's poles are bitwise the input model's poles; only
   residues move.
6. Determinism: identical inputs refit bitwise (no RNG anywhere;
   deterministic starting poles, fixed iteration caps).
7. Scaling covariance: fitting `a*H` scales the linear parameters by
   `a` and moves no pole.

## Error model

Typed enums per module (`VfError`, `LoewnerError`, `PassivityError`,
`DiscretizeError`); degenerate input, order, eigensolve failure, sigma
collapse, repair exhaustion, and beyond-Nyquist are all named
refusals. No silent degradation.

## Determinism class

Deterministic: fixed grids, fixed iteration caps, canonical eigenvalue
ordering from fs-la, no time or RNG dependence. Bitwise repeatability
asserted in the battery.

## Cancellation behavior

Pure synchronous computation; longest call (24-pole auto-order fit on
2400 samples) is seconds-class. No `Cx` integration (workspace
`frankensim-ccmn` effort).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/vfit.rs` (12): realization parity (two routes, on and off
axis); known-answer 6-pole recovery; noisy fit at the noise floor;
fit-of-a-fit idempotence; scaling covariance; bitwise determinism;
Loewner cross-check clean vs aliased; passivity detect + repair on an
executed near-miss ACTIVE fixture (Hamiltonian crossings present
before, empty after; poles unmoved; KKT residual near zero);
discretization (exact prewarp point; authored 8% in-band envelope
derived from Q-amplified bilinear warp; dropped-prewarp mutation 10x
visible; state-space vs sections parity; non-vacuous f32 probe);
auto-order plateau + overfit refusal; weight-preset visibility (the
inverse-magnitude fit must win at the antiresonance); typed refusals.

Inline discretize runtime tests: DF-II first-order recurrence,
reflectance scattering map, delay peel of a pure `e^{-iωτ}`, identity
`DelayedFilter` delay, tabulated-constant realization.

`tests/characteristic_line.rs` (1): an open viscothermal cylinder's
`Z_in` becomes a characteristic `DelayedFilter`; an impulse returns
inverted near `2L/c` with `|R| < 1`. The same port as a muffler.

`tests/clarinet_casebook.rs` (1): fs-duct TMM clarinet-class bore
(USSA-1976 air, Zwikker-Kosten losses, unflanged radiation) -> 24-pole
fit (order curve logged; peaks within 2 cents authored, measured
0.17) -> passivity certification -> 192 kHz prewarped bilinear biquad
bank re-measured from the digital filter (all 10 peaks within 2 cents
authored, measured 1.22 worst) -> f32 sensitivity logged. JSON-lines
at every stage.

## No-claim boundaries

- SISO only; matrix-valued (MIMO) vector fitting and pencil
  identification are follow-ups.
- Scattering-form fitting (`|S| <= 1` enforcement for waveguide
  reflection filters) is a recorded follow-up; the impedance form is
  the v1 passivity law.
- The Loewner front end assumes proper-dominant data; the iterated
  stripping handles polynomial parts d and s*e, but heavier improper
  behavior is out of scope.
- Passivity certification is exact only through the Hamiltonian arm
  (`d > 1e-12`); the `GridOnly` class is honest about its weaker
  guarantee (a dip narrower than the grid can hide). Consumers reading
  only the `passive` bool inherit that weaker guarantee — check
  `class` when the certificate strength matters.
- The QP repair minimizes residue movement, not a band-weighted
  perceptual metric; a perceptually weighted repair is a follow-up.
- Discrete-time (z-domain) DIRECT vector fitting is a follow-up; v1
  discretizes the continuous fit.
