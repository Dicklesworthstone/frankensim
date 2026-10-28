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
  silent drop). Stable continuous poles above the sampled band remain
  admissible because Tustin maps the complete open left half-plane inside
  the unit circle; only sampled/prewarp frequencies are Nyquist-bounded.
- `discretize::DigitalFilter::step` — transposed DF-II runtime of the
  parallel bank. `reflectance` is the scattering map
  `(Z−Zc)/(Z+Zc)`. `modulate_delay` peels or applies a known
  bulk delay. `realize_tabulated` is the foundry-to-runtime door:
  vector-fit + bilinear. `realize_tabulated_impedance` is the same
  door after `repair_passivity` (`Re Z ≥ 0`); repair exhaustion
  keeps the raw stable fit. `DigitalFilter::enforce_abs_bound` is
  the scattering projection `|H| ≤ bound`. `DelayedFilter` is a
  characteristic port (delay ⊕ residual filter). `scattering_peak`
  measures `max |R|` on a caller grid and changes nothing; render callers
  refuse a peak above 1. DEF-03: the global `1/peak` gain projection
  `enforce_scattering_passivity` is a `#[cfg(test)]` diagnostic only.
  `from_impulse_response` is the exact sampled-`R(ω)` or sampled-`Z(ω)`
  FIR filling of the same port (scattering vs impedance). A TMM bore,
  a muffler, a pulse tube, and a BEM radiation load are fillings of
  that port. Music is not a special case.

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

- Dispersion stage (music bead 3ez8g.7.2): `DelayedFilter::with_dispersion`
  adds a cascade of identical first-order allpasses
  `H(z) = (a + z^-1)/(1 + a z^-1)` to the line output — stiffness
  realized as frequency-dependent loop delay instead of modal
  frequencies. `|H| = 1` EXACTLY at every frequency, so the stage is
  passive BY STRUCTURE (it reshapes phase, never energy); `|a| >= 1`
  and zero sections refuse; the stage is inert unless requested, so
  every existing consumer is bit-unchanged.
  `first_order_allpass_group_delay` publishes the closed form the
  design math solves against.

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

`DelayedFilter::history` / `from_impulse_response_with_history` (bead
3ez8g.4.2): the FIR-mode outgoing-wave ring buffer is readable (oldest
first) and a new exact-FIR line can be built PRE-LOADED with a prior
line's history — the D17 state-lift seam for characteristic-line
reflectance swaps (in-flight waves persist; the new taps govern future
reflections). History longer than the new buffer is truncated at the
OLD end; shorter history zero-pads the old end. IIR-mode lines return
an empty history (no lift claim).

## One-way observation FIRs (IRP M2.5-W6)

`observation::{FirKernel,ObservationFir}` carries the complete source-to-receiver
component FIR data from a declared frequency-domain observation bake. Kernels
have explicit input/output indices and exact integer delay prefixes; fractional
propagation and reconstruction latency are in the producer's taps/metadata.
No vector fitting, positive-real test, clipping or passivity enforcement is
applied. This MIMO convolution is an observation outside the reaction loop;
the SISO-only limitation below concerns rational identification, not FIR replay.

The cold caps are 4096 inputs, 256 output components, 4194304 coefficients and
4194304 shared input-history entries, independently admitted before allocation.
Source/receiver units, provenance, pass band, latency and producer uncertainty
travel in the upstream artifact; replay supplies no defaults for them. Histories
start at zero. Traversal and sums are fixed; identical same-ISA replay is bitwise
deterministic. Invalid input or output overflow leaves all history unchanged.
`observation::tests` covers delayed two-input summation, error/retry identity,
empty/index errors and each allocation-free cap/cap+1 admission.

## Multiport impedances (IRP M2.5-W3b)

`mimo::mimo_vector_fit(omega, data, ports, opts, reciprocal)` fits a row-major
`ports x ports` impedance with one shared pole set.
- It uses relaxed vector fitting with a single sigma shared across all fitted
  entries.
- With `reciprocal`, it fits the upper triangle of the data's symmetric part.
- Each frequency is weighted by the inverse Frobenius norm of its data matrix,
  so the fit controls matrix-relative error.
- The report gives the max and RMS matrix-relative misfit against the
  original data (asymmetry included) and the data's reciprocity defect.

`MimoModel::state_space` is the scalar block realization times the port
identity. `MimoStateSpace::eval` is the independent parity route.

`check_mimo_passivity` tests `Z + Z^H >= 0`.
- It uses the dense log grid plus the impedance Hamiltonian
  `[[A - B R0^-1 C, -B R0^-1 B^T], [C^T R0^-1 C, -(A - B R0^-1 C)^T]]` with
  `R0 = D + D^T`.
- A singular `R0` downgrades the certificate to `GridOnly`.
- It refuses unstable models, and models whose `E` is not symmetric PSD or
  whose `R0` is indefinite.

`repair_mimo_passivity` moves residues only.
- Each `(frequency, eigenvector)` pair is an exact linear constraint on the
  Rayleigh quotient, anchored on the original model.
- The constraints accumulate over at most 12 rounds, with re-certification
  after each round.
- It reports the coordinate perturbation and the in-band matrix-relative
  change.

`MimoModel::zoh(T)` steps the model by exact zero-order hold.
- Each pole term advances as `exp(p T) x + (exp(p T) - 1)/p u`, per input
  port. This is the state matrix's exponential, with no bilinear map.
- The improper `E` term uses the SISO first-order differentiator.
- With a held input it reproduces the continuous step response at the
  samples (`mimo_zoh_step_response_is_exact_at_the_samples_g1`).

The shared `mimo::load::{MimoLoadLimits,PreparedMimoLoad}` admits proper causal
loads for both waveguide and field propagators. It retains the source model,
declared band/limits, original MimoZoh, continuous report and separate actual
sampled interval-mean report. Both passivity checks must pass; the sampled
report remains GridOnly. It admits 1..=128 ports, at most 256 pole terms and
at most 256 real states, exact finite matrix shapes, stable finite poles,
real residues for real poles, positively oriented conjugate pairs, zero
improper E, a finite positive clock and a source band reaching Nyquist.
The original network's conservative retained-payload budget is checked before
realization. Caller buffers, allocator overhead and cold eigensolver scratch
are outside that budget.

The mean-history/direct methods expose the same accepted interval-mean map
for a simultaneous port solve. Checked step_mean_into advances original ZOH
coordinates without allocation, refusing wrong shapes, nonfinite input,
output or candidate state before accepted history or caller output changes.
Cloning the prepared owner supports a larger transaction; the caller still
owns coupled work validation and whole-step acceptance. The existing
waveguide terminal uses this admission without changing its port solve or
rational integrator. Tests cover original-runtime parity, resource/shape/
passivity refusals, finite-effort state overflow, refusal/retry and abandoned
candidate history preservation.

## Continuous access to an admitted MIMO load

`PreparedMimoLoad::continuous_operator()` borrows the original admitted source.
Its dimension and port count preserve the complete real pole realization;
`state_space(maximum_dense_bytes, cancelled)` calls the existing
`MimoModel::state_space` without balancing, reduction, refitting or timestep
conversion. The explicit matrix payload cap is
`8*(n*n + 2*n*ports + 2*ports*ports)` bytes, checked before allocation under
the original 256-state/128-port bounds. Cancellation brackets this bounded
cold construction, and nonfinite realized coefficients refuse. The instantaneous
`D` is distinct from the original runtime's interval-mean direct matrix. Static
and singular-`D` sources remain representable; each downstream interconnection
owns its algebraic-solvability requirements.

A real pole contributes one coordinate per input port. A conjugate pair
contributes all `Re(z)` followed by all `-Im(z)` in original input-port order,
matching the original real block `[[Re(p),Im(p)],[-Im(p),Re(p)]]`.
`real_coordinates` reads this chart; `candidate_from_coordinates` returns a
candidate of the same original owner with finite, complete endpoint coordinates
and endpoint input. Publication belongs to the enclosing work-checked
transaction. The source, physical clock, limits and both reports are retained.
The next original ZOH step continues the adopted history without resetting poles.

Continuous frequency queries must remain inside the source's declared band.
The view does not change the original source, accepted pole history or reports.
In particular, GridOnly admission is not promoted to a physical storage or
dissipation certificate; terminal work remains an explicit boundary exchange
unless a separate physical storage owner is supplied.

`delayed_mimo_vector_fit` fits the entries after removing caller-declared
delays (`Z_ab = exp(-s tau_ab) Zhat_ab`).
- A delayed entry never carries the improper `s E` term: under a delay that
  term puts `2 w E sin(w tau)` into the Hermitian part.
- `check_delayed_mimo_passivity` and `repair_delayed_mimo_passivity` work on
  the restored model, grid only.

The tests check the following:
- known-model recovery (`< 1e-8`) and pole recovery;
- realization parity;
- an exact Hamiltonian certificate on a passive model;
- a planted violation located by the Hamiltonian and repaired;
- delayed known-answer recovery;
- refusals.

fs-bem `mimo_rational_order_ladder` realizes a real aperture matrix.

## Causality of sampled data (IRP M2.5-W3b)

`causality::check_causality(omega, data, opts)` tests one response sampled on
a uniform grid that starts at 0 (with a real DC sample) or at half a step.
- It reconstructs every sample below `0.9 W` from the generalized dispersion
  relation with `opts.subtractions` subtraction points (Triverio &
  Grivet-Talocia). The band is extended by Hermitian symmetry, and the
  integrand's poles are integrated analytically.
- A miss larger than the sum of three declared terms refuses with
  `NonCausal`, naming the worst frequency. The terms are: the out-of-band
  tail under `|H| <= tail_bound_ratio max |H|`; twice the trapezoid
  second-difference estimate; and `relative_data_tolerance` propagated by the
  triangle inequality.
- `check_mimo_causality` tests each entry after its declared delay is
  removed. A delay beyond the first arrival leaves an advanced remainder and
  refuses with that entry.

Tests (`causality::tests`): causal rational and delayed responses pass on both
grids; 0.5 ms and 0.2 ms advances and conjugated (`exp(-i w t)`) data refuse; a
pulsating sphere's radiation impedance passes while its conjugate refuses; an
over-extracted matrix delay names its entry; and out-of-range grids, data and
options refuse. A pass is consistency with causality under the declared
assumptions, not a proof.

Resolution: an anti-causal part enters the error as twice its miss from the
subtraction interpolant. Non-causality is resolved only on the band's own
time scale, about `1/f_max`. On the 5 kHz test band a 50 us advance passes
(yto c1, worst ratio 0.17). A delay over-extracted by less than that, such as
the centre distance instead of the gap for finite pistons, is not caught.

## No-claim boundaries

- MIMO covers vector fitting, realization, passivity and the causality test
  of its data. MIMO Loewner pencil identification is a follow-up. The
  closed-loop audit of an in-loop MIMO load is fs-couple
  `closed_loop_audit_reed_bore_radiation`.
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
