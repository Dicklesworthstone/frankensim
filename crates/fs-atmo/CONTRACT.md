# fs-atmo — CONTRACT

Bead: frankensim-wf-root-guzez.4.5 (E3.3a, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md §5.4
(ROUND 6 steady state). Evidence: air-state-v1.json (E1.8) under the frozen
registry (E1.7). Conventions: frame-conventions-v1; heights are altitude h
above the aerodynamic ground plane.

## Purpose and layer

L2. The atmosphere foundation: FlatSiteLogLaw mean wind, wall-compatible
solenoidal modal turbulence with exact analytic derivatives, and the
provenance-bound `sample_air_state` API. E3.3b adds the exact-discrete OU
amplitude evolution (sequential checkpointed state) and a fit of selected
one-point stresses from the Mann-class spectral tensor; E3.3c the optional
FetchAdjustedMassConsistent mode.

## Public types and semantics

- `FlatSiteLogLaw { scenario_effective_z0, displacement_height,
  reference_height, reference_speed }` — ONE scenario-level z₀; `speed(h)`
  = (u*/κ)·ln((h−d)/z₀) above the sublayer, 0 inside; analytic `dspeed_dh`.
- `TurbulenceField::build(seed, n_modes, sigma, length, u_adv)` — philox
  counter-addressed modes (StreamKey{seed, ATMO_KERNEL, tile = mode});
  `sample(x, y, h, tick)` → velocity + EXACT analytic gradient. Per mode:
  ψ_h ∝ sin(k_h·h)·cos θ, ψ_v ∝ cos(k_h·h)·cos θ, u = ∇×ψ, with the
  frozen-phase clock θ = k_x(x − U_adv·t) + k_y·y + φ.
- `AirScenario` / `DEC17_AIR` — ρ/μ/T/p with provenance (E1.8 derivations).
- `Atmosphere::sample_air_state(x, y, h, tick)` → `AirState` whose
  `dynamic_pressure_pa()` and `reynolds(chord)` derive from the SAME
  provenance-bound state as the velocity.
- `Atmosphere::sample_flow` and `EvolvingAtmosphere::sample_flow` expose
  the exact original mean-plus-turbulence velocity and gradient without
  imposing constant scenario thermodynamics. L3/L4 material compositions
  attach their own coherent gas/profile; no upward dependency is introduced.
  Original air-state sampling delegates to this same flow evaluation.
- `evolving::EvolvingTurbulenceField::at_anchor(basis, seed, times)` —
  opt-in space–time field with one declared OU correlation time per spatial
  mode. Two independent OU quadratures per mode multiply the vector potential
  before its curl; both have stationary variance 1/2. This preserves the
  basis realization's complete phase-averaged stress tensor and makes its
  conditional one-point statistics horizontally homogeneous. Per-mode
  covariance decays exponentially while its phase advects with the declared
  mean speed. `sample` reads the current signed tick; `sample_batch` fills
  caller-owned storage without allocations or random draws. `advance_to`
  alone evolves the sequential state. A self-contained checkpoint retains
  basis, times, and both OU paths; resume is bitwise.
- `Atmosphere::into_evolving` / `evolving::EvolvingAtmosphere::sample_air_state`
  — the same mean law and provenance-bound air properties with the evolved
  vector velocity and analytic gradient, so dynamic pressure and Reynolds
  use that same instantaneous field.
- `fit::{fit_amplitudes, build_fitted_field}` fits (u², v², w², uw) at a
  declared reference height, with a positive-semidefinite amplitude covariance
  and explicit finite-basis residuals. The fit retains its seed, Mann target
  and target shear stress. Building a realization checks that binding and
  reuses the retained stress instead of repeating the target quadrature.
  Γ = 0 fits the variances with zero cross-correlation, without dividing by a
  zero shear-stress target. Nonfinite heights/covariances/advection refuse.
- `mann_field::{MannSpectralModel, HomogeneousMannField}` provides a separate
  homogeneous prescribed velocity field from the full Mann spectral tensor.
  Each finite-band wavevector has two transverse polarizations whose complete
  covariance equals its weighted tensor, including cross-components. Four
  independent unit-variance OU paths supply both polarizations' sine/cosine
  quadratures. Geometry and amplitude seeds are separate; vector advection,
  decorrelation times, mode count and both kL cutoffs are explicit. Initialization
  performs the tensor algebra; velocity/gradient samples allocate no storage.
  `covariance` exposes the exact represented spatial or space-time tensor,
  including finite-band and quadrature effects. Checkpoints retain model inputs
  and all four paths and rebuild the deterministic model on resume. The
  original homogeneous observation has no wall condition; z = 0 is an ordinary
  plane. The explicit plane-blocked observation below adds its own boundary
  treatment without changing the retained model or quadratures.
- `mann_field::fixed_step::FixedStepMannField` starts from a retained homogeneous
  snapshot and evolves its four quadratures at an explicitly supplied physical
  time step. Each mode uses the exact OU transition
  `rho = exp(-dt/tau)`, `sigma = sqrt(-expm1(-2dt/tau))`, following Gillespie
  (1996), DOI `10.1103/PhysRevE.54.2084`. Its distinct counter-addressed stream
  is bound to the authored step size and seed; changing the time step does not
  preserve the same realization or the original 120 Hz path. The original
  integer-tick owner and its checkpoint semantics remain intact.
  Sampling reuses the same spectral tensor, transverse polarizations and
  spatial filters, with advection applied once at the retained phase origin
  plus elapsed physical time. The immutable model is shared; accepted/trial
  amplitudes and transition coefficients are prepared once. Transition work
  and actual numerical capacities are admitted explicitly, and cancelled or
  refused advances retain the accepted state. Checkpoint/resume and temporal
  state copying bind the complete model, seed, clock and correlation times.
  The nominal numerical payload is 88 bytes per mode, excluding immutable
  spectral geometry and fixed metadata; shared coefficient storage is counted
  conservatively per owner. Authored OU decorrelation is a synthetic temporal
  model, not resolved nonlinear eddy turnover, and is not temporally bandlimited.

### Full-tensor Mann observations above a plane

`mann_field::wall::MannPlaneWall` declares a finite plane `z = height_m` in
the authored Mann frame. Both temporal owners expose `sample_wall_blocked`
and `covariance_wall_blocked`; the spectral model also exposes equal-time
covariance. Sampling preserves the original accepted quadratures and adds the
harmonic, irrotational correction of Hunt & Graham (1978),
<https://doi.org/10.1017/S0022112078000130>. With horizontal wavenumber magnitude
`q` and height `h`, a modal amplitude `a` receives
`a_z * (i*kx/q, i*ky/q, -1) * exp(i*k.wall_point - i*k.U*t - q*h)`.

The correction enforces zero normal velocity at the plane, preserves
solenoidality and decays aloft. All nine derivatives are analytic. Covariance
transforms both original polarizations and all cross-components at both
observation points, retaining the original OU decay and phase advection.
Both coordinates are required because normal homogeneity is lost. The common
first-point phase is removed before covariance evaluation. Near-wall
differences use half-angle identities and `expm1`; exact vertical modes have
a separate zero-correction branch, and an exactly vanished correction does
not evaluate an irrelevant distant-wall phase.

Only finite points in the exact upper half-space are admitted. The same
retained modal arrays are read without new allocations, draws or clock
changes. This is an outer-flow slip-plane construction, without a viscous
no-slip layer, ground-generated turbulence or a calibrated near-ground
spectrum. Point observations are supported; homogeneous sinc filters do not
supply the evanescent correction's spatial average. Six focused Rust
regressions accompany this source-reviewed implementation; native compilation
and execution were unavailable during the change.

## Invariants

- Wall-compatible solenoidality: u is a curl — ∇·u ≡ 0 analytically; the
  analytic-gradient trace cancels to machine precision (V-04a relative bound).
- Wall-compatible parity: sin(k_h·0) = 0 makes the vertical component vanish
  at h = 0 IDENTICALLY (bitwise zero, battery-asserted).
- Pointwise z₀ insertion is FORBIDDEN while solenoidal claims stand
  (Round-2): the law carries exactly one scenario-level effective z₀.
- The analytic gradient is term-wise exact (central differences converge to
  it at measured order 2).
- Mode draws are a pure function of (seed, kernel, tile): independent of
  mode count and of other fields built in between for the original wall basis.
- Homogeneous Mann modes instead use an explicit finite-band quadrature whose
  geometry depends on its mode count. Every velocity polarization is transverse
  to its wavevector, so its instantaneous analytic divergence vanishes without
  a wall projection. Unresolved energy is omitted, never rescaled away.

- `ou::OuMode::from_correlation_time(sigma_st, tau)` — exact-discrete
  ρ = e^(−Δt/τ), σ_innov = σ_st·√(1−ρ²) (the stationary law is a FIXED
  POINT of the transition, battery-proven).
- `ou::StationaryOuPath` — StationaryOuPathV1: exact stationary draw at
  the FIXED anchor tick −3840; philox counter-addressed innovations
  (pure function of (seed, OU_KERNEL, mode, tick)); sequential
  `advance_to` (never backwards, span-capped), `checkpoint`/`resume`
  with BITWISE reconstruction equivalence; every preroll window is a
  suffix of ONE physical path.

## Error model

Typed `Refusal { code, message, ranked_repairs }`. Codes:
`non-finite-input`, `z0-outside-domain`, `displacement-invalid`,
`reference-height-invalid`, `reference-speed-invalid`,
`mode-count-invalid`, `turbulence-params-invalid`, `below-surface-query`;
OU layer: `ou-params-invalid`, `ou-advance-backwards`,
`ou-advance-span-exceeded`, `ou-checkpoint-invalid`.
Caps tested at cap AND cap+1.
Evolving-field admission also reports `evolving-params-invalid`,
`evolving-checkpoint-invalid`, and `sample-count-mismatch`. The two
quadratures must share a clock and use distinct prescribed stream keys;
invalid advance requests and invalid batches leave state/output unchanged.

## Determinism class

Deterministic: det:: transcendentals, philox counter-addressed draws, no
global state; V-04a golden pinned under
`org.frankensim.fs-atmo.v04a-golden.v1` (golden-bump protocol).

## Cancellation behavior

Synchronous pure functions; nothing to cancel.

## Unsafe boundary

Workspace `deny(unsafe_code)`; no unsafe.

## Feature flags

None.

## Conformance tests

`tests/v04a_battery.rs` — V-04a EXECUTED: relative divergence < 1e-12 over
288 probe/tick combinations; bitwise wall parity; FD-vs-analytic gradient
Richardson order ≈ 2 with finest-step error < 1e-6; air-state consistency
(Re/q hand-recomputed; log law recovers the reference point; E1.8 rho with
provenance string); seed determinism bit-identical + counter-addressed
partition law; refusal caps; log-law identities (U(d+z₀) = 0 exactly,
monotone, instrument-height ordering); pinned golden. JSONL receipts.

`tests/ou_battery.rs` — E3.3b-i / V-17 core: exact-discrete relations vs
hand closed forms + the stationary fixed point; no-spin-up window
variances; the suffix property (pause/resume vs straight advance,
bitwise); checkpoint reconstruction equivalence (bitwise) + checkpoint
refusals; the 39-minute (284k-tick) recurrence battery in checkpointed
chunks; caps at cap AND cap+1; pinned OU golden.

`evolving::tests` — G0/G3/G5: all nine stationary stress components against
an independent curl calculation (2e-14 absolute tolerance), horizontal
homogeneity, advected OU autocovariance and stationary variance over 4096
independent paths (0.07 normalized tolerance), exact wall parity, analytic
divergence below 1e-14 on the fixture, gradient finite-difference order
2 ± 0.05, read-only query/batch equivalence, bitwise checkpoint/resume, and
velocity/gradient/air-property consistency.

`tests/fit_battery.rs` additionally checks the finite Γ = 0 limit, seed/target
binding, input refusal, and the unchanged nonzero-shear fit golden. The
existing per-truncation and bootstrap checks concern the fitted one-point
stresses, not full spectral/coherence agreement.

`mann_field::tests` checks all nine modal tensor components for Γ = 0, 3.9
and 10, independent isotropic radial-energy integration, sheared two-point
covariance at oblique separations over 4096 independent amplitude seeds,
advected OU covariance over 2048 independent seeds, analytic gradients and
divergence, absence of an invented wall condition, exact checkpoint/resume,
and model/state refusals. These validate the represented finite-band model;
they do not certify continuum quadrature convergence or measured wind spectra.

## No-claim boundaries

- The modal amplitude decay is an ANALYTIC-CONSTRUCTION placeholder: no
  statistical claim of any kind (spectra, coherence, TI) — those are
  E3.3b's Mann-class fit under V-04b1/V-04b2 with its own artifacts.
- `TurbulenceField` retains static per-realization amplitudes and pure mean
  advection. `EvolvingTurbulenceField` supplies the declared finite-mode OU
  temporal model; its correlation times are explicit model inputs. Its
  conditional stationarity and covariance tests do not validate atmospheric
  spectra, intermittency, nonlinear energy transfer, or the full Mann target.
  The existing Flyer lifecycle's reduced gust forcing is not changed by
  opting into the new fs-atmo API elsewhere.
- Neutral stability only (declared, per air-state-v1); other classes are
  future additive modes.
- No gust events, no thermals, no site heterogeneity (E3.3c territory).
- `DEC17_AIR` is the ensemble MEAN of the pre-registered distributions;
  H-case members draw their own scenario per the frozen registry.
