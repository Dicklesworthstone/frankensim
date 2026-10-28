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
amplitude evolution (sequential checkpointed state) and the Mann-class
spectral-tensor fit; E3.3c the optional FetchAdjustedMassConsistent mode.

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

## Invariants

- Solenoidality: u is a curl — ∇·u ≡ 0 analytically; the analytic-gradient
  trace cancels to machine precision (V-04a gate, relative bound).
- Wall parity: sin(k_h·0) = 0 makes the vertical component vanish at h = 0
  IDENTICALLY (bitwise zero, battery-asserted).
- Pointwise z₀ insertion is FORBIDDEN while solenoidal claims stand
  (Round-2): the law carries exactly one scenario-level effective z₀.
- The analytic gradient is term-wise exact (central differences converge to
  it at measured order 2).
- Mode draws are a pure function of (seed, kernel, tile): independent of
  mode count and of other fields built in between.

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
