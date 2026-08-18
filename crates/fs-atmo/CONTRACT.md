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

## Error model

Typed `Refusal { code, message, ranked_repairs }`. Codes:
`non-finite-input`, `z0-outside-domain`, `displacement-invalid`,
`reference-height-invalid`, `reference-speed-invalid`,
`mode-count-invalid`, `turbulence-params-invalid`, `below-surface-query`.
Caps tested at cap AND cap+1.

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

## No-claim boundaries

- The modal amplitude decay is an ANALYTIC-CONSTRUCTION placeholder: no
  statistical claim of any kind (spectra, coherence, TI) — those are
  E3.3b's Mann-class fit under V-04b1/V-04b2 with its own artifacts.
- Static per-realization amplitudes: the OU time evolution (and hence any
  temporal-statistics claim) is E3.3b; the phase clock here is pure mean
  advection.
- Neutral stability only (declared, per air-state-v1); other classes are
  future additive modes.
- No gust events, no thermals, no site heterogeneity (E3.3c territory).
- `DEC17_AIR` is the ensemble MEAN of the pre-registered distributions;
  H-case members draw their own scenario per the frozen registry.
