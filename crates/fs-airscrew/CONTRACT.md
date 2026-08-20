# fs-airscrew — CONTRACT

Bead: frankensim-wf-root-guzez.5.12.1 (E4.5-i, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN §5.3 (ROUND 6 steady state). Evidence:
prop-geometry-v1 (E1.6) under the frozen registry.

## Purpose and layer

L3 (sibling of fs-wing; NEITHER depends on the other — fs-flyer owns the
coupling). E4.5-i ships the BEMT rotor kernel in the w-formulation (valid
through J = 0), Prandtl tip/root loss, the warm-started bounded station
solve with per-station convergence receipts, the declared engine torque
curve, and the rotor spin-up step. E4.5-ii adds fs-flyer's
CoupledPropAirframeStep (Aitken candidate A).

## Public types and semantics

- `BladeStation { r_over_r, chord_m, beta_rad }`, `Rotor { radius,
  n_blades, camber, stations }` — provenance lives with the caller
  (prop-geometry-v1 rules).
- `bemt_solve(rotor, rho, V_axial, omega) -> BemtSolution { thrust,
  torque, ct, cq, j, station receipts }` — CT/CQ/J in the
  prop-geometry-v1 conventions.
- `engine_torque_at_prop_nm(omega)` — 12 hp @ 1025 engine rpm through the
  23:8 chain, flat-torque-then-power-limited (DISCLOSED approximation).
- `rotor_spinup_step(I, omega, Q_eng, Q_prop, dt)`.

## Invariants

- The w-formulation momentum balance (V+w)·w = dT/(4πrρF) stays valid at
  V = 0 (static bench states are first-class, not a special case).
- Station marching warm-starts from the inboard neighbor; nonconvergence
  is a TYPED refusal, never a clamped answer.
- Prandtl F ∈ (0, 1], falling toward the tip.

## Error model

`rotor-invalid` (stations cap AND cap+1, ordering, ranges),
`operating-point-invalid` (V < 0, ω ≤ 0), `station-did-not-converge`
(names the station), `rotor-dynamics-invalid`.

## Determinism class

Deterministic (det:: transcendentals, fixed iteration order); golden
pinned under `org.frankensim.fs-airscrew.v03-golden.v1`.

## Cancellation behavior

Synchronous pure functions.

## Unsafe boundary

Workspace `deny(unsafe_code)`; none.

## Feature flags

None.

## Conformance tests

`tests/bemt_battery.rs` — V-03 core: the E1.6 HOLDOUT static anchor
(1903 reconstruction at 350 rpm vs the 285 N bench, factor-2 trend band);
CT strictly falling over J 0.1–0.9; η rising to a 0.5–0.95 peak in the
Wright J-class then falling; Prandtl limits + tip bite; spin-up from rest
to the 250–500 rpm equilibrium class against the declared engine curve;
caps at cap AND cap+1; golden at the rail state (J ≈ 0.71).

## No-claim boundaries

- The 1903 radial table is an ESTIMATED reconstruction (prop-geometry-v1
  1903-absence rule); quantitative CT/CQ validation awaits real tables.
- Swirl enters via the section angle only (small-swirl DISCLOSED);
  dynamic inflow, blade unsteadiness, and inflow harmonics are E4.5-ii+.
- The engine curve is a declared flat-torque approximation; the real
  torque-vs-rpm curve is load-bearing and remains an E1 follow-on.
- No airframe coupling here (fs-flyer owns CoupledPropAirframeStep).
