# CONTRACT: fs-couple

Multiphysics composition through port-Hamiltonian Dirac structures: partitioned
coupling that is passive by construction.

## Purpose and layer

Layer L3 (multiphysics coupling). No dependencies — pure Rust.

## Public types and semantics

- `PortKind` (mechanical force/velocity, fluid pressure/flux, thermal
  temperature/entropy); `Port { effort, flow, kind }` with `power` = effort ×
  flow and `conjugate_to` (same physical type — the composition-time type
  discipline).
- `interconnect(kind_a, kind_b, effort, flow) -> Result<Interconnection,
  CoupleError>` — a Dirac structure (shared effort, opposite flow) whose
  `interface_power` is `0` exactly (power-conserving by construction); refuses
  incompatible ports. `interface_power(&[Port])` = `Σ effort·flow`.
- `EnergyAudit` — `record`, `max_generation`, `is_passive(tol)`: the
  passivity claim is MEASURED every exchange (a nonzero balance is a bug alarm).
- `AitkenRelaxation::new(omega_init, omega_max)` + `next_omega(residual)` — the
  scalar Δ² dynamic relaxation factor, magnitude-capped.
- `iterate_fixed_relaxation` / `iterate_aitken` — the added-mass interface
  fixed point under fixed vs Aitken relaxation → `FsiResult { converged, steps,
  solution, final_residual }`.

## Invariants

- The Dirac interconnection conserves interface power EXACTLY (to roundoff) —
  the G0 law; incompatible ports are refused.
- The energy audit reports non-passive exactly when some exchange generates
  power above `tol`.
- On the added-mass fixture (`μ ≥ 1`): naive staggering (`ω = 1`) diverges while
  Aitken-relaxed coupling converges to `x* = c/(1+μ)`; Aitken never takes more
  steps than a stable fixed under-relaxation.

## Error model

Structured `CoupleError`; no panics.

## Determinism class

Fully deterministic: interconnection, audit, and the iterations are pure
functions of their inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/couple.rs` (8 cases): power-conjugate ports; exact interface power
conservation + incompatible-port refusal; the energy audit measures passivity
+ alarms on generation; the Aitken Δ² factor; naive staggering diverges where
Aitken stays stable (the added-mass claim); Aitken accelerates over a stable
fixed relaxation; light added mass converges even naively; determinism.

## No-claim boundaries

- The FSI fixture is the classic LINEARIZED added-mass interface map
  (`H(x) = −μx + c`) — enough to reproduce the instability and its fix; a full
  nonlinear FSI solve over real fluid/structure subsystems is the consumer.
- `AitkenRelaxation` is the scalar Δ² relaxer; the vector INTERFACE
  QUASI-NEWTON (IQN-ILS) accelerator and MULTIRATE co-simulation are staged.
- `PortKind` encodes the conjugate physical type; full `Qty`-dimensioned
  effort/flow conjugacy checking (fs-qty) and the categorical composition API
  over general port-Hamiltonian systems are staged with the interface here.
- The energy audit's balances are supplied by the caller each exchange; wiring
  them onto the ledger is the coupling driver's integration.
