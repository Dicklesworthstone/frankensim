# CONTRACT: fs-flowcert-e2e

FlowCert — a certified credibility map for a lattice-Boltzmann channel flow.
Layer L4 (ASCENT).

## Purpose and layer

Composes `fs-lbm` (LBM + analytic Poiseuille solution + scaling planner),
`fs-archive` (MAP-Elites), `fs-evidence` (Verified). Deps point downward.

## Public types and semantics

- `OperatingPoint { reynolds, ny, tau, viscosity, profile_error, accurate,
  regime_stable }`.
- `certify_point(reynolds, ny, u_lattice, max_steps, tol) -> OperatingPoint` —
  marches a channel to STEADY STATE (chunked, stopping when the profile
  stabilizes, capped at `max_steps`), compares to the analytic Poiseuille
  profile, and reads the scaling-planner regime certificate.
- `run_campaign(&reynolds, &resolutions, max_steps, tol) -> FlowReport` —
  illuminates the (Reynolds × resolution) atlas.
- `default_sweep()` — the default grid.

## Invariants

- ACCURACY REFLECTS THE REGIME, NOT THE BUDGET: each point is run to steady
  state (`converged`), so `profile_error` is the inherent O(1/ny²) discretization
  error — every converged point matches the analytic solution within `tol`.
- CREDIBILITY MAP: because accuracy is uniform, the differentiation is the REGIME
  certificate — low-Reynolds points sit in a `Verified` (comfortable `τ`-margin)
  regime, while near-`τ=½` points are flagged as risky even where they are
  accurate; `stable_fraction ∈ (0, 1)`.
- The whole-map color is `Verified` only if every point is accurate and stable.
- Deterministic (fixed chunked LBM march; no RNG).

## Error model

Panics only on an empty sweep.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch); production LBM would poll `Cx` per streaming sweep.

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/flowcert.rs` (4): the credibility map is illuminated, every point
converges + is accurate, and the regime certificate separates credible from
flagged points; a single low-Reynolds point is fully verified; a near-`τ=½` point
is flagged by regime even though it is accurate; determinism.

## No-claim boundaries

A 1-D Poiseuille channel (no scalar transport, so no mixing metric); accuracy is
a manufactured-solution comparison, not a grid-convergence order; the regime
certificate is the `fs-lbm` low-Mach/`τ`-margin heuristic, not a full stability
spectrum.
