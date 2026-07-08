# CONTRACT: fs-vpm

Vortex particle method (2-D inviscid core): vorticity-native wake dynamics.

## Purpose and layer

Layer L3 (FLUX). No dependencies — pure Rust.

## Public types and semantics

- `VortexParticle { pos, circulation }` — a point vortex carrying circulation
  `Γ`.
- `induced_velocity(&particles, point, core)` — the desingularized 2-D
  Biot–Savart velocity `Σ (Γⱼ/2π)·perp(r)/(|r|²+core²)`; `core = 0` skips
  coincident particles.
- `advect(&particles, dt, core)` / `simulate(&particles, dt, steps, core)` —
  RK4 advection of the N-body vortex system (circulations invariant).
- `total_circulation(&particles)` — `Σ Γᵢ` (conserved).
- `vorticity_centroid(&particles)` — `Σ Γᵢxᵢ / Σ Γᵢ` (invariant), `None` when
  the total circulation is zero.

## Invariants

- A single vortex `Γ` induces a purely TANGENTIAL field of magnitude `Γ/(2πr)`
  (`u·r = 0`) and does not advect itself.
- A counter-rotating PAIR separated by `d` self-propels in a straight line at
  exactly `Γ/(2πd)`, preserving its separation (the 2-D analog of a vortex
  ring's translation).
- Total circulation is conserved to roundoff; a co-rotating pair conserves its
  vorticity centroid.
- The desingularized kernel keeps the velocity finite at/near a particle.
- Deterministic (fixed particle order, no RNG).

## Error model

Total functions; no panics (`vorticity_centroid` returns `None` on zero total
circulation).

## Determinism class

Fully deterministic: induction and RK4 advection are pure functions.

## Cancellation behavior

None here; the production solver polls `Cx` at step boundaries.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/vpm.rs` (7 cases): a single vortex induces the analytic tangential field;
a single vortex does not move itself; a counter-rotating pair self-propels at
the analytic speed; a co-rotating pair conserves circulation + centroid; a
symmetric pair has no defined centroid; the desingularized core bounds the
velocity; determinism.

## No-claim boundaries

- v0 is the INVISCID 2-D core by DIRECT `O(N²)` Biot–Savart. The FMM /
  treecode acceleration for large particle counts (the whole point of the
  method at scale) is the fuller deliverable, staged.
- VISCOUS diffusion via PARTICLE STRENGTH EXCHANGE (Lamb–Oseen evolution) and
  3-D vortex filaments / rings with stretching are staged.
- The HYBRID BEM+VPM airfoil (pitching/plunging unsteady lift vs Theodorsen —
  the flapping credential) composes this core with fs-bem downstream.
- Spatial adaptivity (particle insertion/merging, remeshing) is a follow-on.
