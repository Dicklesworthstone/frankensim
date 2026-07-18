# CONTRACT: fs-vpm

Vortex particle method (2-D inviscid core): vorticity-native wake dynamics.

## Purpose and layer

Layer L3 (FLUX). Safe Rust with `fs-exec` for the production `Cx` and admitted
budget accountant. Tests use `fs-alloc` to construct scoped contexts.

## Public types and semantics

- `VortexParticle { pos, circulation }` — a point vortex carrying circulation
  `Γ`.
- `VpmBudget` — explicit caps for particle count, RK4 steps, exact attempted
  pair evaluations, and logical peak live payload bytes.
- `simulate_with_cx(cx, particles, dt, steps, core, budget)` — the checked,
  cancellation-aware production transaction. It returns a complete `VpmRun`
  with exact work and ambient-budget consumption, or a typed `VpmError` with no
  partial particle state.
- `induced_velocity(&particles, point, core)` — the desingularized 2-D
  Biot–Savart velocity `Σ (Γⱼ/2π)·perp(r)/(|r|²+core²)`; `core = 0` skips
  coincident particles.
- `advect(&particles, dt, core)` / `simulate(&particles, dt, steps, core)` —
  original unchecked compatibility helpers for analytic fixtures and existing
  callers. They have no explicit work/memory/cancellation authority and are not
  the production entry point.
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
- The checked path attempts exactly `4 S N²` source-target contributions for
  `N` particles and `S` steps, including coincident/self visits that the kernel
  skips numerically.
- Its exact logical-work charge is `N + S(4N² + 4N + 1)`: one combined input
  validation/copy visit per particle, four all-pairs sweeps, three stage
  displacements plus one final combine, and one private step commit.
- For `S > 0`, the reusable-buffer implementation admits a logical peak of
  `N(5P + 4V)` bytes, where `P = size_of::<VortexParticle>()` and
  `V = size_of::<[f64; 2]>()`; for `S = 0`, it admits `NP` bytes. This is
  operation-owned vector payload, not allocator-realized resident memory.
- Input order, target order, source order, and RK4 arithmetic association match
  the compatibility path. Valid admitted fixtures therefore have bitwise-equal
  results on the same build and ISA.

## Error model

The checked path validates finite particle coordinates/circulations and finite
`dt`; signed finite `dt` is allowed. `core` must be finite, nonnegative, and
have a finite square. Count/work/byte arithmetic is checked before allocation.
Local cap, ambient budget, allocation, and non-finite-intermediate failures are
typed. All owned vectors use fallible reservation, and an error publishes no
particle vector.

The compatibility helpers retain their original infallible signatures and
ordinary `Vec` allocation behavior. `vorticity_centroid` returns `None` on zero
total circulation.

## Determinism class

Same-build, same-ISA deterministic: fixed input/target/source order, fixed RK4
association, and no RNG or worker-count dependence. No cross-ISA bit-stability
claim is made.

## Cancellation behavior

`simulate_with_cx` admits the ambient cost quota against its complete checked
work plan, then checkpoints before allocation, before each chunk of at most
`VPM_WORK_CHECKPOINT_STRIDE = 256` logical units, at RK-stage and step
boundaries, and immediately before publication. Completed units are charged
before the next checkpoint. Cancellation/deadline/poll refusal latches through
`AdmittedBudget`; synchronous private buffers drain by stack unwinding, and no
partial state escapes.

The cancellation bound is in logical work units, not microseconds. The
compatibility `advect` and `simulate` helpers do not accept a `Cx` and make no
cancellation claim.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/vpm.rs` contains the seven original analytic/determinism fixtures plus G0
checked-path evidence for bitwise legacy equivalence, exact work/byte receipts,
zero-step identity, cap and arithmetic-overflow refusal, invalid inputs, and
non-finite intermediates; G4 evidence covers pre-cancellation and ambient
cost/poll refusal without partial publication.

## No-claim boundaries

- v0 is the INVISCID 2-D core by DIRECT `O(N²)` Biot–Savart. The FMM /
  treecode acceleration for large particle counts (the whole point of the
  method at scale) is the fuller deliverable, staged.
- VISCOUS diffusion via PARTICLE STRENGTH EXCHANGE (Lamb–Oseen evolution) and
  3-D vortex filaments / rings with stretching are staged.
- The HYBRID BEM+VPM airfoil (pitching/plunging unsteady lift vs Theodorsen —
  the flapping credential) composes this core with fs-bem downstream.
- Spatial adaptivity (particle insertion/merging, remeshing) is a follow-on.
- Logical live bytes exclude the borrowed input, allocator metadata, capacity
  rounding, allocator fragmentation, and process RSS. Admission cannot promise
  allocator success; reservation failure remains typed.
- No integration-error, convergence-order, stability/CFL, conservation
  certificate, wall-clock deadline-latency, resumable-partial-state,
  cross-ISA-bit-stability, or performance claim is made by the checked receipt.
- G4 cancellation evidence is synchronous and deterministic. It does not yet
  establish cancellation under multi-threaded signal races or a wall-time
  latency bound.
