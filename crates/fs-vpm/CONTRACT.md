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
- Particle-strength-exchange viscosity and dynamically stretched 3-D
  filaments/rings remain deferred. The connected 3-D wake's prescribed
  geometry and reduced finite-core model are specified below.
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
- Retired from v1 (owner decision 2026-09-01, bead frankensim-rc-root-q61wp.25):
  FMM-accelerated Biot–Savart, particle-strength-exchange viscosity, remeshing,
  and the BEM+VPM hybrid are NOT scheduled; the promise is preserved in deferred
  epic frankensim-rc-root-q61wp.43 and un-parks only behind a consumer journey.

## Finite-core connected 3-D wake

The connected filament, aggregation, and image lanes use the Rosenhead--Moore
kernel `dl cross r / [4 pi (|r|^2 + a^2)^(3/2)]`. Straight-segment integration
is analytic. Its unit-tangent form avoids comparing incompatible length powers,
and its exterior endpoint difference is rationalized to retain the small field
near the segment's extended axis. There is no absolute length cutoff. Exactly
axial or degenerate segments return zero, including the singular zero-core
convention. Geometry and core radius scaled by s give velocity scaled by 1/s
at fixed circulation.

`WakeCoreEvolutionMode::CoreSpreading` admits finite effective kinematic
viscosity in `[0, MAX_NU_EFF]`, in m^2/s. For a fixed declared viscosity,
`a^2 = a0^2 + 4 nu t`; zero viscosity recovers Frozen velocities. This is a
reduced algebraic core-spreading model, not exact Gaussian Lamb--Oseen
diffusion, particle-strength exchange, or resolved Navier--Stokes diffusion.
The core width is a regularization scale; it is not a finite physical second
moment of the algebraic profile. Coarsening also contributes a geometric
smoothing scale without changing that interpretation.

Each shed row retains elapsed `age_ticks`, starting at zero and increasing
only on a successful shed. Coarsening preserves circulation-weighted age
alongside position and core width. Relative weights are normalized without a
dimensional circulation floor, so changing circulation units cannot alter
geometry, core width, or age. Aggregation uses retained age and its declared
`shed_dt_s`; removing rows never resets the ages of those left near. Public
hybrid evaluation and image snapshots apply the same viscosity and shedding
interval to the remaining near wake and the aggregated sources.

Every far cell retains each original source segment and its own squared core.
Exact evaluation, near-field fallback, and prune audits use those individual
cores. The far expansion uses a circulation-length-weighted common core,
exact line monopole/dipole/quadrupole moments (including the intrinsic
line quadrupole), and the leading core-variation moment. Spatial/core mixed
terms remain in the truncation. Probes inside the geometric bounding sphere,
or with maximum squared-core deviation greater than
`(|r|^2 + common_core^2)/16`, use exact retained sources. This is an expansion
domain guard, not a certified total-velocity error bound. Separated mixed-core
cells retain constant work per cell rather than always scanning their sources.
The existing singular contribution bound remains conservative for regularized
sources; pruning still independently audits the exact contribution.

`FarField::advance_core(dt)` adds continued diffusion to already aggregated
cells in O(number of cells), without scanning their segments or rebuilding
moments. A common offset and cached maximum base core provide transactional
finite checks; zero time, Frozen, and zero viscosity preserve state exactly.
The method does not convect geometry. Callers keep near shedding and far
diffusion on the same elapsed-time schedule, and advance existing far cells
before appending newly aggregated cells to avoid aging those new cells twice.
Images are evaluation snapshots and must be rebuilt after geometry, diffusion,
or mode changes. Mirroring preserves source cores, diffusion offsets and the
axial-vector core moment.

The filament and far-field digest domains are v2 because retained age and
individual source/core assignment now affect physical evaluation. The updated
replay goldens cover that intentional identity and kernel change.
`tests/viscosity3d.rs` checks independent line quadrature, the analytic long
filament, length scaling, thin exterior fields, exact/near/far/image agreement,
mixed-core convergence and fallback, chronology through conversion, continued
diffusion and transactional refusals. These checks do not establish free-wake
stretching, a viscous particle solver, cross-ISA bit stability, or a wall-clock
performance guarantee.
