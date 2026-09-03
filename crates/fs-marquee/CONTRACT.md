# fs-marquee CONTRACT

## Purpose and layer

Layer: L6 (HELM/integration). `fs-marquee` names the P2 marquee study lane:
raw SDF geometry through CutFEM physics, DWR evidence, ledger records, and
renderable artifacts. The default build exposes the smoke-tier study runner for
the raw-SDF/CutFEM/DWR slice. The `marquee` feature remains available as an
explicit opt-out boundary; the full-resolution nightly golden lane remains a
no-claim boundary.

## Public types and semantics

- `MarqueeStatus`: status of the lane. `Disabled` means the optional `marquee`
  feature was disabled. `SmokeRunnerAvailable` is the default and means the
  smoke-tier runner API is available.
- `status()`: deterministic status query derived only from Cargo feature
  configuration.
- `scope_summary()`: static diagnostic text for agents, ledgers, and reports.
- `VERSION`: crate version for provenance stamping.
- With `marquee`: `study::{PlateWithHoles, StudyConfig, StudyReport,
  IterRecord, run_study}`. The runner performs a deterministic projected
  radius optimization over circular cooling holes, records per-iteration
  compliance/certificate fields, and returns a replay hash for the smoke trace.

The default build exposes the smoke-tier simulation entrypoint. It performs
in-process CutFEM solves and does not mutate ledgers or the filesystem.

- `study` module (the smoke-tier runner, bead mye.1): `PlateWithHoles`
  (an EXACT parametric SDF with a certified box enclosure — the CutSdf
  containment law), `run_study` (CutFEM state solves, the self-adjoint
  compliance shape gradient `dJ/dr = −∮(∂u/∂n)²` — sign CAUGHT by the
  FD falsifier during development, the drill earning its keep — with an
  area-budget rescale projection), `IterRecord` (compliance, area,
  gradient, the three certificate components, including an algebraic term
  formed only from CutFEM's typed recomputed-Euclidean residual accessor,
  composed color, solver
  iterations), `StudyReport.trace_hash` (the G5 replay witness).

## Invariants

1. The default build makes the smoke-tier marquee study available.
2. Disabling default features prevents the smoke runner from being exposed.
3. Runner inputs are admitted before CutFEM work starts: at least one hole,
   matching center/radius lengths, finite unit-plate centers, positive finite
   radii, finite area target in `(0, 1)`, nonnegative finite step size, and
   finite positive radius bounds.
4. The exposed runner is deterministic for a fixed source tree and machine.

## Error model

Default status queries are infallible. With `marquee`, invalid study inputs
panic during admission before solver work starts. Valid study runs return
`fs_cutfem::CutFemError` for CutFEM build/solve failures. Shape-gradient
boundary probes read the solved field only through the canonical
fail-closed `Space::sample_scalar` (bead ay40): missing or non-finite
active nodal evidence surfaces as `InvalidFemInput` instead of a
plausible zero, and the only zero read without evidence is a
certified-Outside classification, mapped explicitly to the homogeneous
Dirichlet exterior value at the use site. Marquee additionally refuses
certificate composition if a future CutFEM solver returns anything other than
a recomputed Euclidean residual.

## Determinism class

D0 for the status API. The smoke runner is deterministic for fixed
inputs and code, but it is not yet a cross-ISA golden-proofed lane.

## Cancellation behavior

The default smoke runner is synchronous and currently has no explicit `Cx`
cancellation polling; production runner cancellation remains a no-claim
boundary.

## Unsafe boundary

No unsafe code.

## Feature flags

- `marquee`: enables the smoke-tier raw-SDF/CutFEM/DWR study runner. It is
  enabled by default and may be explicitly disabled for a status-only build.

## Conformance tests

Unit tests check version stamping, default smoke-runner admission,
feature-derived status, and the explicit nightly-golden no-claim boundary. With
`marquee`, tests also check that invalid runner inputs are rejected before
solver work starts, and execute the six marquee falsifiers:
1. `mq_006_falsifier_rung_climb`: Rung climb on the final optimized design at
   level 5 vs level 4 within certified DWR error band.
2. `mq_007_falsifier_cross_representation_solid`: Cross-code body-fitted P1
   triangular linear elasticity in `fs-solid` validates physical compliance
   against CutFEM.
3. `mq_008_falsifier_adjoint_fd_gate_at_stages`: Informative-direction FD gate
   passes at iterates 1, N/2, and N, and rejects sign-flipped adjoint.
4. `mq_009_falsifier_objective_sensitivity_twin`: Volume fraction and geometric
   perturbation twins track respective objectives.
5. `mq_010_falsifier_replay_and_checkpoint_resume`: Replays bit-exact from seed
   and N/2 checkpoint resumes to identical endpoint.
6. `mq_011_falsifier_mutation_proof_monotonicity`: Sign-flipped ascent fails
   Armijo monotonicity check and is rejected.

## No-claim boundaries

- No sphere-traced render output is shipped here.
- No replayable golden ledger is shipped here.
- No full-resolution/nightly golden study lane is shipped here.
- No filesystem/ledger mutation is performed by the smoke runner.
- No performance, convergence beyond the smoke tests, physical-validity beyond
  the estimated DWR/algebraic fields, or rendering-quality claims attach to
  this crate until the full runner and its Gauntlet evidence land.

## No-claim boundaries (study)

- SMOKE tier only: level-4/5 quadtrees, 8-step budgets — the
  full-resolution nightly golden lane and both-ISA runs are the
  remaining P2 exit work, not claims here.
- The composed certificate's headline color is ESTIMATED (DWR constants
  and the conversion from a recomputed Euclidean residual to a goal-error
  contribution are estimates; the recurrence residual itself is never used.
  The refined-reference check passes within a documented 4x effectivity band).
  Equilibrated 2-D
  brackets would upgrade it to Verified — future work.
- The FrankenScript-IR front end and the fs-report notebook are the
  gp3.10/fs-ir integration seams; the study exposes the runner they
  will drive.
- Thermal (Poisson) compliance, not elasticity: the canonical heat-sink
  layout study — the elasticity marquee follows the same seam once
  CutFEM elasticity lands.
