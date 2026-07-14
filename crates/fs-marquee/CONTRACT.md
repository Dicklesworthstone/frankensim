# fs-marquee CONTRACT

## Purpose and layer

Layer: L6 (HELM/integration). `fs-marquee` names the P2 marquee study lane:
raw SDF geometry through CutFEM physics, DWR evidence, ledger records, and
renderable artifacts. The default build remains an admission/status shell. With
the `marquee` feature enabled, the crate exposes a smoke-tier study runner for
the raw-SDF/CutFEM/DWR slice; the full-resolution nightly golden lane remains a
no-claim boundary.

## Public types and semantics

- `MarqueeStatus`: status of the lane. `Disabled` means the `marquee` feature
  is off. `SmokeRunnerAvailable` means the feature flag is enabled and the
  smoke-tier runner API is available.
- `status()`: deterministic status query derived only from Cargo feature
  configuration.
- `scope_summary()`: static diagnostic text for agents, ledgers, and reports.
- `VERSION`: crate version for provenance stamping.
- With `marquee`: `study::{PlateWithHoles, StudyConfig, StudyReport,
  IterRecord, run_study}`. The runner performs a deterministic projected
  radius optimization over circular cooling holes, records per-iteration
  compliance/certificate fields, and returns a replay hash for the smoke trace.

The default build exposes no simulation entrypoint. The feature-gated smoke
runner performs in-process CutFEM solves and does not mutate ledgers or the
filesystem.

- `study` module (the smoke-tier runner, bead mye.1): `PlateWithHoles`
  (an EXACT parametric SDF with a certified box enclosure — the CutSdf
  containment law), `run_study` (CutFEM state solves, the self-adjoint
  compliance shape gradient `dJ/dr = −∮(∂u/∂n)²` — sign CAUGHT by the
  FD falsifier during development, the drill earning its keep — with an
  area-budget rescale projection), `IterRecord` (compliance, area,
  gradient, the three certificate components, composed color, solver
  iterations), `StudyReport.trace_hash` (the G5 replay witness).

## Invariants

1. The default build cannot accidentally execute a marquee study.
2. Enabling the `marquee` feature is required before the smoke runner is
   available.
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
Dirichlet exterior value at the use site.

## Determinism class

D0 for the default status API. The smoke runner is deterministic for fixed
inputs and code, but it is not yet a cross-ISA golden-proofed lane.

## Cancellation behavior

No long-running work exists in the default build. The feature-gated smoke
runner is synchronous and currently has no explicit `Cx` cancellation polling;
production runner cancellation remains a no-claim boundary.

## Unsafe boundary

No unsafe code.

## Feature flags

- `marquee`: frontier gate for the smoke-tier raw-SDF/CutFEM/DWR study runner.
  The default build remains status-only.

## Conformance tests

Unit tests check version stamping, feature-derived status, and the explicit
nightly-golden no-claim boundary. With `marquee`, tests also check that invalid
runner inputs are rejected before solver work starts.

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
  and the CG-residual proxy are estimates; the refined-reference check
  passes within a documented 4x effectivity band). Equilibrated 2-D
  brackets would upgrade it to Verified — future work.
- The FrankenScript-IR front end and the fs-report notebook are the
  gp3.10/fs-ir integration seams; the study exposes the runner they
  will drive.
- Thermal (Poisson) compliance, not elasticity: the canonical heat-sink
  layout study — the elasticity marquee follows the same seam once
  CutFEM elasticity lands.
