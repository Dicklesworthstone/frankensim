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
`fs_cutfem::CutFemError` for CutFEM build/solve failures.

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
