# fs-marquee CONTRACT

## Purpose and layer

Layer: L6 (HELM/integration). `fs-marquee` names the P2 marquee study lane:
raw SDF geometry through CutFEM physics, DWR evidence, ledger records, and
renderable artifacts. The crate is an admission/status shell only; it keeps the
workspace member valid and the frontier feature gate explicit while the actual
runner remains outside the shipped default path.

## Public types and semantics

- `MarqueeStatus`: status of the lane. `Disabled` means the `marquee` feature
  is off. `FeatureEnabledNoRunner` means the feature flag is enabled but the
  crate still exposes only the status surface.
- `status()`: deterministic status query derived only from Cargo feature
  configuration.
- `scope_summary()`: static diagnostic text for agents, ledgers, and reports.
- `VERSION`: crate version for provenance stamping.

No function in this crate launches simulation, optimization, rendering, ledger
mutation, or filesystem I/O.

## Invariants

1. The default build cannot accidentally execute a marquee study.
2. Enabling the `marquee` feature changes the status value only; it does not
   promote frontier functionality into the default path.
3. The crate is deterministic and side-effect free.

## Error model

No fallible operations are exposed. Invalid or unavailable runner behavior is
represented by `MarqueeStatus`, not by panics or silent work.

## Determinism class

D0 for the exposed API: all outputs are compile-time constants or Cargo-feature
derived constants.

## Cancellation behavior

No long-running work exists in this crate. Future runners must execute under the
project cancellation discipline and document their bounded polling points before
this no-claim boundary can be lifted.

## Unsafe boundary

No unsafe code.

## Feature flags

- `marquee`: frontier gate for the future end-to-end study. In the current
  crate it only changes `status()` from `Disabled` to `FeatureEnabledNoRunner`.

## Conformance tests

Unit tests check version stamping, feature-derived status, and the explicit
no-runner scope text.

## No-claim boundaries

- No raw-SDF-to-CutFEM optimization runner is shipped here.
- No DWR certificate composition is shipped here.
- No sphere-traced render output is shipped here.
- No replayable golden ledger is shipped here.
- No performance, convergence, physical-validity, or rendering-quality claims
  attach to this crate until the runner and its Gauntlet evidence land.
