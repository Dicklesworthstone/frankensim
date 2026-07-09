# CONTRACT: fs-schedule-e2e

CampaignSchedule — certified scheduling of a design campaign, driven by value of
information. Layer L6 (HELM).

## Purpose and layer

Composes `fs-tropical` (max-plus critical path), `fs-voi` (EVPI + recommend),
`fs-evidence` (Verified makespan vs Estimated information value). Deps downward.

## Public types and semantics

- `Study { name, latency, deps }` — a node in the precedence DAG.
- `run_campaign(&[Study], &[DesignEstimate], &[Action], stop_threshold) ->
  ScheduleReport` — computes the makespan/bottleneck/slack (tropical) and the
  EVPI/leading design/flip risk/recommendation (VoI).

## Invariants

- The makespan is the EXACT tropical critical-path length → `Verified`
  (`lo == hi`); the bottleneck is the highest-latency critical study; slack
  studies are safe to defer.
- The recommendation is `Act: <study>` when EVPI exceeds the stop threshold and a
  study has positive value-per-cost, else `Stop` (decision already robust).
- Deterministic (fixed spec; no RNG).

## Error model

Panics only on an empty study list; a cyclic DAG surfaces via `expect`.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/schedule.rs` (3): the schedule + decision are both certified (makespan 13,
bottleneck windtunnel-A, Act recommendation); a robust decision recommends Stop;
determinism.

## No-claim boundaries

The DAG latencies are supplied constants (a real cost model comes from `fs-plan`);
EVPI uses `fs-voi`'s Gaussian decision model; the recommendation is a single
next-step, not a full optimal information-gathering policy.
