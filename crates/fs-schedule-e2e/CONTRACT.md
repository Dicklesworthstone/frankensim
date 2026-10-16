# CONTRACT: fs-schedule-e2e

CampaignSchedule — verified makespan arithmetic plus advisory
value-of-information scheduling. Layer L6 (HELM).

## Purpose and layer

Composes `fs-tropical` (max-plus critical path), `fs-voi` (EVPI + recommend),
and `fs-evidence` (Verified makespan enclosure vs Estimated information value).
Deps downward.

## Public types and semantics

- `Study { name, latency, deps }` — a node in the precedence DAG.
- `run_campaign(&[Study], &[DesignEstimate], &[Action], stop_threshold) ->
  Result<ScheduleReport, ScheduleError>` — validates bounded canonical inputs,
  then computes the makespan/bottleneck/nominal slack (tropical) and the
  EVPI/leading design/flip risk/typed disposition (VoI).

## Invariants

- The makespan is admitted through a directed-rounding enclosure only after
  finite non-negative latency, graph, work-cap, and overflow checks →
  `Verified`. A bottleneck is named only when interval bounds prove a unique
  critical path. Positive nominal slack remains scheduling guidance, not a
  certified deferability claim.
- The recommendation is `Act: <study>` when EVPI exceeds the stop threshold and a
  study has positive value-per-cost; `RobustStop` only when EVPI is at or below
  the threshold; and `NoEffectiveAction` when ambiguity remains but the menu is
  deficient.
- Recommendation evidence is explicitly `Estimated`; it is never promoted by
  the tropical makespan enclosure.
- Deterministic (fixed spec; no RNG).

## Error model

Structured `ScheduleError`; empty/oversized campaigns, malformed or duplicate
ASCII-graphic names, non-finite/negative inputs, ambiguous action targets,
graph defects, cycles, Cartesian decision-work excess, and numerical overflow
refuse before positive evidence. No caller input causes a panic.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/schedule.rs`: the makespan enclosure contains 13, the unique bottleneck is
windtunnel-A, and the decision is explicitly Estimated; robust stop and
no-effective-action are distinct; deterministic malformed/cyclic/non-finite,
identity, and exact Cartesian-work boundaries fail closed.

## No-claim boundaries

The DAG latencies are supplied constants (a real cost model comes from `fs-plan`);
EVPI uses `fs-voi`'s Gaussian decision model; the recommendation is a single
next-step, not a full optimal information-gathering policy.
