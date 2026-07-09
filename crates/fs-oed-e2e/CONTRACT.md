# CONTRACT: fs-oed-e2e

SensorForge — optimal experimental design that knows when to stop. Layer L4
(ASCENT).

## Purpose and layer

Composes `fs-assimilate` (Kalman fusion), `fs-voi` (EVPI + recommend),
`fs-toleralloc` (budget), `fs-evidence` (Verified/Estimated). Deps downward.

## Public types and semantics

- `Candidate { name, truth, prior_mean, prior_var, sensor_noise, sensor_cost }`.
- `run_campaign(&[Candidate], threshold, max_sensors) -> OedReport` — greedily
  places VoI-chosen Kalman-fused sensors, stops when EVPI ≤ threshold, and
  allocates the precision budget.
- `demo_candidates()` — four candidates whose close top two drive the decision.

## Invariants

- Sensors land on decision-relevant candidates (the close contenders), never on
  a dominated one; each placement shrinks total posterior variance.
- The campaign STOPS the instant EVPI ≤ threshold (`decision_robust`), and a
  clear winner needs zero sensors.
- The posterior variance is `Verified` (exact scalar Kalman); the EVPI stop is
  `Estimated`. A cost-optimal precision allocation is produced for every
  candidate.
- Deterministic (readings hit each truth; Kalman variance is observation-free).

## Error model

Panics only on an empty candidate list.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/oed.rs` (3): sensors target the decision + stop when robust (A chosen,
D never measured, EVPI falls below threshold); a clear winner needs no sensors;
determinism.

## No-claim boundaries

Scalar (diagonal) beliefs and a Gaussian EVPI model; the greedy VoI policy is
near-optimal, not the optimal information-gathering policy; the precision budget
uses a prior-std sensitivity proxy.
