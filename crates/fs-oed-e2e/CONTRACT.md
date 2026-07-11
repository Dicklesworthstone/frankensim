# CONTRACT: fs-oed-e2e

SensorForge — optimal experimental design that knows when to stop. Layer L4
(ASCENT).

## Purpose and layer

Composes `fs-assimilate` (instrument-bound Kalman candidates), `fs-voi` (EVPI +
recommend), `fs-toleralloc` (first-order budget), `fs-evidence` (Estimated
lineage), and `fs-blake3` (bounded canonical identities). Dependencies point
downward.

## Public types and semantics

- `Candidate::new(...) -> Result<Candidate, CandidateError>` checks and then
  seals the name, finite truth/prior, non-negative prior variance, and positive
  finite sensor noise/cost. Read-only accessors expose the declaration.
- `run_campaign(&[Candidate], threshold, max_sensors) -> Result<OedReport,
  OedError>` validates the threshold, unique candidate identities, and explicit
  synchronous work caps before greedily placing VoI-chosen sensors.
- `demo_candidates() -> Result<Vec<Candidate>, CandidateError>` builds the four
  checked candidates whose close top two drive the decision.
- `OedReport` includes the native EVPI trace, posterior summaries, one
  instrument-bound `assimilation_color` per placement, and input-bound final
  variance/EVPI colors so consumers do not need to transcribe the loop.

## Invariants

- Sensors land on decision-relevant candidates (the close contenders), never on
  a dominated one; each placement shrinks total posterior variance.
- The campaign evaluates STOP before its placement cap, including when
  `max_sensors == 0`. `decision_robust` is true only when final modeled EVPI is
  at or below the checked threshold; a no-useful-action stop above threshold is
  not mislabeled robust.
- The posterior variance is `Estimated`: the scalar Kalman formula is exact for
  its declared linear-Gaussian model, but neither floating-point roundoff nor
  model-form assumptions carry an interval certificate. The EVPI stop is
  `Estimated`. Both bounded estimator identities commit to the complete ordered
  candidate declarations, threshold, placement cap, realized placement and
  posterior sequences, and canonical assimilation colors.
- Zero total prior variance has fractional reduction `0.0`, not NaN. A
  zero-sensitivity candidate receives `+infinity` allocation, the exact
  unconstrained first-order optimum; positive-sensitivity allocations must be
  positive and finite.
- Deterministic (readings hit each truth; Kalman variance is observation-free).

## Error model

No documented input panic. Candidate and campaign rejection is structured as
`CandidateError` / `OedError`; lower-layer assimilation failures retain their
`AssimError` source. Resource admission caps candidates, placements, and their
quadratic `candidate_count^2 * max_sensors` action-design evaluation product
before campaign allocation or iteration (`recommend` scores every action
against the design set at each placement).

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/oed.rs`: decision-relevant placement and stopping; initial STOP at a zero
placement cap; full-report determinism; adversarial candidate/campaign input
rejection; zero-variance behavior; unmeasured-input evidence binding; and
instrument-bound assimilation lineage.

## No-claim boundaries

Scalar (diagonal) beliefs and a Gaussian EVPI model; the greedy VoI policy is
near-optimal, not the optimal information-gathering policy; the precision budget
uses a prior-std sensitivity proxy. `decision_robust` is a modeled EVPI
criterion, not a physical decision certificate or independent validation.
