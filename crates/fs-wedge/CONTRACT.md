# CONTRACT: fs-wedge

Go-to-market wedge selection as data (plan addendum, Proposal 7): the
historical three-vertical ranking, its superseded judgment scores, measured
decision inputs, and the cycle-time kill criterion.

## Purpose and layer

Layer UTIL (pure data + audit; no dependencies). This gates NOTHING technical —
it constrains vertical-specific kernel work and records the commercial bet.

## Public types and semantics

- `WEDGE_DOCTRINE` — the load-bearing NEGATIVE rule ("do not sell against peak
  single-physics fidelity").
- `WedgeCriterion` (4: kernel maturity, iteration pain, quantifiable ROI, low
  regulatory friction) with `ALL` + `label`.
- `Vertical { name, display, rank, scores: [CriterionScore; 4], score_use,
  exercises, rationale }`; `score(criterion)` and
  `weakest_criterion_score()` expose historical values for replay.
  `decision_score` returns `None` because every retained plan score is
  `ScoreUse::SupersededForDecisionUse`.
- `verticals()` — the three historically ranked verticals; `chosen_wedge()` —
  the plan's retained rank-1 proposal (conjugate heat transfer), not current
  decision authority; `four_criteria()`.
- `InputAxis` — kernel readiness, validation-data access, CAD burden, and
  compute cost.
- `Readiness` — `Present`, `Partial`, or `Absent`, with score ceilings 10, 7,
  and 2. These scores mean input readiness, not physics accuracy or commercial
  attractiveness.
- `Measurement { readiness, score, method, evidence, finding }` plus
  `EvidencePointer { kind, reference, locator }`. Methods distinguish direct
  workspace inventory, contract-boundary review, official-dataset review, and
  static complexity analysis. Evidence kinds distinguish tracked workspace
  paths, Beads, and official publisher URLs.
- `KernelReadinessEntry`, `ValidationDataEntry`, `CadBurdenEntry`, and
  `ComputeCostEntry` carry their domain-specific fields around a common
  measurement. Compute envelopes state variables and operation/complexity
  shape; they are not wall-time estimates.
- `MeasuredWedgeInputs { vertical, measured_on, kernels, validation_data,
  cad_burden, compute_cost }`; `measured_wedge_inputs()` returns one record per
  candidate and `measured_inputs_for` resolves by vertical slug.
- `CycleTimeBaseline` + `CHT_BASELINE` — `baseline_days`, `target_reduction`
  (`3.0`), `kill_within_quarters`; `meets_kill_criterion(measured_days)`.
- `audit() -> WedgeAudit` (+ `STRONG_THRESHOLD`); `to_json()`.

## Invariants

- Historical plan scores never authorize a current decision. They remain
  byte-stable inputs for replay only.
- Every candidate has at least one evidence-complete record on all four
  measured axes. Every measurement has a non-empty method, finding, and
  evidence pointer.
- `Measurement.score <= Readiness::score_ceiling()`. In particular, an absent
  capability can never reach `STRONG_THRESHOLD` (`8`).
- Exactly three verticals ranked 1, 2, 3; each names at least one exercised
  proposal (V1→2/1/3/12, V2→1, V3→11/4).
- The kill criterion is measurable: `target_reduction == 3.0`, and
  `meets_kill_criterion` guards divide-by-zero.

## Error model

Total functions; no panics.

## Determinism class

Fully deterministic: pure `const` data; inventory and evidence order are
fixed; `to_json` emits the complete historical and measured ledger
byte-for-byte.

## Cancellation behavior

None.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/wedge.rs` (Proposal 7): the beachhead identity; historical-score
supersession; complete measured inputs on all four axes; status/score ceilings;
three ranked verticals with proposal mappings; measurable cycle-time shape
(including divide-by-zero guard); complete audit; negative doctrine and unique
labels; deterministic JSON. The workspace evidence test reads every
`WorkspacePath`, checks its locator marker, prints a deterministic
`PASS`/`FAIL` table, and fails on drift.

## No-claim boundaries

- The four historical criterion scores are the plan's strategic judgment, not
  empirical measurements, and are explicitly superseded for decision use.
- The measured-input scores classify readiness of the stated evidence at the
  dated inventory snapshot. They do not aggregate into a replacement wedge
  rank, prove model accuracy, or predict adoption/ROI.
- The CHT `correlation-Nu`, `RANS`, and `LES` entries in `fs-ladder` are rung
  declarations with a generic `Refine1d` demonstrator; the ladder contract says
  it does not run solves. No correlation catalog, fan curve, RANS/LES solver,
  or solid-fluid thermal transfer is inferred from those labels.
- `fs-lbm::ThermalLbm` is measured present only for its implemented
  two-dimensional Boussinesq slab. It is not promoted into an electronics CHT
  kernel.
- `fs-adjoint::HeatAdjoint` owns a backward-Euler reference problem over
  caller-assembled matrices. It is not a CHT assembler or coupled adjoint.
- `fs-vpm` is a two-dimensional inviscid direct kernel and `fs-couple`'s FSI
  fixture is a scalar linearized map. Neither is a coupled flutter solver.
- AM Bench data access is recorded from NIST's official data-management pages;
  a specific case/version/file/checksum and dataset-specific reuse terms remain
  to be pinned. The NASA/AGARD and Sandia records similarly remain partial
  where raw packaging or explicit reuse terms are not pinned.
- CAD burden compares each vertical with `fs-io`'s strict faceted STEP subset.
  It does not treat external tessellation as native assembly, units, material,
  NURBS, shell, or process semantics.
- Static compute envelopes describe loop/operation scaling only. They make no
  wall-time, memory-residency, accuracy, convergence, or performance claim.
- The cycle-time baseline (`5` days) is still a scoping placeholder to be
  replaced with a real measured customer baseline per the acceptance criteria.
- The kill criterion (`>= 3×` within two quarters of GA) is a COMMERCIAL gate
  on the wedge, not the architecture — a miss means re-select the vertical, not
  change the platform.
