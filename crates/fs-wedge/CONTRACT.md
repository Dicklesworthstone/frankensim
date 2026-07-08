# CONTRACT: fs-wedge

Go-to-market wedge selection as data (plan addendum, Proposal 7): the
conjugate-heat-transfer beachhead, scored, with its second/third verticals and
the measurable cycle-time kill criterion.

## Purpose and layer

Layer UTIL (pure data + audit; no dependencies). This gates NOTHING technical —
it constrains vertical-specific kernel work and records the commercial bet.

## Public types and semantics

- `WEDGE_DOCTRINE` — the load-bearing NEGATIVE rule ("do not sell against peak
  single-physics fidelity").
- `WedgeCriterion` (4: kernel maturity, iteration pain, quantifiable ROI, low
  regulatory friction) with `ALL` + `label`.
- `Vertical { name, display, rank, scores: [CriterionScore; 4], exercises,
  rationale }`; `score(criterion)` and `weakest_criterion_score()`.
- `verticals()` — the three ranked verticals; `chosen_wedge()` — the rank-1
  beachhead (conjugate heat transfer); `four_criteria()`.
- `CycleTimeBaseline` + `CHT_BASELINE` — `baseline_days`, `target_reduction`
  (`3.0`), `kill_within_quarters`; `meets_kill_criterion(measured_days)`.
- `audit() -> WedgeAudit` (+ `STRONG_THRESHOLD`); `to_json()`.

## Invariants

- The chosen wedge is STRONG (`>= STRONG_THRESHOLD`) on EVERY criterion — a
  wedge weak on any of the four is not a wedge.
- Exactly three verticals ranked 1, 2, 3; each names at least one exercised
  proposal (V1→2/1/3/12, V2→1, V3→11/4).
- The kill criterion is measurable: `target_reduction == 3.0`, and
  `meets_kill_criterion` guards divide-by-zero.

## Error model

Total functions; no panics.

## Determinism class

Fully deterministic: pure `const` data; `to_json` reproduces byte-for-byte.

## Cancellation behavior

None.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/wedge.rs` (Proposal 7, 7 cases): the beachhead is conjugate heat
transfer; the chosen wedge is strong on all four criteria; three ranked
verticals with proposal mappings; the measurable cycle-time kill criterion
(incl. divide-by-zero guard); the complete audit; the negative doctrine +
unique labels; deterministic JSON.

## No-claim boundaries

- The four-criteria SCORES are the plan's strategic judgement encoded as data,
  not empirical measurements; the cycle-time baseline (`5` days) is a scoping
  placeholder to be replaced with a real measured baseline per the acceptance
  criteria.
- The CHT correlation-based bottom rung that makes this vertical concrete is
  implemented in `fs-ladder` (`LadderRegistry::cht()`); this crate records the
  strategic selection, not the kernels.
- The kill criterion (`>= 3×` within two quarters of GA) is a COMMERCIAL gate
  on the wedge, not the architecture — a miss means re-select the vertical, not
  change the platform.
