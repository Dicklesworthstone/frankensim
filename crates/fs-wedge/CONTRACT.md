# CONTRACT: fs-wedge

Go-to-market wedge selection as data (plan addendum, Proposal 7): the
historical three-vertical ranking, its superseded judgment scores, measured
decision inputs, the explicit three-candidate recommendation with one-factor
sensitivity tables, and the cycle-time kill criterion.

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
  static complexity analysis. `DecisionAssumptionReview` explicitly marks
  commercial or schedule judgments whose empirical measurement is pending.
  Evidence kinds distinguish tracked workspace paths, Beads, and official
  publisher URLs.
- `KernelReadinessEntry`, `ValidationDataEntry`, `CadBurdenEntry`, and
  `ComputeCostEntry` carry their domain-specific fields around a common
  measurement. Compute envelopes state variables and operation/complexity
  shape; they are not wall-time estimates.
- `MeasuredWedgeInputs { vertical, measured_on, kernels, validation_data,
  cad_burden, compute_cost }`; `measured_wedge_inputs()` returns one record per
  candidate and `measured_inputs_for` resolves by vertical slug.
- `ScoringFactor` defines nine higher-is-better factors: customer pain, kernel
  readiness, validation tractability, data access, low CAD burden, short time
  to decision, differentiation, low compute cost, and low regulatory risk.
  `DEFAULT_FACTOR_WEIGHTS` records integer percentage weights summing to 100.
- `FactorRating` keeps the comparative `0..=10` rating separate from the
  attached `Measurement`'s evidence authority. `ComparisonCandidate` carries
  exactly one such input per factor, an inventory revision, and a minority
  case. `comparison_candidates()` returns full electronics CHT, SDF
  structural/topology assurance, and thermal design assurance.
- `score_candidates` is the pure scoring function. It validates weights and
  candidates, computes the integer sum `rating * weight`, and orders by total
  descending then stable candidate slug ascending. `tilted_weights` raises one
  weight and proportionally rescales the others to 100 using deterministic
  largest-remainder allocation.
- `ranked_recommendation` and `default_recommendation` return the ranking,
  runner-up minority report, and exhaustive single-factor flip tables. Rating
  sensitivity raises only one challenger rating. Weight sensitivity raises one
  factor weight while proportionally reducing every other weight.
- `render_comparison_report()` emits the complete weights, factor evidence,
  rank, minority report, and both sensitivity tables as a deterministic TSV-like
  artifact.
- `CycleTimeBaseline` + `CHT_BASELINE` — a per-step sourced envelope over the
  six `IncumbentStep`s (CAD preparation, meshing, solver setup, solve,
  post-processing, report assembly), never a point value. Each
  `BaselineStepEstimate` carries `low_hours`/`high_hours`, a bound-derivation
  basis, and `SourcedFigure` rows (figure, verbatim quote, citation, URL,
  `DurationSourceClass`, bias caveat). `BaselineProvenance` distinguishes
  `Placeholder` (refused), `PublishedSourceDerived` (current), and
  `ExecutedRun` (supersedes on arrival). The dossier
  `data/cycle-time-baseline-dossier.md` (marker
  `fs-wedge-cycle-time-baseline-dossier-v1`) is the replay handle for every
  figure.
- `evaluate_kill_criterion(measured_days)` returns a
  `KillCriterionEvaluation` with a conservative three-way `KillVerdict`:
  `Met` only when the LOW baseline bound clears the `3.0x` target, `NotMet`
  only when even the HIGH bound cannot, `Indeterminate` otherwise. The
  evaluation logs which record it used and the full derivation.
  `RETIRED_PLACEHOLDER_BASELINE` retains the old flat 5-day figure solely to
  prove the refusal path.
- `audit() -> WedgeAudit` (+ `STRONG_THRESHOLD`); `to_json()`.

## Invariants

- Historical plan scores never authorize a current decision. They remain
  byte-stable inputs for replay only.
- Every candidate has at least one evidence-complete record on all four
  measured axes. Every measurement has a non-empty method, finding, and
  evidence pointer.
- Every explicit comparison candidate has exactly one complete input for each
  scoring factor. Comparative ratings are in `0..=10`; evidence-authority
  scores independently obey the readiness ceiling.
- Default weights contain each factor exactly once and sum to 100. Weighted
  totals therefore lie in `0..=1000`. Exact total ties resolve by ascending
  candidate slug, independent of candidate or weight input order.
- At inventory revision `e5c8061f4faed986b831b8978d0c8d1812e960fb`, the
  recorded default totals are thermal design assurance `638`, SDF structural
  assurance `623`, and full electronics CHT `502`. Thermal design assurance is
  the provisional recommendation; SDF structural assurance is the runner-up
  whose minority report is retained.
- `Measurement.score <= Readiness::score_ceiling()`. In particular, an absent
  capability can never reach `STRONG_THRESHOLD` (`8`).
- Exactly three verticals ranked 1, 2, 3; each names at least one exercised
  proposal (V1→2/1/3/12, V2→1, V3→11/4).
- The kill criterion is measurable: `target_reduction == 3.0`; the baseline
  record in the decision path is complete and non-placeholder (one estimate
  per step in pipeline order, every source complete, at least one
  non-vendor source per step); and the placeholder refusal drill passes in
  `audit()` (`cycle-time-baseline-measured`, `placeholder-baseline-refused`).

## Error model

Static accessors are total. Scoring returns `ScoringError` for non-normalized,
missing, or duplicate weights; missing, duplicate, or incomplete candidates;
invalid factor inputs; and non-increasing weight tilts. The kill-criterion
evaluation returns `KillCriterionError` for placeholder baselines, incomplete
baseline records, and non-positive or non-finite measured cycle times; it never
coerces or defaults. Internal `expect` calls are applied only after those
structural validations and to infallible writes to `String`.

## Determinism class

Fully deterministic: pure `const` inputs; inventory and evidence order are
fixed; integer scoring avoids floating-point rank drift; ties use stable slugs;
largest-remainder weight normalization breaks equal remainders by canonical
factor order. `to_json` and `render_comparison_report` replay byte-for-byte.

## Cancellation behavior

None.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/wedge.rs` (Proposal 7): the beachhead identity; historical-score
supersession; complete measured inputs on all four axes; status/score ceilings;
three ranked verticals with proposal mappings; explicit factor completeness and
recorded totals; weight normalization/refusals; monotonicity of every factor;
candidate-permutation and exact-tie determinism; exhaustive rating/weight flip
sensitivity including degenerate full-weight ties; deterministic verbose
report; the conservative kill-verdict trichotomy with refusal of
non-measurable inputs; measured-baseline completeness (pipeline step order,
load-bearing sources, https URLs); dossier existence, marker, and URL
cross-citation; envelope arithmetic against the dossier derivation; the typed
placeholder refusal; the printed kill-criterion derivation; manifest numeric
round-trip of the envelope; the audit's baseline and refusal-drill checks;
complete audit; negative doctrine and unique labels; deterministic JSON. The
workspace evidence test reads every `WorkspacePath`, checks its locator
marker, prints a deterministic `PASS`/`FAIL` table, and fails on drift.

## No-claim boundaries

- The four historical criterion scores are the plan's strategic judgment, not
  empirical measurements, and are explicitly superseded for decision use.
- The measured-input scores classify readiness of the stated evidence at the
  dated inventory snapshot. They do not aggregate into a replacement wedge
  rank, prove model accuracy, or predict adoption/ROI. The later explicit
  comparison uses separate comparative ratings rather than repurposing these
  readiness scores.
- The explicit ranking is a reproducible decision model, not empirical proof
  that the recommended market exists or that any candidate will achieve a
  delivery schedule, cycle-time reduction, adoption, accuracy, or return on
  investment. Customer-pain, time-to-decision, and regulatory inputs remain
  declared assumptions wherever no retained measurement exists.
- Evidence authority and comparative desirability are separate axes. A high
  factor rating cannot promote an `Absent` or `Partial` measurement into a
  stronger scientific claim, and the weighted total is not a certificate.
- Weight sensitivity permits a factor to reach 100, which can zero every other
  factor. In that degenerate case a recommendation may flip solely through the
  documented slug tie-break; the report exposes rather than hides that result.
- The recommendation is provisional pending ratification and successor
  customer-baseline work. It records the repository inventory at the pinned
  revision above and does not silently incorporate later workspace changes.
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
- The cycle-time baseline is `PublishedSourceDerived` and therefore an
  ESTIMATED quantity: nobody on this project timed the incumbent run. Its
  authority ends at the cited sources (Sandia DART 2005, the vendor-run
  6SigmaET 2018 survey, Tech-Clarity 2016, peer-reviewed solve timings), whose
  bias classes and caveats travel with each figure. The envelope covers one
  FIRST-PASS iteration of the representative task; re-iterations are cheaper
  and were deliberately not used to narrow it.
- The baseline is NOT customer-pain evidence: it is the kill-criterion
  denominator only, and it must not upgrade the comparison's customer-pain
  factor, which stays a declared assumption until interviews or workflow
  traces exist.
- The report-assembly row is the envelope's weakest (no independent
  electronics-specific figure exists in public sources) and is flagged for the
  quarterly review. An executed incumbent-workflow run supersedes every survey
  row under the dossier protocol.
- An `Indeterminate` kill verdict means the envelope straddles the target: the
  reduction claim may not be stated as met, and only a measured FrankenSim
  cycle time at or below `baseline_days_low / 3` supports "met" against even
  the fastest incumbents.
- The kill criterion (`>= 3×` within two quarters of GA) is a COMMERCIAL gate
  on the wedge, not the architecture — a miss means re-select the vertical, not
  change the platform.
