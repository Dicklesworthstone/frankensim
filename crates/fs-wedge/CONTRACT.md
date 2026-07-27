# CONTRACT: fs-wedge

Go-to-market wedge selection as data (plan addendum, Proposal 7): the
historical three-vertical ranking, its superseded judgment scores, measured
decision inputs, the explicit three-candidate recommendation with one-factor
sensitivity tables, and the cycle-time kill criterion.

## Purpose and layer

Layer UTIL (pure decision data, bounded offline historical replay, and audit).
Its only production crate dependency is `fs-blake3`, used for domain-separated
artifact and source identities. This crate gates no scientific or
kernel-capability claim; it does fail closed when the default recommendation's
retained evidence cannot be authenticated and replayed.

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
- `historical` defines the retained-comparison replay protocol.
  `HistoricalEvidenceSnapshot` binds schema and policy versions, the reviewed
  inventory revision, manifest and bundle paths, domain-separated BLAKE3
  identities, and admitted source/pointer counts.
  `COMPARISON_EVIDENCE_SNAPSHOT` is the descriptor for the embedded
  `b3b5f2c1...` snapshot.
- `comparison_evidence_manifest()` and `comparison_evidence_bundle()` expose
  the exact embedded TSV manifest and TAR bytes.
  `verify_default_comparison_evidence()` replays those artifacts without
  reading Git or the filesystem. `verify_comparison_evidence()` admits
  caller-supplied artifacts for mutation and migration tests only after the
  supplied descriptor authenticates them.
- Successful replay returns `HistoricalEvidenceReceipt`, exposing schema,
  policy, revision, both recomputed artifact identities, and authenticated
  source/pointer counts. `HistoricalEvidenceError` carries the fail-closed
  diagnostic, including per-source expected/observed BLAKE3 identities when
  exact content differs.
- `score_candidates` is the pure scoring function. It validates weights and
  candidates, computes the integer sum `rating * weight`, and orders by total
  descending then stable candidate slug ascending. `tilted_weights` raises one
  weight and proportionally rescales the others to 100 using deterministic
  largest-remainder allocation.
- `ranked_recommendation` is the pure caller-supplied scoring/sensitivity
  routine. `default_recommendation` first requires successful replay of the
  retained historical evidence, then returns the default ranking, runner-up
  minority report, and exhaustive single-factor flip tables. Rating
  sensitivity raises only one challenger rating. Weight sensitivity raises one
  factor weight while proportionally reducing every other weight.
- `render_comparison_report()` emits `FS-WEDGE-COMPARISON v2`: one
  authenticated `HISTORICAL_EVIDENCE` receipt row followed by the complete
  weights, factor evidence, ranks, minority report, and both sensitivity tables
  as a deterministic TSV-like artifact.
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
- At inventory revision `b3b5f2c1c809eec06cde1e40cbc916d6995469b5`, the
  recorded default totals are thermal design assurance `638`, SDF structural
  assurance `623`, and full electronics CHT `502`. Thermal design assurance is
  the provisional recommendation; SDF structural assurance is the runner-up
  whose minority report is retained.
- The default comparison is bound to that revision through an embedded
  canonical manifest and retained source TAR. Replay authenticates the exact
  manifest and TAR bytes with separate domain-separated BLAKE3 identities
  before parsing either artifact.
- The canonical v1 manifest contains one revision row, 13 strictly path-sorted
  `SOURCE` rows, and the exact canonical sequence of 31 historical
  `WorkspacePath` pointer occurrences. Each source row carries byte length,
  extraction-audit SHA-256 and Git blob metadata, and a domain-separated BLAKE3
  over revision, path, and exact source bytes. Every candidate revision must
  equal the snapshot revision; every retained source must be consumed;
  missing, extra, duplicate, reordered, or identity-mismatched evidence fails
  closed.
- The retained TAR must have valid header checksums, zero payload padding, at
  least two zero termination blocks, exactly one leading matching revision PAX
  record, normalized paths, supported entry types, and exactly the
  manifest-declared regular-file set and byte lengths. Every historical locator
  must occur in its identity-verified retained source bytes.
- Historical comparison replay reads no live workspace or Git state. Current
  measured-input drift checks remain intentionally live and are a separate
  authority lane.
- `Measurement.score <= Readiness::score_ceiling()`. In particular, an absent
  capability can never reach `STRONG_THRESHOLD` (`8`).
- Exactly three verticals ranked 1, 2, 3; each names at least one exercised
  proposal (V1→2/1/3/12, V2→1, V3→11/4).
- The kill criterion is measurable: `target_reduction == 3.0`; the baseline
  record in the decision path is complete and non-placeholder (one estimate
  per step in pipeline order, every source complete, at least one
  non-vendor source per step); and the placeholder refusal drill passes in
  `audit()` (`cycle-time-baseline-measured`, `placeholder-baseline-refused`).
- `audit()` requires `comparison-history-bound` in addition to comparison
  completeness, normalized weights, ranking, and sensitivity checks.

## Error model

Static data accessors are total. Pure scoring returns `ScoringError` for
non-normalized, missing, or duplicate weights; missing, duplicate, or
incomplete candidates; invalid factor inputs; and non-increasing weight tilts.
`default_recommendation` and `render_comparison_report` additionally return
`ScoringError::HistoricalEvidenceUnavailable { source }`, preserving the exact
historical replay refusal.

Historical replay returns `HistoricalEvidenceError` for invalid descriptors,
artifacts beyond their admitted caps, whole-artifact or per-source identity
mismatches, malformed manifests or TARs, oversized caller evidence,
candidate-revision or pointer mismatches, unreferenced sources,
missing/extra/wrong-length source files, and missing locator markers. Admission
caps are 128 KiB per manifest, 1 MiB per bundle, 512 source or pointer rows, and
4096 UTF-8 bytes per manifest scalar.

The kill-criterion evaluation returns `KillCriterionError` for placeholder
baselines, incomplete baseline records, and non-positive or non-finite measured
cycle times; neither subsystem coerces or defaults invalid evidence. Internal
`expect` calls are applied only after structural validation, in tests, and to
infallible writes to `String`.

## Determinism class

Fully deterministic for identical compiled artifacts and inputs: embedded
evidence bytes and pointer order are fixed; manifest, TAR, and per-source
identities are domain separated; parsing uses canonical row/file ordering;
integer scoring avoids floating-point rank drift; ties use stable slugs; and
largest-remainder weight normalization breaks equal remainders by canonical
factor order. Historical replay receipts, `to_json`, and
`render_comparison_report` reproduce byte-for-byte.

## Cancellation behavior

None. Replay is synchronous and has no cancellation points; its parser work is
bounded by the documented artifact, row, scalar, and caller-evidence caps.

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
complete audit; negative doctrine and unique labels; deterministic JSON.

The live workspace-evidence test reads every `WorkspacePath` in
`measured_wedge_inputs()`, checks its current locator marker, prints a
deterministic `PASS`/`FAIL` table, and fails on drift. Historical comparison
pointers are deliberately excluded from that live scan and instead replay from
the retained snapshot.

Historical replay tests prove deterministic 13-source/31-pointer receipts and
independence from today's workspace. G3 mutations refuse unauthenticated
manifest/TAR bytes, candidate-revision drift, a renamed/missing TAR source, an
empty snapshot, an authenticated pointer-path change, a same-length non-marker
content change after whole-bundle reauthentication, and an
identity-reauthenticated source whose locator was removed. Unit tests
independently recompute all 13 per-source identities and both top-level
artifact identities. Report tests require the v2 header and exactly one
`HISTORICAL_EVIDENCE` row; audit tests require `comparison-history-bound` and
13 total checks.

A separate exhaustive kernel-matrix test independently probes the required
implementation markers for all 15 kernel rows, derives `present`, `partial`,
or `absent` from the number found, prints a labeled matrix/observed diff, and
fails when a new or removed module makes the recorded readiness stale. The
cross-crate `fs-govern/tests/wedge_audit.rs` e2e lane composes these measured
inputs, scores, exhaustive sensitivities, sourced baseline, and kill
derivation with the fail-closed ratification record. Its `wedge-audit` binary
emits one content-addressed artifact plus reconstructable JSON-lines logs and
proves a seeded missing-cycle-time-evidence pointer fails without partial
output.

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
- Artifact authentication proves that replay consumed the descriptor-bound
  manifest, TAR, and revision/path/source identities and that their recorded
  paths, lengths, pointer sequence, revision label, and locator substrings
  satisfy the v1 protocol. It does not prove that the retained statements are
  scientifically true, complete, current, or correctly interpreted, and it
  does not validate the comparative ratings or weights.
- The Git revision, Git blob OIDs, extraction-audit SHA-256 values, and TAR PAX
  comment are retained provenance metadata, not an independent Git trust
  proof. Offline replay does not contact a repository, recompute Git/SHA-256
  identities, validate commit ancestry or signatures, or establish that the
  revision was reviewed by an authorized human. Exact runtime content binding
  comes from the separately domain-bound manifest, TAR, and per-source BLAKE3
  identities.
- Locator replay proves substring presence only. It does not certify the
  surrounding contract's semantics, implementation behavior, model validity,
  or evidence authority.
- `verify_comparison_evidence` authenticates bytes against its
  caller-supplied descriptor, so its receipt proves protocol self-consistency,
  not canonical review authority. Only the embedded constant plus
  `verify_default_comparison_evidence` names this crate's reviewed default
  trust root.
- `ranked_recommendation` remains a pure scoring primitive for caller-supplied
  records. Only the default recommendation/report/audit path requires the
  retained-history receipt, and that receipt does not turn a weighted score
  into a scientific or commercial certificate.
- Weight sensitivity permits a factor to reach 100, which can zero every other
  factor. In that degenerate case a recommendation may flip solely through the
  documented slug tie-break; the report exposes rather than hides that result.
- The recommendation is provisional pending ratification and successor
  customer-baseline work. Its comparison factors replay exclusively from the
  retained `b3b5f2c1...` source bundle and therefore do not silently
  incorporate later workspace changes. Separately dated measured-readiness
  inputs continue to inspect the live workspace and must not be confused with
  the historical comparison snapshot.
- `fs-convection` and `fs-airflow` now provide the typed correlation catalog,
  fan curve, flow-network operating point, and evidence-preserving correlation
  handoff used by the low-cost CHT rung. Their synthetic fixtures and
  conditional mathematical certificates do not establish manufacturer or
  enclosure accuracy. The `RANS` and `LES` entries in `fs-ladder` remain rung
  declarations, not solvers, and no solid-fluid thermal field transfer is
  inferred from them.
- `fs_airflow::conjugate` (bead f85xj.5.7) adds the partitioned conjugate
  exchange the `conjugate-heat-transfer` vertical's `kernel-maturity` rationale
  now cites: a stream-wise 1-D air path and a solid/air fixed point over the
  per-region Robin reference temperature. That rationale previously claimed
  "forced-convection CFD" and credited `fs-ladder`'s `cht()` bottom rung
  directly; both were overstatements and are corrected in place with the prior
  text retained. The correction does not raise the vertical: the exchange is a
  correlation-rung coupling with a frozen `h`, its fixtures are synthetic, and
  no CFD kernel, RANS rung, or field transfer exists. `ScoreUse::
  SupersededForDecisionUse` continues to deny every retained score decision
  authority.
- `fs-lbm::ThermalLbm` is measured present only for its implemented
  two-dimensional Boussinesq slab. It is not promoted into an electronics CHT
  kernel.
- `fs-adjoint::HeatAdjoint` owns a backward-Euler reference problem over
  caller-assembled matrices. It is not a CHT assembler or coupled adjoint.
- `fs-vpm` is a two-dimensional inviscid direct kernel and `fs-couple`'s FSI
  fixture is a scalar linearized map. `fs-flutter-e2e` adds a deterministic
  two-degree-of-freedom stability-boundary campaign, not a real aeroelastic
  objective or verified coupled gradient.
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
