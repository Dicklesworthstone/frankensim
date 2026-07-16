# CONTRACT: fs-plan

> Status: ACTIVE (models v1). Owns per-operator cost/error models and the
> Error/Time attribution ledgers. The budget ALLOCATOR that optimizes over
> these models is the gp3.9 bead's; DWR/conformal/MLMC estimator inputs
> arrive with their owning crates.

## Purpose and layer

Plan §11.4 (Bet 12), Decalogue P4: every operator publishes an error model
and a cost model; composition composes the models, so "how accurate is
this number and where did the error come from" — and "where did the
seconds go" — are queries over attribution trees. Layer: L6 (HELM).
Runtime deps: `std`, fs-blake3, fs-geom, fs-ledger; feature-gated VoI uses
asupersync for cancellation-aware decision-oracle evaluation and fs-eproc for
the anytime-valid audit threshold used by authenticated authority.

## Public types and semantics

- `CostModel` — log-log power-law fit (`cost ≈ exp(a)·size^b`) with
  EMPIRICAL residual quantiles supplying P10/P50/P90 bands; predictions
  carry observation counts and an extrapolation flag (an estimate is
  itself evidenced). Below `MIN_OBS` it refuses (`CostRefusal`) instead of
  guessing. Fits use centered, finite-checked regression; batch fitting sorts
  once and online insertion preserves total order. Observation and held-out
  work are explicitly bounded. `observe()` is transactional; invalid input,
  an unstable/nonfinite result, or capacity exhaustion leaves the prior model
  unchanged. `calibration()` audits held-out band coverage;
  `median_rel_error()` is the improvement metric. Empty probe sets refuse
  rather than returning vacuous perfect scores.
- `ErrorLedger` — attribution tree over the plan's canonical sources
  (geometry, discretization, algebraic, surrogate, statistical,
  model-form), each entry carrying a non-blank operator identity and its
  `Rigor` class (certified /
  estimated / rate-model). First-order ADDITIVE composition; `total()`
  = Σ contributions + declared residual; `lint()` refuses NaN/negative
  mass and aggregate overflow (no silent error, ever); `dominant()` names
  what escalation should attack; `explain()` renders strict JSON and emits an
  explicit `valid:false` diagnostic instead of non-JSON NaN/Infinity.
- `TimeLedger` — per-stage predicted quantile bands vs measured wall
  seconds, with a band-coverage calibration audit, a lint for finite
  nonnegative ordered quantiles/measured times and aggregate overflow, and a
  strict-JSON `explain()` (`null`, never Rust `Some`/`None` syntax).
- `PlanCostOracle` — implements fs-geom's `CostOracle`, so the Rep Router
  plans with THIS machine's measured history (models are bound to the complete
  `ConverterSpec` at registered reference sizes; sealed execution actuals feed
  refits; the maximum observed error upper bound backs
  `measured_error_abs`). Raw scalar recording is private. Reads and writes are
  fallible and bounded; a record is atomic across cost and error histories.
  Unregistered or same-name/different-spec edges, nonfinite values, capacity
  exhaustion, and reference-size changes after observation refuse without
  changing routing authority. Router never uses observed errors to tighten a
  declared hard bound.
- `sealed::SealedCostModel` / `SealedCostPrediction` / `CostModelScope` /
  `CostEvidenceClass` (bead 2pmb) — authority carriers for cost models.
- `FreshnessPolicy` / `StalenessVerdict` / `SealedCostModel::assess` /
  `predict_assessed` (bead jle3m) — the class DEGRADES with age and
  machine drift under a caller-preregistered horizon: receipt-backed
  evidence assessed past the horizon, under a different machine
  fingerprint, or with a future timestamp becomes
  `StaleRooflineReceipt`, which ranks BELOW a fresh provisional fit in
  the total lattice (`rank`/`weakest`: fresh exact > fresh provisional
  > stale exact) — a year-old sealed roofline can no longer outrank a
  fresh fit forever. Assessment is pure (caller supplies `now_ns` and
  the current machine fingerprint; this crate never reads a clock or
  host identity); the mint-time stamp is immutable and the numeric
  bands never change — only the evidence label degrades.
  `CostModel` stays freely constructible provisional MATH; the sealed
  type is the AUTHORITY: its `ExactRooflineReceipt` class is mintable
  only by `cost_model_from_tune` after full validation, with the
  validated scope (kernel, shape class, machine key, run receipt, op id,
  build identity, recorded-at) retained behind private fields so a
  caller-fitted model can never impersonate receipt-backed evidence.
  `provisional_unaudited` wraps caller fits as visibly provisional; every
  prediction carries the scope fields and mint-time class, and
  composition never upgrades the class.
- `cost_model_from_tune` — rebuilds a SEALED model from one EXACT
  `(kernel, shape_class, machine)` fs-ledger key. It accepts only the current
  fs-roofline receipt-v3 / row-v4 / `roofline-v8` production-v3 schema
  through a bounded strict JSON parser that rejects duplicate and unknown
  keys. Every timed sample is
  retained as a same-size observation; median/P25/P75/dispersion and rate are
  rederived from those samples. Embedded kernel/version/run/op/repetition
  identities, the baseline-bound 40-byte machine key, completed production-v3
  op envelope, build identity, result-manifest membership, payload artifact
  bytes/metadata/OUT edge, and dependency-receipt bytes/domain digest/IN edge
  must all agree before the model sees any evidence. `baseline_admission` is
  not an opaque object: the consumer requires the byte-canonical
  `fs-roofline-axis-admission-v2` attested envelope, independently rederives
  its baseline-record hash and trusted axis/age math, requires exact
  baseline/identity/source-list agreement plus authorized/trusted verdicts,
  binds pre/post axis receipt hashes and row post-axis bits, and requires the
  frozen decision day to equal the operation's UTC completion day. Historical
  replay is host-independent: OS and architecture must agree inside the bound
  baseline/admission identity, but are never compared with the audit process's
  ambient target. The exact caller-supplied 40-byte machine key still prevents
  cross-machine model lookup. It then reconstructs finalized-run v3 over the
  exact admission bytes, every ordered manifest sibling's currently stored
  measured bytes, and the manifest child digest. Production-v3 accepts only the
  sealed ordered `simd-axpy-f64/1`, `simd-dot-f64/1`, `simd-sum-f64/1`, and
  `gemm-f64/2` registry. All siblings must bind one repetition/warmup profile;
  vector element counts must equal the sealed `n`, GEMM elements must equal
  `max(isqrt(n), 256)^2`, and the consumer independently enforces `n <= 2^24`,
  at most 64 total warmup-plus-timed runs, the 4096-worker ceiling, and the
  producer's checked `2^39`-FLOP / `2^33`-logical-byte aggregate caps. A
  missing/rekeyed/tampered sibling, mismatched sibling configuration, impossible
  producer workload, policy receipt, source list, baseline, or manifest refuses
  the complete run. Source-pin tests make duplicated L6 consumer constants and
  sealed-registry algebra fail when fs-roofline rotates its schema.
  fs-plan does not re-run the external promotion authority or prove retained
  source bytes; it consumes the exact authorized policy decision sealed into
  the producer's recomputed finalized-run identity. The op row uses
  fs-ledger's metadata-preflighted bounded read; the measured artifact is
  capped at its exact tune-row byte length and the fs-la dependency receipt at
  its producer-owned 1 MiB ceiling before either payload is materialized; a
  source-pin unit test fails if fs-la changes that declaration independently.
  Foreign machine/shape
  rows are never scanned or mixed.

- `voi` module (addendum Proposal C, bead knh1.6; [F], behind
  `voi-queries`): the ignorance market, v0 as a RANKED LIST.
  `UncertaintyNode` (interval + nominal), `DecisionBudget` (explicit
  evaluation and declared-work caps), and `DecisionOracle` (fallible
  evaluation under an asupersync `Cx` and a library-issued private permit).
  The permit binds canonical ordinal, exact query size, charged work,
  remaining envelope, and the caller's full budget; callers cannot forge one.
  `LiveDecision` is the cached synchronous adapter: it charges exactly one
  declared work unit and brackets the callback with library-owned checkpoints,
  but makes no time or memory claim for an arbitrary closure.
  `flip_probability` sweeps one node's interval at grid cost,
  `Probe`/`ProbeKind` (the priced menu UNIFYING computational and
  physical evidence), `rank_purchases` (MYOPIC one-step sampled VoI:
  flip-fraction reduction per dollar, deterministic tie-breaks) returning a
  sealed complete `RankedMenu`. Its private rows cannot be omitted, spliced,
  or reordered. Its `DecisionComputationReceipt` retains the exact evaluation
  count, exact declared-work charge, and caller-supplied budget.
  `source_context_id` binds the caller-declared policy and decision snapshot,
  validated nodes, complete supplied source menu, grid, oracle metadata, exact
  declared-work charge, and full budget;
  final `context_id` additionally binds every canonical output row and score.
  These roots are registered as `fs-plan:voi-ranked-source` and
  `fs-plan:voi-ranked-menu`, both at identity v2 under their existing BLAKE3
  domains. Artifact domain and producer version are explicit semantic inputs
  for all three identities. The ranked-source grammar binds explicit counts,
  node input order, and exact `f64::to_bits()` payloads. Caller probe input
  order is deliberately nonsemantic because validation requires unique names
  and the encoder sorts by name; that canonical sort is schema behavior, not a
  separate logical semantic field. The ranked-menu grammar binds the exact
  source root, row count/order, ranked probe names, and exact sampled/score
  bits. `fs-plan:voi-audit-context` v2 separately binds policy, fixed alpha,
  audit cap, chronological count/order, and every exact matched-audit field.
  The declarations pin each local domain/version constant, plus the audit
  alpha and record cap, and fingerprint the exact in-tree fs-blake3 constants
  and function bodies. The ranked-menu schema depends explicitly on the
  ranked-source schema whose root it embeds.
  Retained consumers must admit the producer-declared version through
  `RankedMenu::admit_retained_identity_versions` or
  `AuditReport::admit_retained_identity_version` before comparing a retained
  root; stale and future versions fail closed.
  `QueryHint` is structured and grid-qualified, with escaped lossless text and
  strict JSON renderers. A sampled zero is explicitly non-authoritative.
  `MatchedAuditRecord` validates private matched-cost observations with bounded
  identities/provenance and rebuilds caller-owned strings before retention so
  spare capacity cannot bypass the record cap. A non-cloneable `VoiScheduler` owns one append-only
  chronological fixed-alpha fs-eproc `PairwiseRace`, the remaining total
  budget, and a bounded set of consumed decision snapshots. `observe_audit`
  updates that one live process transactionally; `schedule(&mut self, menu)`
  rechecks the current verdict, refuses a foreign policy or repeated snapshot,
  and returns at most one highest-value affordable `ScheduledPurchase` receipt
  retaining source/final menu roots, policy/snapshot, audit root/e-value
  support, and the exact internal budget transition. `audit_scheduling` is a
  reporting-only replay helper and returns no spending capability.
  The raw `observe_audit` / `schedule` surface remains explicitly advisory.
  Ledger authority uses the separate `VoiLedgerContext`,
  `record_ranked_menu_receipt`, `record_prospective_audit_receipt`, independent
  `verify_*` functions, opaque `AuthenticatedRankedMenu` /
  `AuthenticatedAuditRecord` capabilities, and
  `schedule_authenticated`. Producer-facing artifact-kind/metadata constants,
  `voi_decision_artifact_metadata`, `voi_model_artifact_metadata`,
  `voi_probe_catalog_bytes`, and `voi_audit_outcome_bytes` define the exact
  admitted envelopes. A context binds session, decision/model/policy
  versions, VCS branch and snapshot root, issuance day, and inclusive expiry.
  Ranked receipts require decision, model, and canonical complete-catalog
  artifacts to be finished outputs in the current bound snapshot. Prospective
  audit receipts bind their zero-based order and exact recommended/alternative
  outcome artifacts, whose canonical bytes include candidate, matched cost,
  provenance, and realized decision-change bit. Receipt operations are
  deterministic finished ledger operations with exact Five Explicits, exact
  input/output lineage, and a sole sealed output producer. Verification
  rederives all bytes and envelopes rather than trusting a locator. Authenticated
  receipt minting is idempotent: retrying an exact canonical request reverifies
  and returns its existing sole-producer locator, while changed content has a
  different identity. The logical issuance day is checked for a lossless
  day-start nanosecond conversion before it is recorded in the ledger op
  envelope. Authenticated
  scheduling reverifies the ranked receipt against the current branch head and
  every audit receipt against its immutable known snapshot before any budget or
  consumed-snapshot mutation; all evidence must share one session,
  decision/model/policy version, branch, and exact decision/model artifact
  identities.
  Decision, sweep, ranking, and scheduling entry points are fallible:
  arity, nonempty bounded collections/names, unique identities, finite
  ordered intervals, nominal containment, callback margins, target
  resolution, grid and aggregate evaluation work, probe economics,
  ranked values, and budgets are validated before evaluation or spend.
  Probe-menu identity and score ordering are canonical. Matched-audit order is
  authoritative chronological input and is bound exactly; completed outcomes
  are never sorted through the adaptive e-process. Duplicate source/audit
  identities refuse, and a positive scheduled cost must strictly decrease the
  scheduler-owned finite budget.

- `alloc::{Knob, KnobSetting, AllocProblem, Plan, allocate, Allocator,
  AllocationError, PlanInputError, BudgetInfeasible, oracle_min_error}`
  (bead gp3.9 V1): the
  GREEDY-PLUS-LOOKUP budget allocator — Pareto-ladder upgrades by
  marginal utility Δerror/Δwall under a TROPICAL wall-clock (max over
  parallel tracks of within-track sums, §14.3: slack upgrades are free
  and taken first); online re-planning via `Allocator::observe_error`
  (a-posteriori estimates override model values transactionally and
  re-prune dominance). Constructors, observations, public evaluators,
  and the fixture oracle are fallible. Finite nonnegative scalars,
  Pareto ordering, choice vectors, tracks, knob/setting counts, aggregate
  arithmetic, and oracle Cartesian work are bounded before allocation or
  enumeration. An empty problem is the explicit zero-error/zero-wall
  identity. `allocate` returns only an in-budget plan, otherwise a typed
  malformed-input refusal, `MinimumPlanExceedsBudget`, or a structured
  `BudgetInfeasible` with ranked, VERIFIED relaxations. Its in-budget error is
  the deterministic V1 greedy stall, not an exact optimum or lower bound. For
  a reachable target, its budget suggestion is the wall clock of the same
  uncapped greedy prefix: re-planning at exactly that budget follows the same
  prefix and succeeds (gated), without claiming that the suggestion is the
  globally minimum feasible budget. Measured on the fixture matrix:
  worst greedy/oracle error ratio 1.134 over 60 random problems; the §11.4
  "drag to 2% in 2h" scenario plans to 1.68% at 3600 s with an 8-line
  rationale, deterministically.
- `moonshot::{optimize_exact, waterfill, cma_continuous, RateKnob,
  ScoreRow}` (bead gp3.9 V2, feature `moonshot-planner`, [M], ships
  OFF): the co-optimizer — exact per-track multiple-choice-knapsack DP
  (the tropical budget DECOMPOSES per track; within 2% of the brute
  oracle from bucket rounding only, never loses to V1 — 80 fixtures);
  convex water-filling (KKT bisection) for rate-based models,
  cross-checked by CMA-ES to 1e-16; BIPOP-CMA-ES for non-convex models
  (MEASURED rejection: single-run CMA-ES converges AWAY from
  activation cliffs — surrogate-threshold fixture stuck at spend 0
  with 3.6× the error until BIPOP restarts crossed it). Scoreboard:
  V2 beats-or-ties V1 AND hand allocation on 25/25 fixtures —
  NECESSARY promotion evidence; the flagship-set gate is huq.15's
  Gauntlet call and this feature stays off until it passes.

## Invariants

1. Fits are pure functions of the observation multiset — identical ledger
   snapshots give identical models (P2; arrival order is irrelevant).
2. Predictions always carry uncertainty (bands + n + extrapolation flag);
   thin data refuses structurally.
3. `ErrorLedger::total()` bounds the sum of true per-stage errors whenever
   each contribution bounds its stage (additive conservativeness — the
   fixture-pipeline law).
4. Attribution is complete by construction: unattributed mass must be
   declared residual, and the lint refuses invalid mass.
5. Explanation payloads are valid JSON for valid and invalid internal state.
   Invalid ledgers produce a fail-closed diagnostic object; they never emit
   `NaN`, infinities, or Rust `Debug` option syntax as evidence.
6. Tune evidence is scoped and transactional: no sample from a different
   machine, shape key, producer schema, operation, build, or artifact lineage
   can influence a returned model, and any failed validation returns no model.

## Error model

`CostRefusal` (insufficient data, bad input, observation/evaluation limits,
nonfinite arithmetic, empty evaluation), `PlanOracleError`, `TuneModelError`,
`LedgerDefect` (bad contribution/residual/aggregate), `TimeLedgerDefect` (bad
stage or aggregate), `PlanInputError`/`AllocationError`, and feature-gated
`VoiError` are structured and teaching; none panic across the boundary.
`TuneModelError` preserves ledger I/O errors and distinguishes absence, schema
failure, scope mismatch, and numerical refusal.

## Determinism class

Deterministic: sorted cost-model observations, nearest-rank quantiles with
deterministic tie-breaking, and BTreeMap iteration. No RNG anywhere. VoI uses
a fixed sweep and ranking order, canonical evaluation-permit ordinals, and
exact computation receipts while preserving prospective audit append order;
its result is deterministic when the supplied oracle is deterministic and its
bounded metadata and declared work are truthful. VoI identity floats bind by
exact IEEE-754 bits rather than display text; node and audit order are
semantic, ranked-row order is semantic, and source probe-menu input order is
the documented canonicalized exception. Ledger receipt bytes use fixed tagged
binary encodings; duplicate receipt inputs canonicalize through `BTreeSet`,
and receipt producer operations run in `ExecMode::Deterministic`.

Numeric scope: SAME-ISA (determinism-tier policy, bead frankensim-lyms).
Cost-model fitting and the moonshot allocator call platform libm
(`ln`/`exp`/`powf`), which is not correctly rounded and may drift last-ULP
across ISAs or libm versions. Identity/receipt BYTES above are exact and
ISA-independent; fitted cost predictions and allocations are bitwise
reproducible only on one ISA/libm build. Promote by routing through
`fs_math::det` and registering in `check-libm` if cross-ISA cost replay is
ever claimed.

## Cancellation behavior

Library-owned work is bounded pure computation or bounded ledger reads. Cost-model
history/evaluation, oracle edges/errors, receipt bytes/depth/nodes/container
items, op fields, retained artifacts, and tune sample counts each have explicit caps. VoI sweeps are bounded
by `MAX_VOI_EVALUATIONS`; audits and consumed snapshots are bounded by
`MAX_VOI_AUDIT_RECORDS` and `MAX_VOI_SCHEDULED_CONTEXTS`; the fixture oracle has
its own Cartesian-work cap. Canonical receipt and outcome artifacts are bounded
by `MAX_VOI_LEDGER_RECEIPT_BYTES`, and a complete probe catalog is bounded by
`MAX_VOI_PROBE_CATALOG_BYTES`; decision/model identity artifacts are fully
content-hash-checked under `MAX_VOI_AUTHORITY_ARTIFACT_BYTES`. Bounded ledger
reads precede receipt materialization. VoI admits the exact evaluation count and exact
declared-work charge against both public caps and the caller's `DecisionBudget`
before the first oracle call. Every call receives a canonical private permit
and has a `Cx` checkpoint before and after it; direct `DecisionOracle`
implementations must also checkpoint long-running internal work and refuse an
insufficient charge. `LiveDecision` remains invocation-bounded only: it cannot
preempt or enforce a time or memory limit on an arbitrary synchronous closure.

## No-claim boundaries

- A structurally and cryptographically linked roofline row proves provenance
  and internal measurement consistency, not that the machine was free of every
  unmodeled environmental disturbance. The producer's baseline/admission
  policy remains the authority for that claim.
- Repetitions in one receipt are independent timed observations for empirical
  residual bands, but they cover one problem size and one run environment. The
  model marks other sizes as extrapolation; it does not claim cross-shape or
  cross-machine transfer.
- `PlanCostOracle` costs are predictive empirical estimates, not worst-case
  wall-time certificates. Its observed error maximum is retrospective; the Rep
  Router therefore uses it only to enlarge an uncertified declaration.
- The sealed carrier proves WHO minted a model and from WHAT validated
  row; it does not re-verify the ledger at prediction time (staleness
  beyond `recorded_at_ns` is the consumer's freshness policy, not yet
  folded into the class), and `provisional_unaudited` is a labeling
  mechanism, not a sandbox — a consumer that ignores
  `CostEvidenceClass` forfeits the distinction (fs-ir admission and
  fs-session estimate do not). External issuer signatures over sealed
  scopes remain future scope, coordinated with the admitted-scientific-
  color lane (bead 6pf9).
- VoI work units are oracle-declared accounting, not measured wall time,
  instructions, allocations, or proof of cooperation. A direct oracle that
  ignores `Cx`, understates work, blocks, or panics, and an arbitrary
  `LiveDecision` callback with those behaviors, remain outside the typed
  cancellation and resource guarantee; refusals produce no ranking authority.

## Unsafe boundary

None. Safe Rust only.

## Feature flags

- `voi-queries` [F] (default OFF) — value-of-information query planning
  (knh1.6, Proposal C); gates the `voi` integration target and its optional
  fs-eproc authority dependency.
- `moonshot-planner` [M] (default OFF, bead gp3.9) — the co-optimizing
  allocator + fs-dfo dependency; promotion gated on the huq.15 flagship
  Gauntlet (the fixture-matrix scoreboard here is necessary, not
  sufficient).

## Conformance tests

`tests/conformance.rs`: fixture-pipeline conservativeness with tightness
tracked (2000 adversarial draws), completeness-lint refusals, held-out
quantile-band calibration (coverage logged), online-update improvement
under cost drift, exact-scope/legacy-schema tune refusals, LIVE Rep Router
replanning from fitted models, Time Ledger attribution + calibration. The
fs-roofline producer integration records a real production receipt-v3 row and
rebuilds a three-sample exact-key model through this crate at the exact 1 MiB
dependency-receipt ceiling; the paired retained cap+1 case returns
`LedgerArtifactReadLimit` before fitting or payload materialization. Unit tests cover
refusals, stable degenerate fits, band ordering, extrapolation, transactional
invalid observations, exact caps plus limit+1, empty evaluation, hostile
duplicate-key receipts, and producer-statistic rederivation.

`tests/alloc_battery.rs` covers budget-safe allocation, typed input/work
refusals, online re-planning, oracle bounds, evaluator safety, and tropical
composition. Feature-gated `tests/voi.rs` covers exact boundaries and limit+1,
callback/domain refusals before evaluation, cancellation at every evaluation
position, ambient `Cx` budget exhaustion before callback entry, zero-callback
evaluation/work preflight refusals, permit ordinal/remainder accounting, exact
computation receipts and budget-bound provenance, exact target resolution,
probe economics, sealed menu/context identity, asymmetric subset contraction,
structured estimated hints, chronological-order e-process counterexamples,
bounded append-only matched-cost audits, live activation/demotion, policy and
snapshot isolation, concurrent duplicate-spend refusal, and cumulative monotone
budget arithmetic. The ledger-authority cases create a real VCS snapshot and
sealed receipt operations, authenticate a complete catalog and prospective
matched outcomes, reverify all authority at scheduling, and fail closed for a
forged locator, model/policy drift, expired context, swapped outcome artifacts,
and superseded ranked head. Owner-unit identity batteries additionally change every
registered semantic field, including counts, sequence order, and exact one-ULP
float-bit mutations whose chosen fixed-precision display text remains
identical; an explicit non-movement case reverses the source menu.
The integration lane admits current retained versions and refuses both stale
and future ranked-source, ranked-menu, and audit-context versions.

## No-claim boundaries

- Error contributions are attribution bookkeeping over estimates or
  certificates; RIGOROUS enclosure composition lives in
  fs-evidence/fs-ivl (each entry's `Rigor` says which you have).
- Cost features are single-scalar (size); multi-feature quantile
  regression is future scope.
- The strict receipt decoder intentionally accepts only current
  fs-roofline production receipt-v3 / row-v4 evidence. It is not a generic
  JSON API and does not provide compatibility authority for older schemas.
- No DWR / conformal / MLMC estimator inputs yet — they wire in as their
  owning crates land (the enum slots exist).
- Greedy V1 is not an exact optimizer; the fixture oracle is bounded and is
  evidence, not a production planner. Co-optimizing V2 remains `[M]`, default
  OFF, and makes no flagship superiority claim before the huq.15 gate.

## No-claim boundaries (voi)

- MYOPIC one-step VoI only — sequential VoI is intractable and
  deliberately not offered (there is no tree API to misuse).
- The reported flip fraction is a MIDPOINT-GRID ESTIMATE under the uniform
  interval measure, not a certified probability. `RankedMenu::grid` and every
  hint retain that qualifier; a sampled zero cannot authorize a universal
  zero-value claim. Adaptive or regularity-certified quadrature arrives later.
- One-node-at-a-time sweeps ignore uncertainty INTERACTIONS; joint
  flips are underestimated — documented, not hidden.
- Probe `shrink` factors are menu declarations, not measured posteriors. Their
  endpoint-wise contraction remains a subset with the declared width factor,
  but the prospective audit is what tests their empirical value.
- The library bounds and validates surrogate calls but cannot prove an
  arbitrary callback pure; VoI determinism and replay require the caller to
  supply the declared cached deterministic margin. A callback panic is also
  outside the typed-refusal contract and propagates to its owner.
- Raw `RankedMenu`, `MatchedAuditRecord`, `observe_audit`, and `schedule` calls
  remain caller-supplied advisory data. Only independently verified wrappers
  passed through `observe_authenticated_audit` and `schedule_authenticated`
  carry this crate's ledger-backed authority.
- An authenticated ranked receipt proves exact agreement with the retained
  decision/model/version envelopes, exact supplied complete-catalog artifact,
  exact ranked roots, and the VCS registry's current branch head. It does not
  prove arbitrary callback code pure or correct, nor prove that an external
  universe contains no probe omitted before the producer declared the catalog
  complete. The current `Vcs` commit registry is process-local rather than a
  portable external signature service.
- An authenticated audit receipt proves that the exact prospectively ordered
  record agrees with the two retained outcome artifacts and sealed ledger
  lineage. It does not prove that a producer reported physical reality
  truthfully, that trials were statistically independent, or that no
  postselection occurred before the ledger artifacts were created.
- Expiry is enforced against a caller-supplied logical day (including a
  not-before check), not an authenticated wall clock. One scheduler prevents
  duplicate observation IDs and duplicate snapshot spend locally, but receipts
  can still be replayed into a second process; cross-process exactly-once spend
  and external issuer signatures remain future session/ledger policy work.
