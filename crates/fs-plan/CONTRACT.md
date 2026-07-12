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
fs-eproc for anytime-valid audit authority.

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
- `cost_model_from_tune` — rebuilds a model from one EXACT
  `(kernel, shape_class, machine)` fs-ledger key. It accepts only the current
  fs-roofline receipt-v3 / row-v4 production schema through a bounded strict
  JSON parser that rejects duplicate and unknown keys. Every timed sample is
  retained as a same-size observation; median/P25/P75/dispersion and rate are
  rederived from those samples. Embedded kernel/version/run/op/repetition
  identities, the baseline-bound 40-byte machine key, completed production-v2
  op envelope, build identity, result-manifest membership, payload artifact
  bytes/metadata/OUT edge, and dependency-receipt bytes/domain digest/IN edge
  must all agree before the model sees any evidence. Foreign machine/shape
  rows are never scanned or mixed.

- `voi` module (addendum Proposal C, bead knh1.6; [F], behind
  `voi-queries`): the ignorance market, v0 as a RANKED LIST.
  `UncertaintyNode` (interval + nominal), `LiveDecision` (a cached
  surrogate threshold verdict; `flip_probability` sweeps one node's
  interval at grid cost — near-free by construction),
  `Probe`/`ProbeKind` (the priced menu UNIFYING computational and
  physical evidence), `rank_purchases` (MYOPIC one-step sampled VoI:
  flip-fraction reduction per dollar, deterministic tie-breaks) returning a
  sealed complete `RankedMenu`. Its private rows cannot be omitted, spliced,
  or reordered. `source_context_id` binds the caller-declared policy and
  decision snapshot, validated nodes, complete supplied source menu, and grid;
  final `context_id` additionally binds every canonical output row and score.
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
  `BudgetInfeasible` with ranked, VERIFIED relaxations (re-planning at the
  suggested budget succeeds — gated). Measured on the fixture matrix:
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
a fixed sweep and ranking order while preserving prospective audit append
order; its result is deterministic when the supplied cached surrogate callback
is itself deterministic.

## Cancellation behavior

Library-owned work is bounded pure computation or bounded ledger reads. Cost-model
history/evaluation, oracle edges/errors, receipt bytes/depth/nodes/container
items, and tune sample counts each have explicit caps. VoI sweeps are bounded
by `MAX_VOI_EVALUATIONS`; audits and consumed snapshots are bounded by
`MAX_VOI_AUDIT_RECORDS` and `MAX_VOI_SCHEDULED_CONTEXTS`; the fixture oracle has
its own Cartesian-work cap. The VoI limit bounds callback INVOCATIONS, not the
time or memory of arbitrary caller callback bodies. No `Cx` crosses this pure
planning API; callback owners must provide already-budgeted, cached,
cancellation-correct surrogate evaluation.

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
  Router therefore uses it only to enlarge an uncertified declaration. Fully
  authenticated admitted cost-model authority remains tracked separately.

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
rebuilds a three-sample exact-key model through this crate. Unit tests cover
refusals, stable degenerate fits, band ordering, extrapolation, transactional
invalid observations, exact caps plus limit+1, empty evaluation, hostile
duplicate-key receipts, and producer-statistic rederivation.

`tests/alloc_battery.rs` covers budget-safe allocation, typed input/work
refusals, online re-planning, oracle bounds, evaluator safety, and tropical
composition. Feature-gated `tests/voi.rs` covers exact boundaries and limit+1,
callback/domain refusals before evaluation, exact target resolution, probe
economics, sealed menu/context identity, asymmetric subset contraction,
structured estimated hints, chronological-order e-process counterexamples,
bounded append-only matched-cost audits, live activation/demotion, policy and
snapshot isolation, concurrent duplicate-spend refusal, and cumulative monotone
budget arithmetic.

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
- `RankedMenu` identities bind caller-declared policy/snapshot, nodes, source
  menu, grid, and canonical output rows, but they cannot authenticate the
  declared policy/snapshot, identify arbitrary callback code, prove external
  catalog completeness, or prove that a ledger/session snapshot is current.
- Matched audit records are structurally validated, canonically content-bound,
  bounded, and enter one live scheduler in append order, but they are
  caller-supplied rather than ledger-authenticated. A dishonest producer can
  lie, postselect before constructing a scheduler, or replay the same evidence
  into a second scheduler. Ledger-enforced unique audit streams, signatures,
  snapshot freshness, expiry, and independent outcome authentication remain
  required follow-up work (`frankensim-wk4m`); this crate makes no
  authenticated or cross-process exactly-once audit claim yet.
