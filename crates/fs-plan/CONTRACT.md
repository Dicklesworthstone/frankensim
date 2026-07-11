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
Runtime deps: `std`, fs-geom, fs-ledger.

## Public types and semantics

- `CostModel` — log-log power-law fit (`cost ≈ exp(a)·size^b`) with
  EMPIRICAL residual quantiles supplying P10/P50/P90 bands; predictions
  carry observation counts and an extrapolation flag (an estimate is
  itself evidenced). Below `MIN_OBS` it refuses (`CostRefusal`) instead of
  guessing. `observe()` refits online; `calibration()` audits held-out
  band coverage; `median_rel_error()` is the improvement metric.
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
  plans with THIS machine's measured history (per-edge models at
  registered reference sizes; recorded actuals feed refits; observed
  error p90 backs `measured_error_abs`).
- `cost_model_from_tune` — rebuilds a kernel's model deterministically
  from fs-ledger `tune` rows (fs-roofline's recorded rates).

- `voi` module (addendum Proposal C, bead knh1.6; [F], behind
  `voi-queries`): the ignorance market, v0 as a RANKED LIST.
  `UncertaintyNode` (interval + nominal), `LiveDecision` (a cached
  surrogate threshold verdict; `flip_probability` sweeps one node's
  interval at grid cost — near-free by construction),
  `Probe`/`ProbeKind` (the priced menu UNIFYING computational and
  physical evidence), `rank_purchases` (MYOPIC one-step VoI:
  flip-probability-per-dollar, deterministic tie-breaks),
  `hint_for_query` (the Proposal-8 anytime hint, now decision-priced),
  `schedule_probes` (greedy affordable top-k — the discrepancy-probe
  scheduler), `audit_verdict` (the kill criterion as code: VoI keeps
  scheduling authority only while recommended purchases measurably
  outperform agent-chosen alternatives; no evidence → no authority).
  Decision, sweep, ranking, and scheduling entry points are fallible:
  arity, nonempty bounded collections/names, unique identities, finite
  ordered intervals, nominal containment, callback margins, target
  resolution, grid and aggregate evaluation work, probe economics,
  ranked values, and budgets are validated before evaluation or spend.
  `RankedPurchase` construction is sealed behind ranking with read-only
  getters. Hinting and scheduling reapply the canonical
  score-descending/cost-ascending/name-ascending comparator, so caller
  order carries no authority. Duplicate ranked identities refuse rather
  than buying twice, and a positive cost must strictly decrease a finite
  remaining budget.

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

## Error model

`CostRefusal` (insufficient data, bad input), `LedgerDefect` (bad
contribution/residual/aggregate), `TimeLedgerDefect` (bad stage or aggregate),
`PlanInputError`/`AllocationError`, and feature-gated `VoiError` are structured
and teaching; none panic across the boundary. Ledger I/O errors propagate as
`fs_ledger::LedgerError`.

## Determinism class

Deterministic: sorted observations, nearest-rank quantiles with
deterministic tie-breaking, BTreeMap iteration. No RNG anywhere. VoI uses a
fixed sweep and ranking order; its result is deterministic when the supplied
cached surrogate callback is itself deterministic.

## Cancellation behavior

All calls are short pure computations or single ledger reads. VoI sweeps are
bounded by `MAX_VOI_EVALUATIONS`; the fixture oracle has its own Cartesian-work
cap. No `Cx` integration is needed at this layer.

## Unsafe boundary

None. Safe Rust only.

## Feature flags

- `voi-queries` [F] (default OFF) — value-of-information query planning
  (knh1.6, Proposal C); gates the `voi` integration target.
- `moonshot-planner` [M] (default OFF, bead gp3.9) — the co-optimizing
  allocator + fs-dfo dependency; promotion gated on the huq.15 flagship
  Gauntlet (the fixture-matrix scoreboard here is necessary, not
  sufficient).

## Conformance tests

`tests/conformance.rs`: fixture-pipeline conservativeness with tightness
tracked (2000 adversarial draws), completeness-lint refusals, held-out
quantile-band calibration (coverage logged), online-update improvement
under cost drift, deterministic rebuild from ledger tune rows, LIVE Rep
Router replanning from fitted models, Time Ledger attribution +
calibration. Unit tests cover refusals, band ordering, extrapolation
flags, and arrival-order determinism.

`tests/alloc_battery.rs` covers budget-safe allocation, typed input/work
refusals, online re-planning, oracle bounds, evaluator safety, and tropical
composition. Feature-gated `tests/voi.rs` covers exact boundaries and limit+1,
callback/domain refusals before evaluation, exact target resolution, probe
economics, duplicate scheduling identities, and monotone budget arithmetic.

## No-claim boundaries

- Error contributions are attribution bookkeeping over estimates or
  certificates; RIGOROUS enclosure composition lives in
  fs-evidence/fs-ivl (each entry's `Rigor` says which you have).
- Cost features are single-scalar (size); multi-feature quantile
  regression is future scope.
- `json_f64_field` is a scanner for OUR canonical flat payloads, not a
  JSON parser (fs-ir owns real JSON).
- No DWR / conformal / MLMC estimator inputs yet — they wire in as their
  owning crates land (the enum slots exist).
- Greedy V1 is not an exact optimizer; the fixture oracle is bounded and is
  evidence, not a production planner. Co-optimizing V2 remains `[M]`, default
  OFF, and makes no flagship superiority claim before the huq.15 gate.

## No-claim boundaries (voi)

- MYOPIC one-step VoI only — sequential VoI is intractable and
  deliberately not offered (there is no tree API to misuse).
- The flip probability uses the UNIFORM measure over the node's
  interval (v0's declared prior); posterior-weighted sweeps arrive
  with the Proposal-3 inventory integration.
- One-node-at-a-time sweeps ignore uncertainty INTERACTIONS; joint
  flips are underestimated — documented, not hidden.
- Probe `shrink` factors are menu declarations, not measured
  posteriors; the prospective audit is what keeps the menu honest.
- The library bounds and validates surrogate calls but cannot prove an
  arbitrary callback pure; VoI determinism and replay require the caller to
  supply the declared cached deterministic margin. A callback panic is also
  outside the typed-refusal contract and propagates to its owner.
- Ranked rows cannot be forged or made authoritative by reordering, but a
  caller passing a raw slice can omit valid rows or splice unique sealed rows
  from different ranking requests. Scheduling authority assumes one complete
  `rank_purchases` output (or ledger evidence retaining its context and
  completeness); omission/splice detection requires a future sealed
  ranked-menu receipt.
