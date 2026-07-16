# CONTRACT: fs-voi

Value-of-information and active validation: the strategic layer deciding what
information to acquire next — spend where it can CHANGE A DECISION.

## Purpose and layer

Layer L4 (decision/optimization). No dependencies — Gaussian decision algebra
with an in-house normal CDF (`erf`). Pure, deterministic.

## Public types and semantics

- `Uncertainty { numerical, statistical, model }` — `total_std` (quadrature),
  `dominant` component.
- `DesignEstimate { name, mean, uncertainty }` — an Evidence-carrying objective
  estimate (minimizing).
- `ranking_flip_probability(chosen, other)` — `P(other actually better)` (Φ).
- `decision_posture(&[DesignEstimate])` — best, runner-up, flip probability.
- `evpi(&[DesignEstimate])` — expected opportunity loss of the current top-two
  decision (0 when robust; positive when close).
- `evpi_by(len, mean_at, std_at)` — the allocation-free accessor-driven
  form of `evpi` (bead sj31i.62): the SAME top-two scan and pairwise
  opportunity-loss computation through one shared code path, so results
  are bitwise-identical to `evpi` over the equivalent slice. Non-finite
  means are skipped and equal-mean ties break toward the LOWER index;
  callers own their index order (a canonically ordered menu gets
  canonical tie-breaking with no clone and no sort). `std_at` is
  consulted only for the final top-two pair. `EVPI_SEMANTICS_VERSION`
  is NOT bumped by this addition: the scan and arithmetic are shared
  with, and bit-identical to, the existing `evpi`.
- `ActionKind` (Surrogate / Simulate / Refine → numerical; Sample → statistical;
  Test → model) + `Action { name, kind, target_design, reduction, cost }`.
- `action_value(&[DesignEstimate], &Action) -> ActionValue` — the EVPI reduction
  per cost; ~0 for a decision-irrelevant target. A zero-cost action with
  positive decision value has infinite value-per-cost; negative or non-finite
  costs are not recommended.
- `recommend(&[…], &[Action], stop_threshold) -> Recommendation` — the best
  decision-value-per-cost action, or STOP when EVPI ≤ threshold.
- `heuristic_choice(&[…], &[Action])` — the uncertainty-proportional baseline
  ([M]) VOI must beat.

## Invariants

- ESTIMATOR-vs-DECISION: an action on a design outside the decision boundary
  buys ~0 EVPI reduction, however uncertain that design is.
- STOP is returned exactly when the current EVPI ≤ the stop threshold (a robust
  decision), never prematurely below it.
- When the decision-boundary uncertainty is MODEL-dominated, a model-reducing
  action (Test) beats a statistical one (Sample) — decision-aware escalation.
- Ranking-flip probability is the Gaussian `Φ` of the standardized mean gap.
- Non-finite means are excluded from top-two boundary selection so malformed
  estimates cannot displace finite decisions; if fewer than two finite means
  remain, no decision boundary is reported.

## Error model

Total functions; no panics (degenerate zero-σ cases handled explicitly).

## Determinism class

Fully deterministic: every quantity is a pure function of the estimates +
actions.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/voi.rs` (20 cases): ranking-flip probability vs separation; full
opportunity loss zero when robust / positive when close + posture; the
high-variance third-alternative falsifier (surrogate ~0, full material, no
robust STOP); two-design full-vs-closed-form agreement to quadrature
resolution; full dominates the surrogate on every menu (incl. the three-way
exact tie against the order-statistic constant); permutation invariance;
dominated alternatives contribute nothing; a seeded deterministic Monte Carlo
oracle within its confidence bound; near-√MAX uncertainty composition stays
finite end to end; bitwise power-of-two scale equivariance across ~120 orders
of magnitude; subnormal menus; information on a decision-irrelevant design is
worthless; STOP for a robust decision; VOI beats the uncertainty-proportional
baseline; model-fidelity escalation; zero-cost per-cost ranking; invalid-cost
refusals; non-finite mean exclusion; determinism; semantics-version lock.

## No-claim boundaries

- v2 (bead sj31i.5): [`expected_opportunity_loss`] includes EVERY alternative
  via a survival-product quadrature (`E[min] = L + ∫ Π S_j`) at a fixed
  [`EOL_QUADRATURE_PANELS`]-panel Simpson rule with `12σ` tail truncation. It
  is an ESTIMATED value: quadrature resolution is `(U−L)/panels`-scale
  curvature error, NOT a certified enclosure, and menus mixing deviation
  scales wider than ~2 orders of magnitude inside one window are resolved at
  the coarsest scale. The renamed [`top_two_evpi_surrogate`] survives only as
  an action-ranking baseline and must never gate a global robustness claim.
- Uncertainty composition and pairwise deviation sums are overflow-safe scaled
  norms; power-of-two rescaling of a menu is bitwise-exact through the full
  evaluator.
- Objectives are treated as Gaussian; heavy-tailed / correlated estimates are a
  refinement.
- Action cost models are SUPPLIED (fs-plan-models); this crate arbitrates value
  per cost. Sequential-decision validity (fs-eproc) and the HELM planner
  consuming the rankings are downstream integrations.
- `reduction` is the fractional uncertainty cut an action achieves; calibrating
  it from real action outcomes is later work.
