# CONTRACT: fs-probe

Discrepancy probes + the budget pie (plan addendum, Proposal 3): makes the
third epistemic color (`estimated`) computable, and makes model-form error
operator-legible.

## Purpose and layer

Layer L3. Depends on `fs-ladder` (L3, the rung registry) and `fs-evidence`
(UTIL, the `Color` lattice). Pure, deterministic functions — no numerical
solver runtime; the physics kernels produce the per-rung states that probes
compare.

## Public types and semantics

- `probe_adjacent(&Ladder, from_rung, coarse, fine) -> Result<DiscrepancyField,
  ProbeError>` — prolongate the `from_rung` coarse state onto `from_rung+1`
  (via the ladder's transfer), take `|fine − prolongate(coarse)|` per
  subdomain. The result is a `DiscrepancyField { kernel, from_rung, to_rung,
  per_subdomain, l_inf, mean, color }` whose `color` is always
  `Color::Estimated` (estimator id + dispersion = `l_inf`). A near-zero gap
  yields a near-zero dispersion (no manufactured error).
- `ErrorContribution { source, color, magnitude }` and `BudgetPie` —
  `BudgetPie::of(&[..])` sums magnitudes by color rank into
  `{ total, verified, validated, estimated }`. `fraction(rank)`, `dominant()`
  (ties resolve to the WEAKER color, conservatively), and `verdict()` — the
  operator-legible string that attributes a model-form-dominated case to the
  closure ("refining the mesh will NOT help") and a numerical-dominated case
  to the mesh/order.
- `ProbeBudget::new(fleet_budget, cap_fraction)` (`cap_fraction` clamped to
  `[0,1]`, non-finite budget inputs fail closed to a zero cap) with `cap`,
  `spent`, `remaining`, and `try_spend(cost)` — a HARD ceiling: spending up to
  EXACTLY the cap is allowed, beyond is refused; bad costs (negative /
  non-finite) are refused.
- `ProbeError` — `Ladder` / `DimMismatch` / `BudgetExceeded` / `BadCost`, each
  with a teaching `Display`.

## Invariants

- A discrepancy field is ALWAYS estimated-color (a model-form probe is never a
  certified bound).
- `dominant()` never over-claims: equal magnitudes report the weaker color.
- `try_spend` mutates the budget only on success; the cap is never exceeded.

## Error model

Structured `ProbeError` values (refusals that teach), never panics.

## Determinism class

Fully deterministic: probes and pies are pure functions of their inputs
(no RNG, no I/O); outputs are bit-identical on replay (probe outputs feed
estimated-color evidence).

## Cancellation behavior

None — all operations are bounded, synchronous pure functions (no `Cx`). A
consumer that probes with an expensive solve runs THAT under its own `Cx`.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/probe.rs` (Proposal 3, 11 cases): near-zero gap manufactures no error;
localized estimated-color discrepancy with a naming estimator; dim-mismatch
and top-rung structured errors; probe determinism; budget pie points at the
closure when model-form dominates and at the mesh when numerical dominates;
empty budget has no dominant; ties resolve conservatively; probe-budget cap is
a hard ceiling at the exact limit; bad costs refused; fraction clamping.

## No-claim boundaries

- Probes CONSUME the fidelity-ladder registry (`fs-ladder`) and the color
  lattice (`fs-evidence`); they do not own rung declarations or run solves.
  The coarse/fine per-rung states are supplied by the physics kernels.
- The discrepancy is a between-rung MODEL-FORM estimate, not a certified error
  bound; it is deliberately estimated-color.
- Probe SCHEDULING by expected information (Proposal C value-of-information) is
  a separate bead; this crate provides the budget CAP mechanism, not the
  scheduler that picks which probes to run.
- The `Refine1d`-based comparison inherits `fs-ladder`'s demonstrator transfer;
  real per-kernel prolongation is injected there.
