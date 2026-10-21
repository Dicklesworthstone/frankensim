# Check thermal design limits without requiring gradients

`ValueConstrainedScanStudy` in `fs_ascent::scan_design::constraints::value`
wraps an existing `ScanPowerStudy` with the same `ScanControlBound`,
`ScanEnergyCap`, `ScanStateBound`, and `ScanPeakBound` declarations used by SQP.
It evaluates feasibility without creating an SQP checkpoint or assuming that
the temperature/enthalpy/fraction trajectory is differentiable.

This closes a distinction important at a solidus or liquidus: the physical
state and its constraint values can be valid even when a classical derivative
does not exist. The derivative-based APIs retain their original refusals.

## One trajectory, separate objective and constraints

`evaluate_at(cx, point)` returns `(ScanPowerValueEvaluation,
ScanConstraintValues)`. The fitting objective and all point/peak values share
one original forward trajectory and one physical-evaluation reservation. The
returned constraint report has **no Jacobian or KKT fields**. No source,
material, phase, radiation, energy, or cancellation checks are weakened.

All selected physical powers and geometry controls retain their original
coordinate maps. Authored-energy windows include **unselected and fixed source
values**, not just active optimization coordinates. Source schedules are
integrated with the existing authored-energy implementation; these integrals
remain distinct from the time-discretized thermal solver's heat accounting.

Rows have exactly the SQP signs, scales, and ordering: automatic nonnegative
power rows; physical control lower/upper rows; energy caps; sampled state
lower/upper rows; and sampled-peak caps. Each signed residual must be <= 0.
`is_feasible(tolerance)` accepts a finite nonnegative scaled allowance. For a
physical row with scale `s`, tolerance `t` permits up to `t*s` violation. It is
not a physical solver-error allowance and does not certify continuous extrema.

Robust losses change fitting influence, **not constraint responses**. A low
Huber or least-squares loss does not establish feasibility. Conversely, a
feasible design may still miss its fitting targets.

## Resource and continuation semantics

`ScanValueConstraintBudget` bounds the number of scalar inequalities, stored
authored-energy weights, response count, and retained peak samples. No dense
constraint Jacobian, trajectory-wide tangent tape, transpose solve, or replay
is allocated by this path. The constructor validates declarations without
invoking the physical factory. Actual mesh membership is checked when a run is
constructed, before its first endpoint solve.

For example, a two-power design with one peak over eight vertices and two
endpoints needs `max_rows: 3` (two automatic domain rows plus the peak),
`max_energy_weights: 0` when there are no energy windows, and response bounds
of one response and sixteen peak samples.

`constraints_at` avoids physics only when all attached limits concern inputs.
With any state/peak limit it refuses instead of silently omitting the limit.
`evaluate_at` charges original physical attempts through the wrapped study;
`extend_evaluation_budget` increases that same lifetime allowance.
`sources_at` exports a candidate without a feasibility or optimality claim.

## Lower-level response API

`evaluate_scan_values_with_responses` lives in the conduction objective's
`value::responses` module. `ScanPowerStudy::evaluate_value_with_responses`
exposes it through the same factory/source-coordinate owner. Both return point
and peak values without interpreting them as inequalities. They use the exact
same ordered peak reduction as the full/checkpointed adjoint evaluators.

For sixteen zero liquid-fraction samples and smoothing scale 0.01, the raw
maximum is zero while the envelope is approximately 0.0277258872. A cap of
0.02 therefore fails the envelope constraint. The implementation does not
discard that margin or clamp a fraction envelope to one to report feasibility.

## Scope and native tests

This is a physical constraint evaluator, **not** a feasibility-ranked
differential-evolution adapter. Passing only its fitting loss into an
unconstrained population search still does not enforce these constraints.
The API is suitable for explicit candidate screening and for optimizers that
consume separate inequality values. No infeasible or failed simulation is
converted into an invented objective penalty.

```bash
cargo test --manifest-path examples/scan_power_fit/Cargo.toml \
  --test spatial_design value::
```

Regressions compare all three physical observable families with the adjoint
path, full physical constraint rows with SQP, fixed-heater energy, zero-power
phase-corner values, smoothing margins, declaration limits, cancellation and
shared work accounting. The implementation session could not start DSR, its
documented fallback, RCH, Cargo or rustc because those executables were absent;
native compilation and test success remain unverified.
