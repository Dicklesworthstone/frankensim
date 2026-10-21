# Population exploration followed by physical SQP refinement

`ValueConstrainedScanStudy::begin_sqp_at` borrows an existing physical scan
study and starts the existing `SqpState` at a supplied design. The source
templates, coordinate maps, initial material state, observations, robust losses,
all physical constraints and cumulative evaluation counter remain unchanged.
It does not rebuild an uncharged study or approximate the objective gradient
by finite differences. The original fixed-grid thermal and response adjoints
provide the derivatives, including nonlinear conductivity and radiation.

The original `ScanPowerObjective` must admit the required derivative storage.
Give it nonzero retained-vertex/source/radiation budgets when constructing a
study intended for hybrid search. Value-only evaluations do not allocate that
derivative tape. The handoff separately takes a `ScanResponseBudget` and a
maximum dense KKT dimension; these are derivative-work limits, not new physical
constraints. The study's previously selected checkpointing policy is retained.

```rust,ignore
let mut polish = study.begin_sqp_at(
    cx, &discovered_point, &optimizer_coordinate_bounds,
    response_budget, max_kkt_dimension,
)?;
let report = polish.advance(cx, 1e-8, 80, 100)?;
let replay = polish.replay(cx)?;
drop(polish);
// The original study is available again, with all spent work still charged.
```

Nonempty coordinate bounds must cover the complete original optimizer point.
They are appended as lower/upper SQP inequalities and remain hard evaluation
domain limits. This matters when a search box was not also declared as a
physical `ScanControlBound`: local refinement must not silently lose it.
Original physical rows keep their signs, order and scales. The replay reports
box residuals separately in optimizer-coordinate units.

Starting the session spends one original physical evaluation and its adjoint
work. At least one additional physical evaluation is required for replay.
Each `advance` accepts an additional callback allowance and accepted-iteration
allowance, and also respects the owner's remaining physical allowance minus
one replay slot. Failed derivative probes remain charged. Analytically
out-of-box or negative-power trials are unavailable without running physics;
other physical and derivative errors propagate rather than becoming penalties.

The session retains the native SQP point, multipliers and BFGS approximation
across calls. Cancellation or a failed trial leaves the last accepted SQP point
available. The physical owner cannot be used or have its model replaced while
the session borrows it. Dropping the session releases the owner without resetting
its work counter or rewriting a population optimizer's retained state.

`replay` independently uses the original value-only physical evaluator and
checks the objective and **every** inequality bit-for-bit against the accepted
SQP sample before returning source declarations. Full earlier fields are not
stored in SQP, so this is not a comparison of every old nodal value. Deterministic
factories must retain their original physical meaning. Replay spends its own
physical evaluation; repeated replays do not refund work.

Classical derivatives can be undefined at a phase-chart corner. Such a handoff
returns the original derivative error; it does not invent a one-sided slope.
The value-only owner remains available afterward. Feasibility, matching the
measurement targets, and native KKT convergence are distinct results. The
constraints and derivatives remain fixed-grid and sampled; this is not a global
optimum or continuous-temperature-peak certificate.

Focused regressions cover original SQP row/gradient parity, box rows, a
constrained energy-limited physical fit, work reservation and continuation,
cancellation, and a phase-corner refusal followed by a valid value-only solve.
Native Rust compilation/execution was unavailable during authoring; no compiler,
runtime, or performance pass is asserted.
