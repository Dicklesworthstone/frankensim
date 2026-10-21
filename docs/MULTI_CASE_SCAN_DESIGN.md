# One heating design across multiple physical cases

`fs_ascent::scan_design::constraints::value::scenarios::ScanScenarioStudy`
evaluates one authored source design against an explicit, finite case set.
Each `ScanScenario` owns an existing `ValueConstrainedScanStudy`, a unique name,
and a positive weight. Cases can differ in their mesh, original material charts,
initial state, boundary history, measurement targets and physical requirements.
They must have identical source templates, selected control entries, coordinate
maps, reference scales, initial coordinates and source-validation time. Matching
coordinate counts alone is not sufficient.

`ScanScenarioObjective::WeightedMean` returns the normalized weighted mean of
the **original case losses**; `WorstCase` returns their maximum. Weights are
caller-declared preferences, not inferred probabilities. All losses must be
finite and nonnegative. Relative weights which underflow are refused rather
than silently deleting a case. The mean scales both losses and weights before
summing so finite large values need not overflow their unnormalized sums.

Every constraint row survives, unweighted and in original case/row order.
`constraint_offsets` includes a final sentinel identifying each case's slice.
Do not constrain an average temperature or average violation as a substitute:
one specimen at 157 J/kg and another at 147 J/kg average to 152 J/kg, but the
first still violates a 155 J/kg limit. Robust measurement losses likewise do
not soften the physical requirements. Reports retain each case's original
physical report, predictions, raw residuals, constraints and smoothing margins.

## Execution and work accounting

Construction invokes no physical factory. `evaluate_at` calls the existing
producer independently for every case; no thermal surrogate or derivative tape
is introduced. A failed physical solve is a typed `Physical { case, source }`
error. No partial aggregate escapes, no failed case is dropped, and no error is
converted to a penalty. Work already performed remains charged to its original
case. Retrying a failed candidate reruns its cases, rather than caching a
partially completed case vector.

`value_batch` borrows all case declarations for concurrent **candidates**. Cases
within each candidate run in declaration order. Each uses its worker's own Cx
and physical state. The original per-case batches reserve their allowances
before dispatch, atomically count started physical reservations, and settle
unused reservations on drop after workers drain. Mutable-only factories remain
available through the serial interface; concurrent callers require `Fn + Sync`
factories and thread-shareable observations. Forgetting a batch conservatively
retains its reservations. No case receives a fresh budget merely because it
runs on another thread.

`evaluations` and `evaluation_limits` return one counter per case. Failed late
cases can make those counters differ. `remaining_evaluations` is the minimum
remaining capacity across cases. One complete candidate costs one trajectory
**per case**, not one trajectory for the ensemble. Budget extension validates
all requested totals before changing any of them. Existing probes are retained.

The aggregate budget limits case count, name bytes, constraint rows and retained
final-field vertices. A vertex count is not a byte/RSS bound; each physical
producer still owns its mesh, nonlinear work and source-projection allowances.
The retained-vertex limit is checked as final reports become available; it does
not preflight allocation inside a physical callback.

## Scope and evidence

This is finite-case, fixed-grid design. Satisfying every supplied case does not
prove feasibility between cases, coverage of a probability distribution, or
bounds on unsampled spatial/time peaks. Neither objective policy supplies a
stationarity or global-optimality certificate.

Two unit tests cover scaled aggregation and malformed scalar inputs. Four
physical regressions in `examples/scan_power_fit/tests/scenario_design.rs`
cover original-report parity, unaveraged violations, mismatched control maps,
failed/cancelled work, budget extension, and concurrent batch equivalence.
These Rust tests were authored but not executed in the authoring environment:
Cargo/rustc and DSR/RCH were unavailable. Independent Python/Decimal aggregation
checks are not execution of the Rust implementation.
