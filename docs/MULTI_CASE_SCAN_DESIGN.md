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

## Population optimization and executable workflow

`scan_power_fit::scenarios::ScenarioPopulationStudy` connects the evaluator to
`ConstrainedDifferentialEvolutionState`. The existing DE implementation still
owns sampling, mutation, feasibility ranking, frozen generations, parallel
publication and successful-sibling caches. Every case's rows enter selection.
A failed physical case remains an error, not an infeasible/penalty score.

The driver requires enough remaining physical work in **each case** for every
remaining search attempt plus one replay. Failed partial candidates retain their
individual counts. Serial and parallel advances bound dispatch by the least
remaining case allowance, excluding that replay. Extension validates the whole
set of search/physical totals before mutating them. Parallel execution uses the
existing TilePool, worker-local Cx and all original case batches; candidates,
not cases within a candidate, run in parallel. There is no extra nested pool.

`replay_best` runs every case again and requires the aggregate objective and
feasibility score to match their selected values bit-for-bit before returning
the common sources and full case reports. The DE score does not store earlier
full fields, so this is not a bitwise comparison of every old nodal value.
Factories must retain their original deterministic model semantics.

`aggregate_target_met` refers to the search's declared aggregate loss threshold.
`every_case_target_met` additionally requires every original case loss to meet
that threshold. Under mean loss, the former must not be called recovery of every
case. Neither says the sample set covers a continuous uncertainty distribution.

Run the synthetic example from the repository root:

```bash
# Single starting state: recover the familiar approximately 30 W / 20 W fit.
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example scenario_scan_fit -- --cases 80

# Two starting states, with ALL physical limits retained; parallel candidates.
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example scenario_scan_fit -- --cases 80,90 --workers 4 60 0.38 2400

# Minimize worst-case original loss instead of weighted-mean loss.
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example scenario_scan_fit -- --cases 80,90 --worst-case 60 0.32 2400
```

`--cases` gives initial specific enthalpies in J/kg, `--weights` gives matching
positive case preferences, and positional arguments are the per-case absorbed
energy cap in joules, sampled peak liquid mass fraction, and search-attempt cap.
Default cases are 80 and 90 J/kg with equal weights. Each original specimen is a
unit-volume/density cube with a 100 J/kg latent plateau, heated uniformly for
two one-second intervals. Both cases target liquid mass fractions 0.1 then 0.3.
Their exact targets cannot be met simultaneously by one common schedule.

Consequently the two-case default prints its feasible compromise and exits
unsuccessfully for missed individual targets. This is deliberate: exhausting
a search, producing a feasible candidate, and exactly recovering every target
are three different outcomes. The executable reports all of them without an
optimum claim. A single-case invocation can recover its targets.

Independent analytical phase-balance/Python DE checks give approximately:

| Policy / limits | Common watts, first then second second | Aggregate loss |
| --- | --- | --- |
| Equal-weight mean, loose limits | 25, 20 | 0.25 |
| Mean with cold:warm weights 1:3, loose limits | 22.5, 20 | 0.1875 |
| Worst case, loose limits | 25, 20 | 0.25 |
| Worst case, peak envelope capped at 0.32 | 28.2079, 13.5841 | 0.35291 |

The cold-only 30/20 W schedule produces final fractions 0.3 and 0.4 in the cold
and warm cases. Its warm envelope is approximately 0.40207944, violating a 0.38
cap even though the cold case passes. The joint 25/20 W design produces 0.25
and 0.35, with warm envelope approximately 0.35207944. These are independent
uniform-P1 analytical references, not native Rust execution or experimental data.

One driver unit test covers replay reservation. Four executable regressions
cover the compromise, worst-case hot-spot limits, serial/parallel continuation,
and command/cancellation admission. Together with the evaluator, eleven Rust
tests are authored. Native compilation/execution remains unverified because
DSR, its documented fallback, RCH, Cargo and rustc were unavailable.
