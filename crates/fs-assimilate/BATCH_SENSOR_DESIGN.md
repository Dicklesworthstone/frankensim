# Joint sensor selection and real acquisition

`SharedNoiseStream::plan_probe_batch` selects several measurements while
accounting for what earlier selections would already teach the model. It uses
the retained physical/nuisance posterior, signed calibration sensitivities and
each sensor's independent instrument noise. It does not take measured values.

This extends the one-step `rank_probes` API rather than changing its semantics.
Two independently high-scoring sensors may be redundant. For a synthetic prior
`diag(2, 1)`, target `x0 + x1`, noise variance `0.01`, and three sensors observing
`x0`, `x0`, `x1`, independent top-two selection leaves target variance about
`1.00499`. Conditional selection chooses one sensor per component and predicts
`0.01985124`. These are numerical Gaussian-model results, not accuracy claims
about real instruments.

## Selection

Wrap a `SharedNoiseProbe` in `CostedProbe::new(probe, cost)`, using a positive
`u32` cost. All costs and the `u64` acquisition limit must share one explicit
caller-chosen unit. `ProbeBatchBudget::new(max_readings, max_cost)` supplies
both limits. Neither is the numerical-compute budget carried by `Cx`.

```rust,ignore
use fs_assimilate::groups::stream::design::batch::{CostedProbe, ProbeBatchBudget};

let plan = stream.plan_probe_batch(
    &target,
    &candidates, // checked CostedProbe values; no measured outcomes
    ProbeBatchBudget::new(4, 20)?,
    cx,
)?;
for step in plan.steps() {
    println!("{}: {} -> {}",
        step.candidate().probe().row().as_str(),
        step.prior_target_variance(),
        step.posterior_target_variance());
}
```

At each step the largest conditional target-variance reduction per cost wins,
with row-ID ties. Remaining candidates are reconditioned, not merely re-ranked
using their original scores. The plan retains selected declarations, source
posterior history, goal, acquisition budget, conditional variances, and actual
stop reason. Reading count, available stream row capacity, acquisition cost,
exhausted candidates and zero computed information are distinct stops.

The numerical method projects through `sqrt(D) U^T` once, then propagates
square-root Joseph residual rows. Each selected measurement contributes one
new independent-noise coordinate shared by every residual row. Omitting or
mis-signing these coordinates loses conditional cross-covariance. No fabricated
observations, mean updates, dense covariance inverse, diagonal approximation,
or prior-minus-gain subtraction is used to produce the final variance.

The selector is a GREEDY scalar-goal heuristic. It does not certify a globally
optimal subset, an approximation ratio, EVPI, calibration quality, or attainment
of an accuracy requirement. Noise must be mutually independent conditional on
the retained nuisance model. Free measurements, nonlinear relinearization,
dynamic propagation and arbitrary additional pairwise noise are outside this
API. Replan when the posterior or measurement model changes.

## Acquire separately, then publish atomically

Call `plan.check_current(&stream)` before external acquisition. It reserves
nothing. Supply actual values as `ProbeReading::new(row_id, value, instrument)`;
then call `plan.assimilate_readings(&mut stream, &model, &readings, cx)`.

Reading order need not match selection order: IDs attach the correct value to
the exact planned response, Jacobian and independent variance. Returned update
indices refer to selection order, not the incoming reading slice. Missing, extra,
duplicate or unselected rows reject the whole arrival. The API rechecks the
source posterior history, model, physical dimension and remaining row capacity.
The existing checked `JacobianBatch` and transactional stream update perform the
actual assimilation; predicted covariance never substitutes for measured data.

Cancellation, numerical refusal or invalid input does not consume stream rows
or change its posterior/history. Borrowed readings remain available for retry.
A successfully executed nonempty plan is then stale; it cannot silently apply
its old forecast to the new posterior. Exact stream-level replay still exists
in the ordinary stream API. Empty plans accept only empty readings.

The atomic boundary is numerical, not physical: external measurements and their
costs cannot be undone by a failed numerical update. This API does not acquire
hardware, persist observations, authenticate sensors, or issue refunds.

## Bounds and tests

At most 256 candidates and 256 selections; existing joint-state and nuisance
limits apply. Replay protection is never evicted. Numerical selection uses one
cumulative work/poll accountant. Its core work is `O(m d² + k m (d+k))` and
scratch is `O((m+1)(d+k))` for `m` candidates, `k` choices and joint dimension
`d`; bounded identity admission additionally checks candidate pairs. Acquisition
preparation, covariance admission and stream assimilation retain their existing
separate checked compute scopes. Unresolved range/contraction refuses rather
than clamping or claiming certainty.

```sh
cargo test -p fs-assimilate --test shared_design --test shared_noise
cargo run -p fs-assimilate --example shared_batch_design
# Two caller-supplied readings for the example's selected a and c rows:
cargo run -p fs-assimilate --example shared_batch_design -- 3.2 0.9
```

The example uses a labeled synthetic unitless model. Its observation values are
kept separate from selection, and it compares predicted and acquired target
variance through the public APIs. Tests cover redundancy, acquired-prefix
agreement, shared-error effects, permutation/unit changes, count/cost/row caps,
near-exact measurements, row-addressed acquisition, stale plans, and retry.
