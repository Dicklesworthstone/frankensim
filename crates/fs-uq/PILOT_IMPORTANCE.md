# Learning rare-event proposals from physical pilots

`PilotImportanceExecution` runs a fixed number of complete adaptation pilots,
then a separate fixed-size production experiment. It reuses `ImportanceExecution`
for every physical call, likelihood ratio, conditional Gaussian map and probability
reduction. It adds no physical surrogate, covariance inverse or new random generator.
The rare-event ceiling, physical model, laws, units, defensive fraction, round
counts and fitting policy must remain fixed for the entire study.

```rust
use fs_uq::{CorrelationModel, ImportanceConfig, ImportanceFitConfig,
    ParameterUncertainty, PilotImportanceConfig, PilotImportanceExecution,
    PropagationMethod, UqPlan, UqStatus};

let plan = UqPlan::new("standardized response", PropagationMethod::MonteCarlo, 10240)
    .with_parameter(ParameterUncertainty::gaussian("load", 0.0, 1.0, "N"))
    .with_correlation(CorrelationModel::Independent)
    .with_compliance_threshold(5.0);
let config = PilotImportanceConfig {
    pilot_rounds: 4,
    samples_per_pilot: 512,
    initial_proposal: ImportanceConfig {
        germ_shift: vec![0.0],
        defensive_fraction: 0.2,
    },
    fit: ImportanceFitConfig {
        elite_fraction: 0.1,
        relaxation: 1.0,
        maximum_shift_norm: 12.0,
    },
};
let mut execution = PilotImportanceExecution::new(&plan, config)?;
let report = execution.advance(usize::MAX, || false, |x| Ok::<_, &'static str>(x[0]));
assert_eq!(report.status, UqStatus::Complete);
assert_eq!(report.model_calls_evaluated, 10240);
// Only the 8192 independent production responses enter this estimate.
let estimate = report.estimate.ok_or("production statistic unavailable")?;
assert!(estimate.failure_probability.is_finite());
# Ok::<(), Box<dyn std::error::Error>>(())
```

This analytical normal-tail fixture illustrates the API, not qualification of a
physical solver or a guaranteed sampling improvement. A real callback can run
an expensive deterministic physical model at the supplied physical parameters.

## Fitting and statistical separation

`ImportanceExecution::fit_proposal` can also be used independently. It requires
a completed successful pilot. Initially it selects at least
`max(2, ceil(N * elite_fraction))` upper-tail responses, including every tie at
the rank cutoff. Once that many strict failures exist, it instead selects the
actual event `QoI > ceiling`; values equal to the ceiling are not failures.
This intermediate cutoff is a heuristic and need not increase between rounds.

The updated unit-covariance Gaussian component mean is
`sum((p/q) * elite * germ) / sum((p/q) * elite)`. The weights correct for the
pilot proposal; an unweighted mean would fit the tilted law instead. Common
scaling avoids overflow in the weight and coordinate sums. The stored scaled
weights have the same common defensive factor, which cancels from this ratio.
Relative contributions can still underflow or round. Declared relaxation blends
the fitted mean with the previous effective shift, then a numerical Euclidean
norm cap limits it. When the proposal equals the target, the previous effective
shift is zero. Diagnostics retain cutoff, selected count, effective sample count,
pre-limit norm, and whether limiting occurred. No covariance fitting is claimed.

Normalization is for fitting only. The final probability remains the ordinary
mean of target/proposal-weighted failure indicators, not a self-normalized or
clipped estimator. No pilot response enters that production mean or error.
Production retains its learned proposal even when its observed error is large.
`Complete` means the entire fixed design finished, not that an accuracy target,
rare-event discovery, variance reduction or an optimal proposal was established.

Phase j uses `seed + j * 0x9e3779b97f4a7c15` modulo 2^64: j=0 is the initial
pilot; j=pilot_rounds is production. The odd multiplier makes these bounded
phase keys distinct even at wraparound. Independence follows the existing
pseudorandom-stream model, not merely the syntax of different integers. Callers
must avoid reusing these keys across separate studies and must not select a
study after seeing favorable production outcomes. Original Gaussian laws,
including singular covariance, are replayed in independent germ coordinates.

The production owner's existing `assess_failure_probability` remains available
through `production()`. Apply it to production alone, with alpha chosen in
advance. Its conservative range-based confidence sequence is distinct from the
descriptive standard error and does not bound solver/model error, finite-grid
bias or floating-point rounding. No new confidence theorem is claimed here.

## Budgets, continuation and limitations

The lifetime cap of one million physical calls includes all pilots. Admission
requires 1..=16 rounds, at least two calls in every pilot and in production,
<=16777216 sample-coordinate pairs per fit and <=67108864 across all pilots.
Ranking and two germ replays use O(N+d) scratch under those limits. Only current
observations and bounded prior-round diagnostics remain retained; a transition
briefly also allocates the next producer. These are observation/replay budgets,
not a bound on the caller's physical solver memory or wall time.

Cancellation brackets physical calls, fitting, owner transitions and final
publication. An unfinished callback retries identical inputs and weights. A paid
pilot survives cancellation during its fit; a paid final production run can be
published on resume without another model call. Zero allowance is a no-op;
positive allowance permits this pure completion work. Cloning resumes the same
experiment. Errors, nonfinite physical outputs and numerical proposal refusals
are terminal and never become zero-valued or skipped observations. No durable
checkpoint codec, adaptive production stopping, or population validation is added.

A single mean shift can stagnate on flat responses, cancel between opposite
failure modes, or target a less important mode. The retained defensive component
protects support under the ideal Gaussian law but cannot guarantee that a finite
budget discovers every mode. Zero observed failures do not prove zero risk.

`tests/product_importance_proposal.rs` checks weighted fitting, ties, singular
maps and cancellation. `tests/product_importance_pilot.rs` checks rare-tail
learning, standalone reconstruction, call/phase continuation, failures and
admission. These tests require native execution before claiming Rust validation.
