# Multiple physical responses from one uncertainty sample

`MultiQoiExecution` evaluates one deterministic physical model per sampled input
and retains 1..=32 paired outputs. A temperature, stress and pressure obtained
from the same solve do not consume three physical calls. The existing scalar
executor owns input-law admission, Philox addressing, interruption and counts;
its original scalar summaries supply every marginal without a second numerical
implementation. Independent Gaussian/uniform and joint Gaussian input laws,
including singular covariance, retain their original sampling semantics.

The input `UqPlan` names the study and declares its input law and total call cap.
Leave its scalar `compliance_threshold` unset. Each `QoiDefinition` separately
names an output, its physical unit and optional inclusive upper limit. Output
order and model meaning must stay fixed across every callback and resumption.

```rust
use fs_uq::{CorrelationModel, MultiQoiExecution, ParameterUncertainty,
    PropagationMethod, QoiDefinition, UqPlan, UqStatus};
let plan = UqPlan::new("thermal operating study", PropagationMethod::MonteCarlo, 1024)
    .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
    .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
    .with_correlation(CorrelationModel::Independent);
let outputs = [
    QoiDefinition::new("junction temperature", "K").with_upper_limit(311.0),
    QoiDefinition::new("power", "W").with_upper_limit(5.8),
];
let mut run = MultiQoiExecution::new(&plan, &outputs)?;
let report = run.advance(1024, || false, |x| {
    Ok::<_, &'static str>(vec![300.0 + x[0] * x[1], x[0]])
});
assert_eq!(report.status, UqStatus::Complete);
assert_eq!(report.samples_evaluated, 1024);
assert_eq!(report.outputs.len(), 2);
let joint = report.joint_compliance.ok_or("no constrained observations")?;
assert!((0.0..=1.0).contains(&joint.mean));
# Ok::<(), Box<dyn std::error::Error>>(())
```

This algebraic thermal example illustrates the interface, not qualification of a
cooling solver. Units are retained labels; dimensional and model validity remain
caller obligations. Joint compliance counts rows satisfying ALL declared limits,
not the product of marginal compliance frequencies. Outputs without a limit are
not constraints. No limits means no joint event or probability is reported.

A whole response is validated before any column is accepted. Wrong vector length,
nonfinite components and physical errors permanently refuse the whole experiment.
No per-column missing data, skipped draws or replacement observations are allowed.
After refusal, paid counts and raw paired data remain accessible, but summaries
are hidden. A completed marginal's numerical-statistics refusal also refuses the
joint report, consistently with the scalar owner.

Chunking, cancellation, interrupted callbacks and cloning retain completed rows.
`Ok(None)` denotes unfinished work and retries the identical physical input; it
must never censor an inconvenient response. Zero allowance and terminal advances
are no-ops. Callback panics propagate; long physical solvers own cancellation.
Observation storage is reserved before callbacks and capped at 4,194,304 scalar
entries across outputs, within the existing one-million-physical-call cap.
Reports use bounded synchronous scalar sorting/reductions and temporary columns;
these caps do not bound callback-owned solver memory or time.

Marginal and joint standard errors are descriptive, not confidence intervals or
valid stopping criteria. Zero observed joint failures does not prove zero risk.
The study-level evidence dispersion is infinite because output units can differ.
No durable checkpoint, multivariate quantile, automatic optimization integration
or physical-model validation is supplied. The focused Rust regressions are in
`tests/product_multi.rs` and require native execution before claiming a test pass.

## Paired dependence and uncertainty in estimated means

`assess_covariance(max_products, cancelled)` uses the original paired rows to
report the sample covariance, the estimated covariance of the vector mean, and
Pearson correlations. The first matrix divides centered cross products by n-1;
the second divides by n*(n-1). These have different meanings and are not
interchangeable. Matrix indices retain output definitions; covariance (i,j)
has units unit_i*unit_j, whereas correlation is dimensionless.

The implementation independently centers and scales each output, using compensated
reductions. Small differences on a large offset are retained by subtracting an
anchor before scaling. Opposite finite extremes instead scale first. Mixed-unit
rescaling avoids premature products of huge/tiny scales. An unrepresentable
nonzero covariance is None, not zero. Mean covariance can remain available when
sample covariance overflows; correlations can remain available when dimensional
moments overflow or underflow. Constant outputs have zero covariance but undefined
correlations. Subnormal results retain their ordinary f64 quantization limits.

The n*d*(d+1)/2 scalar-product count is checked against both the caller allowance
and a hard 67,108,864 cap before scratch allocation. Normalization uses O(n*d)
scratch within the existing observation cap. Cancellation is checked throughout
normalization and pair reductions; no partial matrix is published. Assessments
are read-only. Fewer than two rows returns None, not a zero matrix. Refused
executions cannot publish covariance. Finite rounding can perturb positive
semidefiniteness; no PSD certificate, repair or distribution fitting is claimed.

## Confidence for joint physical compliance

`assess_joint_compliance(alpha, target_half_width, cancelled)` sends each paired
ALL-limits indicator to the existing `fs-eproc` Bernoulli-mixture confidence
sequence. It does not assume outputs are independent or reconstruct a Gaussian
law from their covariance. With the fixed-event sampling assumptions of that
owner, this is one time-uniform event assessment. Multiple response constraints
inside the fixed conjunction do not require multiplying marginal probabilities.
Searching across different events, limits, models, seeds or confidence levels
still requires separate selection/multiplicity control.

The event and alpha must be chosen before looking at results. Bounds can be
asymmetric; `converged` requires both endpoint distances from the empirical mean
to meet the requested half-width. No observed failures does not return a
probability-one certificate. Cancellation brackets the bounded interval solve
and every row; raw observations, work counts and execution status are unchanged.
The fixed-budget executor does not automatically stop on this assessment.

These are sampling statements for the declared model, not outward numerical
certificates, population transfer or physical qualification. Neither covariance
nor confidence invokes a physical callback. `tests/product_multi_assessment.rs`
adds paired affine/offset/extreme-scale checks, reconstruction through the scalar
Bernoulli owner, boundary decisions and every cancellation checkpoint. Native
execution is required before claiming these Rust tests pass.
