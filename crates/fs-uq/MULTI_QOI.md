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
