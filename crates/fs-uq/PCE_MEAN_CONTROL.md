# Polynomial surrogates as physical-mean controls

`PceMeanControl` freezes a completed `PceExecution` for a separate, initially
empty `UqExecution`. Both must name the same physical parameters, units, laws,
dependence and QoI. Production must have a different Philox seed. The original
production plan, including its count and compliance threshold, remains immutable.

The estimator averages `Y - weight * (P(germ) - E[P])`. Here `Y` is an actual
production model response, not a surrogate prediction. Under independent
production draws from the declared ideal law, the correction has expectation
zero for any fixed finite weight and any fixed admitted polynomial. An inaccurate
predictor can increase variance; it does not substitute its biased mean for the
physical-model expectation. Independence from production outcomes remains a
caller obligation; a different seed is not proof against data-dependent selection.
Finite-generator grids, numerical rounding and solver/model errors are not bounded.

Gaussian, independent mixed Gaussian/uniform, and explicitly joint Gaussian
fits are supported, including singular correlation and deterministic marginals.
The original independent germs are replayed using the canonical sampler and
production ordinals. No inversion of a singular physical transformation is used.
The constant coefficient is removed before evaluation rather than subtracting
two potentially large predictions. The caller chooses the weight before sampling;
zero retains ordinary Monte Carlo and a variance-increasing weight is not clipped.

`assess` performs no physical callbacks and changes neither raw response bits
nor compliance, quantile or risk summaries. It returns a mean-only comparison,
including raw and adjusted descriptive standard errors and their variance ratio.
Adjusted observations must never be interpreted as the physical QoI distribution.
A singleton has no standard error; zero raw variance has no variance ratio.
Failures do not publish a selected-prefix estimate. Partial estimates are
only diagnostics, not confidence sequences or optional-stopping guarantees.

`assess_interruptible` polls before every regenerated sample and around the
shared bounded moment reductions. Cancellation returns no partial assessment;
the immutable control and raw execution remain reusable. At construction, the
full production allowance must fit at most 67,108,864 polynomial term-coordinate
products. No new dependency, RNG, regression solver or statistical reduction is
introduced. There is no checkpoint codec for the combined control in this API.

```rust
use fs_uq::{CorrelationModel, ParameterUncertainty, PceConfig, PceExecution,
    PceMeanControl, PropagationMethod, UqExecution, UqPlan, UqStatus};
let plan = UqPlan::new("temperature", PropagationMethod::PolynomialChaos, 128)
    .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
    .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
    .with_correlation(CorrelationModel::Independent);
let physics = |x: &[f64]| Ok::<_, &'static str>(300.0 + x[0] * x[1]);
let mut fit = PceExecution::new(&plan, PceConfig::new(2, 96))?;
assert_eq!(fit.advance(128, || false, physics).status, UqStatus::Complete);
let mut production_plan = plan.clone();
production_plan.method = PropagationMethod::MonteCarlo;
production_plan.seed ^= 0x5043_454d_4541_4e31;
production_plan.budget_max_samples = 256;
let mut production = UqExecution::new(&production_plan)?;
let control = PceMeanControl::new(&production, &fit, 1.0)?;
production.advance(256, || false, physics);
let estimate = control.assess(&production)?.ok_or("no production observations")?;
assert!((estimate.mean - 310.0).abs() < 1e-9);
# Ok::<(), Box<dyn std::error::Error>>(())
```

This algebraic thermal fixture illustrates the API, not qualification of a
cooling solver. `tests/product_pce_control.rs` exercises exact and imperfect
controls, variance increases, singular maps, binding refusals, cancellation,
chunk equivalence and unchanged raw physics. Native execution remains required
before claiming these Rust regressions pass.
