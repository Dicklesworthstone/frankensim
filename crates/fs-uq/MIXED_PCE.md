# Polynomial chaos for uniform tolerances and Gaussian scatter

`PceExecution` fits one total-degree product polynomial to independent Gaussian
and uniform physical inputs, or to an explicitly joint Gaussian law. Gaussian
germs use normalized probabilists' Hermite polynomials. Uniform germs live in
`[0,1]` and use `sqrt(2k+1) P_k(2u-1)`, the normalized shifted Legendre basis.
Epistemic intervals are still not probability distributions. Unknown multivariate
dependence, generic correlation matrices and non-Gaussian correlated laws refuse
before physical work; this API does not infer a copula.

```rust
use fs_uq::{CorrelationModel, ParameterUncertainty, PceConfig, PceExecution,
    PropagationMethod, UqPlan, UqStatus};

let plan = UqPlan::new("junction temperature", PropagationMethod::PolynomialChaos, 128)
    .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
    .with_parameter(ParameterUncertainty::gaussian("resistance", 2.0, 0.1, "K/W"))
    .with_correlation(CorrelationModel::Independent)
    .with_compliance_threshold(310.0);
let mut run = PceExecution::new(&plan, PceConfig::new(2, 96))?;
let report = run.advance(128, || false, |x| Ok::<_, &'static str>(300.0 + x[0] * x[1]));
assert_eq!(report.status, UqStatus::Complete);
let model = run.model().ok_or("fit did not complete")?;
let (temperature, germ_gradient) = model.try_eval_with_gradient(&[0.5, 0.0])?;
assert!((temperature - 310.0).abs() < 1e-9);
// Physical inputs at these germs are power=5 W and resistance=2 K/W.
assert!((germ_gradient[0] - 4.0).abs() < 1e-9); // dT/du; d(power)/du = 2 W.
assert!((germ_gradient[1] - 0.5).abs() < 1e-9); // dT/dz; d(resistance)/dz = 0.1 K/W.
# Ok::<(), Box<dyn std::error::Error>>(())
```

The illustrative thermal response has exact mean 310 K and variance
`4/3 + 25*0.01 + 0.01/3` K^2 under the stated laws. It is not a qualification of
a cooling solver or an assertion that Gaussian resistance is physically suitable
for every operating regime.

## Shared execution and explicit model semantics

Sampling, law admission, scaled QR fitting, conditioning checks, held-out
residuals, compliance indicators, and interruption handling retain their
existing owners. Standardizing a uniform marginal consumes the same draw as
ordinary product Monte Carlo. Its physical value uses the same convex endpoint
interpolation, avoiding `hi-lo` overflow. Singleton intervals remain exact while
still consuming their draw. Pure-Gaussian sequences, design arithmetic, and
successful evaluation/gradient order remain unchanged.

`PceExecution::model()` now returns `&ProductPceModel`, not `&PceModel`. Its
immutable `basis()`, `indices()` and `coefficients()` accessors retain the
reference laws alongside the polynomial. Callers with an explicit old return
type annotation should use `ProductPceModel`; callers reading public coefficient
fields should use the corresponding accessors. The standalone legacy
`PceModel` and `fit_pce` remain Hermite-only and unchanged.

The new model provides scalar evaluation, analytic germ gradients, coefficient
moments, and Sobol decomposition under the retained independent germ measures.
Fallible `try_eval` and `try_eval_with_gradient` refuse invalid coordinates and
nonfinite arithmetic. Uniform germs outside `[0,1]` refuse instead of silently
extrapolating beyond their declared law. Derivatives use the differentiated
Legendre recurrence, so endpoints and polynomial zeros need no singular division.
The scalar convenience methods panic on those same refusals, as documented.
Physical derivatives additionally require the germ-map chain rule. Under a
correlated Gaussian map, germ Sobol indices are not physical-coordinate indices.

The fitting prefix is never contaminated by held-out observations. Only held-out
physical indicators enter compliance frequency; surrogate coefficients do not
certify threshold decisions. Cancellation or an unfinished callback retains
completed work, and a pause at fitting never repeats already paid physical calls.

## Limits and evidence

The existing bounds remain: 1..=16 germs, degree <=8, <=128 total-degree terms,
<=65536 physical calls, at least twice as many training samples as basis terms,
and at least two disjoint held-out calls. Dense regression retains its existing
work and numerical-conditioning admission. Immutable standalone product models
use the same dimension/term/degree limits and validate finite coefficients,
unique terms and a first constant term before evaluation.

Moments and sensitivities describe the fitted polynomial, not its agreement with
physics. Training/validation residuals are measured diagnostics, not error bounds
or confidence sequences. Zero-width physical uncertainty can leave tiny numerical
nonconstant coefficients at positive fit degree; these are not silently clipped.
No basis adaptation, durable checkpoint, correlated-uniform copula, confidence
bound, or automatic physical-model validation is added.

`tests/product_pce_model.rs` checks the polynomial algebra and Gaussian identity.
`tests/product_pce_mixed.rs` checks physical moments, input-draw identity,
continuation, held-out separation and refusal boundaries. Native execution of
these tests is required before claiming a production Rust test pass.
