# Outlier-resistant thermal scan fitting

`ScanObservationWithLoss<O>` attaches a fixed residual-loss policy to the
existing vertex or mesh-bound spatial observation. Use the same observations
in `ScanPowerObjective`, `ScanPowerStudy`, `ConstrainedScanStudy`, or the generic
`evaluate_scan_observations` function. The original quadratic default is unchanged.

```rust
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    ScanObservable, ScanObservation,
    loss::{ScanObservationWithLoss, ScanResidualLoss},
};

let observation = ScanObservationWithLoss::new(
    ScanObservation {
        endpoint: 1, vertex: 0,
        observable: ScanObservable::SpecificEnthalpy,
        target: 157.0, scale: 10.0,
    },
    ScanResidualLoss::Huber { threshold: 1.5 },
)?;
```

Thresholds are dimensionless: this declaration changes from quadratic to linear
loss at an absolute residual of 15 J/kg. A wrapped spatial observation first
forms its actual weighted temperature, enthalpy or liquid-mass-fraction readout.
The loss is applied to that ONE measurement, not separately to its contributing
vertices. Different observations can choose different fixed thresholds or keep
quadratic scoring. An outer wrapper replaces, rather than composes, an inner loss.

For `r = (prediction - target) / scale`, Huber scoring uses
`rho = r*r/2` within the threshold `delta`, and
`rho = delta*(abs(r)-delta/2)` outside it. The derivative with respect to the
physical prediction is `r/scale` or `copysign(delta,r)/scale`. This is the
Huber (1964) and `scipy.special.huber` convention. The shared reverse owner
injects this derivative before the spatial transpose and material-chart chain.
Full-storage and checkpointed runs use the same policy and accumulation order.
No new forward or transpose solve is needed merely to select Huber scoring.

## CSV executable

The CSV format and readout semantics are unchanged. Quadratic fitting remains:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/spatial_measurements_contaminated.csv 0 0
```

Select a standardized threshold explicitly:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/spatial_measurements_contaminated.csv 0 0 --huber 1
```

The fixture is synthetic. Four readings at one second give 152 J/kg; four at
two seconds give 157 J/kg; one deliberately contaminated final reading gives
197 J/kg. All are full-volume means with a 1 J/kg residual scale. In the
analytical latent-plateau reduction, the quadratic optimum is 2/13 W, while
Huber with threshold one gives 2/5.25 W. The uncontaminated schedule is 2/5 W.
Huber reduces the outlier's influence; it does not delete the reading or recover
the clean schedule exactly. The final raw outlier residual remains -39.75 J/kg.
These are analytical reference values, not a claim of native executable results.

The executable reports the selected policy, each unmodified physical prediction,
raw standardized residual, loss contribution and derivative with respect to that
residual. Existing source export, exact shutdown, physical replay, cumulative
budgets and convergence/nonnegative-domain checks remain in force.

## Scope and checks

Hard point, hot-spot, power, geometry and energy constraints do not use the fitting
loss. A lower robust objective cannot hide a constraint violation. Invalid or
nonfinite measurements, policies, mesh identities, physical states or solver
results still refuse; robustness is not a substitute for physical admission.
Thresholds and scales are fixed declarations, not inferred noise estimates.
This is not a correlated likelihood, outlier classifier, uncertainty interval,
sensor-dynamics model, global optimizer or experimental calibration certificate.
Huber is C1 but not C2 at the threshold; no second derivative is exposed.

Focused Rust regressions live in the loss module, the spatial-design integration
target, and the CSV example's tests. Native commands are:

```bash
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --test spatial_design
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --example spatial_measurement_fit
```
