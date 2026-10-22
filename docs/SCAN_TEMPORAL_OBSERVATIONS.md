# Time-averaged thermal measurements

`ScanTemporalObservation<O>` adds finite-exposure means to the existing scan
objective. `O` may be a vertex observation, a mesh-bound spatial readout, or a
loss-wrapped observation. Temperature and liquid mass fraction are resolved
from each node's actual material chart before spatial and temporal averaging.

Use `instantaneous(row)` for an unchanged point-in-time reading. For an
exposure, `trapezoidal_mean(cx, row, &endpoints_s, first_endpoint, max_samples)`
averages through `row.endpoint`, inclusive. The selected window must contain
at least two physical endpoints. Its target and residual scale retain kelvin,
J/kg or dimensionless mass-fraction units. A time mean is not an integral;
multiply a reported mean by its exposure duration when an integrated readout
is needed, with the corresponding change of target/scale units outside this API.

The window begins at an accepted endpoint, not the initial state of the run.
All interior endpoints contribute with nonuniform trapezoidal weights. The
exact clock is retained: modifying or refining the selected grid requires
rebuilding the observation. This integrates a piecewise-linear reconstruction
of the simulated readout; it neither changes the thermal stepper nor certifies
continuum quadrature accuracy. A single-sample exposure, nested time average,
nonfinite duration or underflowed positive coefficient is refused.

The fitting objective forms ONE residual after the whole time/spatial mean.
Huber loss is also applied once, not independently to the exposure's frames.
The loss wrapper can be inside or outside the temporal wrapper. Instantaneous
and averaged readings can coexist in a homogeneous observation vector.

The same observation type works in `ScanPowerStudy`, L-BFGS, SQP, value-only
population evaluations and immutable shared batches. It uses the original
trajectory and adjoints. Each endpoint receives its weighted portion of the
measurement seed; downstream enthalpy sensitivities still traverse all earlier
steps. Peak/state constraints remain independent physical requirements: a low
exposure-average loss does not permit an instantaneous hot spot.

`max_observations` in `ScanObjectiveBudget` and `ScanValueBudget` bounds the
SUM of endpoint contributions, not merely the observation-vector length. For
example, one instantaneous reading and one four-endpoint mean need allowance
five. Each mean also has its own constructor sample cap. Value-only evaluation
retains one scalar accumulator per row, not a trajectory of thermal fields.
Checkpointed differentiation retains one raw readout per contribution, bounded
by the same allowance, and checks every sample bit-for-bit on replay. Matching
averages cannot hide changed intermediate samples.

## Example interpretation

For an adiabatic uniform-mass 1 kg specimen initially at 150 J/kg, suppose
mean enthalpy is 152 J/kg at t=1 s and its time mean over [1,2] s is 154.5 J/kg.
For the two-step latent-plateau fixture, these imply powers 2 W and 5 W.
Misreading the second number as the state at t=2 s instead gives 2 W and 2.5 W.
Those are discrete-fixture identities, not experimental calibration claims.

## Verification scope

G1/G3 tests cover nonuniform quadrature, one-residual semantics, Huber wrapper
order, original instantaneous bits, physical/control gradients, full/checkpoint/
value agreement, constrained-row invariance and inverse recovery with a negative
control. G0/G4 tests cover clock/sample-budget refusals and cancellation.
Native Rust compilation/execution was unavailable in the authoring environment;
independent NumPy P1/SciPy checks do not replace those native tests.

## CSV-driven fitting

The existing `spatial_measurement_fit` executable accepts `start..end` in the
`time_s` column to request one exposure mean. A single number keeps the original
instantaneous interpretation; the six-column header and spatial selection syntax
are unchanged. For example:

```csv
time_s,observable,target,scale,readout,selection
1,specific_enthalpy_j_kg,152,10,volume,0;1;2;3;4;5
1..2,specific_enthalpy_j_kg,154.5,10,volume,0;1;2;3;4;5
1..2,liquid_mass_fraction,0.545,0.1,volume,0;1;2;3;4;5
```

Run the shipped fixture with four substeps in every measurement/event interval:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/temporal_measurements.csv 0 0 --substeps 4
```

`--huber positive_threshold` remains available with either interpretation. Flags
may appear before or after positional arguments. `--substeps` defaults to one
and accepts 1 through 128. The grid always includes exact exposure endpoints,
instantaneous measurement times and source switches. Substeps that cannot be
represented distinctly in binary64 are refused, not silently merged.

This executable still uses its declared synthetic specimen and accepts only
`0 < start < end <= 2` seconds. It does not infer a laboratory setup from CSV,
include the initial state in an exposure, or calibrate sensor dynamics. Refining
the fixed grid is explicit; it is not an automatic time-error certificate.
File, measurement-row, physical-endpoint and total measurement-sample caps are
128 KiB, 256, 4096 and 65536 respectively. Overlapping windows each consume their
own contribution count. The constructor and physical solver retain their other
budgets and rejection checks.

Output identifies exposure starts/ends separately from the raw physical
prediction, residual and loss. The optimizer retains nonnegative powers, exact
shutdown, checkpointed adjoints and its reserved final physical replay.
Convergence and target mismatch remain distinct; robust or inconsistent data
can yield a converged fit with nonzero residuals.

Three additional authored G0/G1/G3 regressions exercise exposure-aware fitting
versus the incorrect-endpoint interpretation, exact grid/order/sample accounting,
invalid windows and resource caps. They have not been executed natively in the
authoring environment.
