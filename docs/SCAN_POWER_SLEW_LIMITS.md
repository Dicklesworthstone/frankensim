# Continuous absorbed-power rise and fall limits

`ScanPowerSlewBound` limits the rate of a heater's **authored absorbed-power
schedule**, in W/s. Power ceilings and energy caps alone cannot enforce this.
It is an input/design constraint, not a modification of the thermal equations
or a hardware-control loop.

For every original linear segment from `(t0, P0)` to `(t1, P1)` the constraint
owner evaluates

```text
rate = (P1 - P0) / (t1 - t0)
rise row = (rate - maximum_rise_w_s) / scale_w_s
fall row = (-rate - maximum_fall_w_s) / scale_w_s
```

Both rows must be nonpositive. A declared feasibility tolerance `tau` permits
`tau * scale_w_s` W/s violation, independently of thermal solver error.
Rise/fall limits are finite and nonnegative; the scale is finite and positive.

## Supported schedules and scope

Constant schedules have no nonzero-rate segments. Ramps and linear tables
produce one rate per original segment, including short segments skipped by a
coarser thermal time grid and the final ramp to a fixed zero shutdown value.
Clamped tails are constant. At a knot the two one-sided slopes are covered;
the schedule need not have a classical derivative exactly at the knot.

**Held tables refuse.** Dividing a held power jump by the time to the next knot
does not bound its actual instantaneous change. Smooth function objects also
refuse; the implementation does not infer an extremum bound from samples.
The method requires finite representable segment durations, rates and reciprocal
durations. It covers declared absorbed power, not an electrical actuator's lag,
incident optical power, scan velocity, or continuous thermal maxima.

Every value in a named source participates, whether or not it is selected as an
optimization variable. A wholly fixed heater can therefore make a candidate
infeasible. The optimizer never clips, smooths, retimes or changes interpolation
to obtain a feasible result. Declare a continuous schedule explicitly.

## Library entry points

`MovingHeatSource::absorbed_power_slew_segments(cx, max_segments)` exposes
original times, power-value indices, endpoint watts, signed rates, and the two
endpoint derivatives. No mesh or thermal solve is needed.

The optimization paths share the same compiled physical rows:

```rust
use fs_ascent::scan_design::constraints::ScanPowerSlewBound;

let rates = vec![ScanPowerSlewBound {
    source: 0,
    maximum_rise_w_s: 3.0,
    maximum_fall_w_s: 2.0,
    scale_w_s: 1.0,
}];
```

Pass these to either constructor:

```text
ConstrainedScanStudy::new_with_power_slew_bounds(
    cx, study, control_bounds, energy_caps, state_bounds, peak_bounds,
    response_budget, max_peak_samples, rates, max_slew_segments,
    max_kkt_dimension)

ValueConstrainedScanStudy::new_with_power_slew_bounds(
    cx, study, control_bounds, energy_caps, state_bounds, peak_bounds,
    rates, max_slew_segments, value_constraint_budget)
```

The value-only owner can be passed directly to `PopulationScanStudy`; its
existing serial/parallel callbacks and final replay retain the new rows.
It still needs no derivative tape and accepts valid physical phase corners.
SQP receives derivatives through the actual logarithmic or nonnegative-linear
power coordinates. Unselected powers contribute to values, not gradients.
Geometry variables do not acquire spurious power-rate derivatives.

`max_slew_segments` limits the combined segment count across all declared
sources. Each segment adds two scalar constraints, counted against the original
value-row or dense-KKT admission. Duplicate source declarations refuse.
Old constructors attach no rate limits and retain their prior row ordering.

## Reports and focused example

Both physical constraint reports expose `power_slew_values_w_s` in declaration
order, then chronological segment order. Signed inequality ordering is:
nonnegative-domain rows, control bounds, energy caps, rate rise/fall pairs,
point-state bounds, and peak bounds.

The `g1_population_fit_enforces_the_fixed_shutdown_rate_not_just_power_or_energy`
regression uses the actual P1 enthalpy model with power values `[0, P, 0]` at
times `[0, 1, 2]` seconds. Only `P` is optimized. Under the declared two-step
backward-Euler sampling, its two enthalpy predictions are both `150 + P` J/kg.
The targets are both 160 J/kg with scale 10 J/kg. Without rate constraints the
discrete objective is minimized at 10 W. A 3 W/s rise limit and 2 W/s shutdown
limit require `P <= 2 W`; the analytical constrained optimum has objective 0.64.
That feasible result misses the targets. It is not labeled target recovery.
These are equations of the stated discrete fixture, not a time-continuum
accuracy claim for the coarse backward-Euler trajectory.

Native checks for this change:

```bash
cargo test -p fs-conduction --test moving_heat_slew
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --test spatial_design slew
```

The new tests cover source rates and derivatives, fixed heaters and shutdown,
coordinate Jacobians, shared-batch agreement, admission, cancellation and the
end-to-end constrained population fit. Their execution status belongs in the
delivery report; this document makes no native-test, performance, physical
calibration or global-optimality claim.
