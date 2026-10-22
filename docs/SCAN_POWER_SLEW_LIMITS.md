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

## Executable ramped-heater fitting

The existing `constrained_scan_fit` example accepts an explicit continuous-power
mode. Commands without these flags retain the original held-power specimen.

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example constrained_scan_fit -- --power-slew 50 50 --workers 4

cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example constrained_scan_fit -- --power-slew 50 30 --ramp-steps 32
```

`--power-slew RISE FALL` uses W/s and permits finite nonnegative limits.
`--ramp-steps N` selects 8–128 equal steps per second (default 16) and requires
the continuous-power mode. It changes the declared time discretization, not the
authored rate constraints. Each original knot remains exactly on the grid.

This mode uses the original cold, uniform one-kilogram P1 specimen and a linear
power table `[0, P1, P2, 0]` at `[0, 1, 2, 3]` seconds. Only `P1` and `P2` are
optimized, each in `[0,50]` W. The targets are **10% liquid at 2 s and 30% at
3 s**, unlike the held-mode observations at 1 and 2 s. Both fixed zero values
remain unchanged on export. The full shutdown interval is simulated, contributes
to the energy cap, and participates in the sampled melt-fraction envelope.
Default limits remain 60 J and a 0.5 peak-fraction envelope.

Output identifies the interpolation, requested rate limits and actual signed
segment rates. It emits the four authored `(time_s,power_w)` knots separately
from the two measured fraction predictions. A finite-rate limit is not silently
applied to the old held schedule. A feasible result that misses either target
still returns an unsuccessful exit; exhausting the search is not optimality.

For this uniform fixture only, the exact **discrete** balances for `N` steps per
second are

```text
h(2) = 80 + P1 + (N+1)/(2N) * P2       J/kg
h(3) = 80 + P1 + P2                   J/kg
authored absorbed energy = P1 + P2    J
rates = [P1, P2-P1, -P2]             W/s
```

These equations are reference checks, never substituted for the finite-element
solver. At `N=16`, loose 50/50 W/s limits admit the target schedule approximately
`P1=7.333333 W`, `P2=42.666667 W`. A 30 W/s fall limit instead gives the analytical
best fit `P1=17.03125 W`, `P2=30 W`, predictions `0.1296875` and `0.2703125`, and
objective `0.088134765625`. Omitting shutdown would admit the target schedule
despite a 12.666667 W/s fall-rate violation.

Refining the grid changes the intermediate-time quadrature. The unconstrained
target schedule approaches 10/40 W as `N` increases; the finite-grid numbers above
are not time-continuum accuracy claims. The envelope bounds only sampled states
and adds its explicitly reported smoothing margin.

Checkpoint/resume, parallel evaluation and final physical replay remain available:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example constrained_scan_fit -- --power-slew 50 30 --workers 4 \
  --checkpoint ramps.checkpoint --pause-after 48
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example constrained_scan_fit -- --power-slew 50 30 --workers 8 \
  --resume ramps.checkpoint --checkpoint ramps.checkpoint
```

Recovery binds the executable, energy/peak limits, power mode, both rate limits
and ramp-grid resolution. Changing workers is permitted; changing the physical
rates or grid while reusing completed scores is refused. Stored physical/search
work is not refunded. The new executable regressions check option admission,
the actual shutdown trajectory against these discrete balances, all signed rate
rows, and restart with changed workers and rejected changed physical policies.
