# Physical scan-path speed limits

`ScanMotionLimit` controls numerical timestep resolution relative to a heating
footprint. It does not restrict the speed of an optimized source. The separate
`ScanPathSpeedBound` constrains the physical three-dimensional center velocity
in metres per second, including heater-off travel and fixed return positions.

## Continuous paths and asynchronous clocks

`MovingHeatSource::path_velocity_segments` accepts constant, ramp and
linearly-interpolated position schedules. All three axes may have different
timestamps. The union of their clocks partitions the path into intervals with
one constant velocity vector. An active axis keeps its *own original segment
duration* and authored endpoint indices, even when another axis subdivides it.
Clamped tails have zero velocity. There is no midpoint sampling or thermal-grid
sampling in this calculation.

Held position tables and Chebyshev profiles are unsupported here and return an
error. A chord through a held jump is not a finite velocity; a chord through a
smooth returning path is not its maximum velocity. The existing timestep-motion
limiter continues to support its previous signal types independently.

`MovingHeatPathBudget` bounds the number of original timestamps and the number
of merged intervals. The optimizer's `ScanPathSpeedBudget` bounds timestamps per
named source and the combined retained intervals across all named sources.

## Use in design

Both `ConstrainedScanStudy::new_with_actuator_bounds` and
`ValueConstrainedScanStudy::new_with_actuator_bounds` accept `ScanActuatorBounds`.
This declaration can contain path-speed bounds, the existing power rise/fall
bounds, or both. Existing constructors retain their previous behavior.

```rust
use fs_ascent::scan_design::constraints::{
    ScanActuatorBounds, ScanPathSpeedBound, ScanPathSpeedBudget,
};

let actuators = ScanActuatorBounds {
    path_speed: vec![ScanPathSpeedBound {
        source: 0,
        maximum_m_s: 0.25,
        scale_m_s: 0.5,
    }],
    path_budget: ScanPathSpeedBudget {
        max_knots_per_source: 64,
        max_segments: 128,
    },
    ..Default::default()
};
```

The fitting objective and all thermal requirements retain their original
forward solver and adjoints. Each path interval adds one inequality,

    g = (vx/scale)^2 + (vy/scale)^2 + (vz/scale)^2 - (maximum/scale)^2 <= 0.

This is a coupled Euclidean-speed constraint, not three independent coordinate
bounds. It is differentiable at rest. A scaled feasibility tolerance `tau`
permits `speed^2 <= maximum^2 + tau*scale^2`; it does **not** permit a simple
`tau*scale` additive speed error. A zero cap is admitted, but its squared
constraint is degenerate for first-order KKT checks. Unselect stationary path
coordinates instead when practical.

All original fixed endpoints participate. Only selected center values receive
Jacobians, through their actual metre-per-optimizer-coordinate scales. Power
and footprint-width controls have zero path-speed derivatives. Clock times are
fixed. Values and derivatives use the same compiled segment rows.

Reports expose `path_speed_values`: the source, start/end times, velocity
vector and Euclidean speed. Row ordering is nonnegative-power domain rows,
control bounds, energy caps, power slew, path speed, point-state bounds, then
peak bounds. Both row and dense KKT admission include the new intervals.

Value-only serial/batch evaluations pass all these rows to existing population
selection. `begin_sqp_at` preserves the complete actuator declaration when
switching to adjoint polishing, including power slew. It must not rebuild only
thermal limits and silently omit actuator requirements. Replay retains the
original constraint vector, and existing physical-evaluation accounting applies.

## Why the vector and full path matter

A simultaneous `vx=vy=0.4 m/s` move travels at approximately `0.565685 m/s`.
Separate coordinate checks against `0.5 m/s` would incorrectly admit it.
Likewise, returning to the starting point does not imply zero travel speed.
Neither turning off the heater nor widening its footprint removes a speed limit.

The integration fixture optimizes an intermediate x position with fixed
`x(0)=x(2)=0.25 m`. Heating measured at one second prefers `x(1)=0.75 m`, but a
`0.25 m/s` cap admits at most `x(1)=0.5 m`. A feasible best fit therefore retains
measurement mismatch. The regression executes the original P1 enthalpy model
and checks replay of the constrained result, not an after-the-fact clipped path.

Focused checks:

```sh
cargo test -p fs-conduction --test moving_heat_velocity
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --test spatial_design path_speed
```

These are piecewise-linear speed constraints evaluated in ordinary binary64,
not outward-rounded certificates. They do not bound acceleration, jerk,
collisions, actuator tracking error, thermal discretization error, or unsampled
thermal peaks. In particular, bounded speed does not remove velocity jumps at
piecewise-linear corners. No global-optimum or physical-machine safety claim is
attached to a successful simulation or a feasible optimizer report.
