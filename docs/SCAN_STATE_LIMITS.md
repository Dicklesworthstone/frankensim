# Constrain the simulated temperature and enthalpy history

`ConstrainedScanStudy::new_with_state_bounds` adds lower/upper limits on actual
thermal responses to the existing power, geometry and absorbed-energy limits.
A limit names a supplied endpoint, a mesh vertex and either temperature in
kelvin or specific enthalpy in J/kg. It is not a least-squares fitting target.
The same physical model may fit an observation at one location while limiting
heating somewhere else, including at an earlier time.

## Select the responses

```rust
use fs_ascent::scan_design::constraints::{
    ConstrainedScanStudy, ScanResponse, ScanResponseBudget, ScanStateBound,
};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::ScanObservable;

let limits = vec![ScanStateBound {
    response: ScanResponse {
        endpoint: 1,
        vertex: 0,
        observable: ScanObservable::SpecificEnthalpy,
    },
    lower: None,
    upper: Some(155.0),
    scale: 10.0,
}];
let constrained = ConstrainedScanStudy::new_with_state_bounds(
    cx, study, control_bounds, energy_caps, limits,
    ScanResponseBudget {
        max_responses: 1,
        max_gradient_entries: 128,
        max_transpose_solves: 4,
    },
    32,
)?;
```

Here `study` is the original `ScanPowerStudy`; its endpoint grid, source factory,
materials and solver controls are unchanged. The illustrated transpose allowance
fits a two-endpoint objective and one response at endpoint index 1. Each response
adds `endpoint + 1` reverse solves; a lower/upper pair shares one response sweep.
The gradient-entry allowance covers responses times (all authored source controls
plus mesh vertices). It is a scalar shape cap, not a byte or global memory bound.

Each candidate uses one original nonlinear forward trajectory. The fitting
objective and the independently seeded limits reuse its checked endpoint tape.
Response gradients include earlier heating, source normalization, the actual
caloric chart and optional radiation feedback. They do not differentiate the
squared observation residual: matching a target must not erase a limit's gradient.
On a melting plateau, temperature derivatives can vanish while enthalpy derivatives
still constrain absorbed heat.

## Run, inspect and replay

Use the existing `advance_optimizer`, evaluation-budget extension and `sources`
methods. `evaluate_at` and `replay_with_constraints` return both the physical
objective and complete constraint values/Jacobians, consuming one trajectory
reservation plus the admitted reverse work. They do not replace the accepted
SQP checkpoint. With state limits attached, the cheap `constraints_at` method
refuses rather than silently omitting the physics-dependent rows.

Reports order input-control bounds first, then energy caps, then state bounds;
each bound contributes its lower side before its upper side. Scaled residuals
are feasible at or below zero. A feasibility tolerance `tau` allows `tau*scale`
physical violation, separately from numerical error in the thermal solution.

## Executable synthetic example

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --state-limits
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --state-limits 151.5 156
cargo test -p fs-conduction --test scan_response_rows
cargo test --manifest-path examples/scan_power_fit/Cargo.toml state_limits
```

A unit-mass body starts at 150 J/kg. Fitting 152 and 157 J/kg at one and two
seconds would use 2 then 5 W. The default output limits of 151 and 155 J/kg instead
have the analytical optimum 1 then 4 W, with residual objective 0.025. Both remain
at the 310 K latent plateau. These are synthetic same-model cases, not calibrated
material/process predictions. The executable prints its actual solved powers,
responses, KKT checks and replay diagnostics, and exits unsuccessfully when its
requested convergence or feasibility conditions are not met.

Limits cover only the specified vertices and times. They do not bound unsampled
peaks, ensure continuum accuracy, or provide a machine-control safety system.
SQP may evaluate infeasible intermediate designs. Original physical/material
refusals still propagate. All derivatives hold the grid, geometry, phase charts
and event sides fixed. Native compilation and Rust execution were unavailable
in the delivery environment; independent equation checks are not native test
results.
