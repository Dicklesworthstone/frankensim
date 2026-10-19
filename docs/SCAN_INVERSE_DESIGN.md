# Fit absorbed-power schedules to complete thermal histories

`evaluate_scan_objective` in
`fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective`
runs the original nonlinear scan solver and reverses the complete fixed-grid
history. A `ScanObservation` identifies an endpoint, a mesh vertex, temperature
or specific enthalpy, a target and a positive scale in the same physical units.
The objective is one half the sum of squared scaled residuals. Scales are not
implicitly interpreted as measurement standard deviations. Duplicate observations
add independent terms.

The result includes predictions in the original observation order, normalized
residuals, the original last physical report, derivatives of initial nodal
specific enthalpy, and every source's authored-control gradients. Source controls
retain their original constant/ramp/table/Chebyshev representation. Intermediate
observations inject their own reverse seeds; downstream history sensitivities
are never replaced by those local seeds. A latent plateau has zero temperature
slope, but its enthalpy target still responds to deposited heat.

`EventLeft` requires an explicit grid that already contains all power, path,
boundary, reservoir and selected motion stops. An omitted stop is an error,
not permission to change the design discretization. `RightEndpoint` retains the
original explicit-endpoint convention. Neither differentiates event times,
adaptive decisions or motion-limit selection. Source position, width and power
values can be differentiated on this fixed grid. Mesh, phase/material charts,
reference states and boundary histories remain fixed; classical derivative
admission at material/chart corners is unchanged.

## Existing optimizer, actual physical objective

Enable `fs-ascent`'s existing `conduction-design` feature and use
`fs_ascent::scan_design::{ScanPowerStudy, ScanPowerControl, ScanPowerObjective}`.
Each control identifies a source and an authored power-value index, with an
explicit positive watt reference. Its coordinate is
`x = ln(P / reference_w)`. Returned coordinate gradients multiply the complete
physical per-watt gradient by P. Power never becomes negative by clipping a
trial or zeroing a gradient. Selected initial values must be positive; leave
exact zero/off values unselected so they remain exactly zero.

This positive-value parameterization supports constants, ramp endpoints and
table values. It refuses Chebyshev power coefficients: making coefficients
positive does not establish a nonnegative polynomial. The lower-level full
objective still returns polynomial coefficient derivatives for callers with a
separate feasible parameterization.

The factory is only a physical-run constructor. It receives the rebuilt source
declarations and must use them unchanged, starting from the same original mesh,
materials, thermal history, boundaries and optional radiation each time. The
study checks the returned source declarations for equality. The factory must
also retain any selected motion limit. Do not capture mutable non-source model
parameters that change between probes.

`begin_optimizer(cx, memory)` initializes the existing fallible L-BFGS engine.
`advance_optimizer(cx, &mut state, &stop_rule, max_iterations)` continues its
checkpoint. The native Wolfe search, stopping algebra, and error handling are
unchanged. The study separately counts reserved physical evaluations, including
failed probes; optimizer cloning does not refund them. `extend_evaluation_budget`
raises the total allowance without changing the physical problem. Original
source/material/phase/solver errors propagate instead of becoming invented
penalties. Supply this study's own optimizer checkpoint on continuation.

`sources_at(cx, &state.x)` exports the complete fitted physical sources without
another solve or evaluation reservation. This uses the same source rebuilder
as optimizer probes, preserves all unselected shutdown entries, and permits a
subsequent application run to reuse the calibrated declarations directly.
The study admits the sum of all source schedule controls before copying them;
its factory must return an unadvanced run as well as the exact supplied sources.

## Executable inverse fixture

```sh
cargo run --manifest-path examples/scan_power_fit/Cargo.toml
cargo test --manifest-path examples/scan_power_fit/Cargo.toml
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- 153 159 2 2
```

The standalone Cargo example is a consumer of the ordinary workspace libraries;
it does not modify their dependency graph or copy a solver. On a unit-volume,
unit-reference-density synthetic body, two one-second held power levels must
produce 152 and 157 J/kg from 150 J/kg. Both targets lie inside a 310 K latent
plateau. Their unique uniform-input solution is 2 W followed by 5 W. The example
starts both levels at 1 W, stops at an explicit evaluation budget, and resumes.
Residuals use an explicit 10 J/kg scale. Of its 128 physical evaluation slots,
one is reserved for replay of the accepted fitted schedule before output.
The command prints diagnostic comment lines with predictions, loss, gradient
norm, actual evaluation count and final nodal temperatures, followed by a CSV
held schedule (`time_s,absorbed_power_w`). The shutdown sample stays exactly zero.
It exits unsuccessfully if the optimizer does not reach its gradient threshold;
an iteration cap, stall or exhausted budget is not relabeled convergence.

The optional first two arguments are the requested enthalpies after one and two
seconds. The next two are starting powers, not answers. The alternate command
above requests 153 and 159 J/kg from 2 W initial levels, corresponding to the
independent uniform-input solution 3 W then 6 W. Targets must increase strictly
inside (150,200) J/kg; initial positive powers must sum to less than 50 W.
These restrictions describe this latent-plateau fixture, not the general study
API. `--help` lists units and arguments. The executable is a numerical inverse
fixture, not calibration of a real heater, laser or material.

Full-storage objective evaluation has explicit endpoint, observation,
source-node, schedule-control and radiation-feedback caps. These are shape/work
allowances, not peak byte measurements. Retained sparse tangents follow the
admitted mesh/contact topology. Long trajectories need a separate checkpoint
strategy; this interface does not silently grow an unbounded tape or claim a
global optimum, unique real-world parameter identification, a continuum error
bound, or an experimental prediction certificate.
