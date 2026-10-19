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

## Joint path, footprint and power design

`ScanPowerStudy::new_joint` uses the same forward solver, full-history adjoint,
fallible L-BFGS engine, cumulative evaluation budget and source export as the
power-only constructor. It additionally accepts `Vec<ScanGeometryControl>`.
There is no second optimizer or approximate geometry-to-heat model. Power or
geometry controls may be empty individually, permitting geometry-only fits.
The original `new` entry point remains power-only.

```rust
let geometry = vec![
    ScanGeometryControl::CenterValue {
        source: 0, axis: 0, value: 0, reference_m: 0.3, scale_m: 0.25,
    },
    ScanGeometryControl::CenterValue {
        source: 0, axis: 0, value: 1, reference_m: 0.65, scale_m: 0.25,
    },
    ScanGeometryControl::HalfWidth {
        source: 0, axis: 0, reference_m: 1.4,
    },
];
let mut study = ScanPowerStudy::new_joint(
    cx, sources, power_controls, geometry, objective, factory, evaluation_limit,
)?;
```

The optimizer vector contains selected log powers first, then geometry entries
in their declaration order. Center values use `c = reference_m + scale_m * x`:
negative and zero world coordinates are allowed. Widths use a positive log
coordinate. Accordingly the full-history adjoint is multiplied by `scale_m`
for centers and by the current width for width controls. Geometry is differentiated
through the existing **normalized** source, including the changing normalization
at clipped boundaries. A center or width change therefore redistributes fixed
absorbed power rather than creating additional watts.

`geometry_at` returns selected physical metre values; `powers_at` still returns
only selected powers. Both take the complete optimizer point. `sources_at`
exports reusable source declarations containing both changes without another
simulation. `ScanPowerEvaluation::geometry_m` and `powers_w` retain the two
physical coordinate groups; `gradient` follows the combined optimizer vector.

Path edits address original constant, ramp (`from`, `to`) or table values.
Times, held/linear interpolation, unselected coordinates, power shutdowns,
material laws and reference mass remain fixed. Duplicate coordinates and invalid
scales or widths refuse. Polynomial coefficient editing is not supported by
this constructor; untouched polynomial histories are retained. The reusable
lower-level declaration editor is `MovingHeatSource::with_geometry_values`,
with `MovingHeatGeometryControl` in `moving::sensitivity::geometry`.

### Executable spatial-recovery fixture

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --geometry
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --geometry 0.25 0.75 1.2 2 3
cargo test --manifest-path examples/scan_power_fit/Cargo.toml
```

This explicitly synthetic, same-model fixture generates spatial enthalpy targets
on a 27-node cube at four times. It fits both ends of an x-axis ramp, the x half
width, and two held power levels. The default targets are a path from 0.2 to
0.8 m, width 1.15 m, and powers 2 then 3 W. Initial values are 0.3/0.65 m,
1.4 m, and 1.7/2.6 W. The shutdown at 0.8 s stays zero. CSV output reports the
actual fitted path, width and powers, not the target parameters. The example
reserves a final evaluation for replay and exits unsuccessfully on nonconvergence.

The targets are generated by the forward model, not measured data: recovery is
an integration test, not independent physical validation. Power-only fitting
cannot generally correct a misplaced spatial profile. Multiple spatial and
temporal observations are needed to separate path, width and power effects;
this API does not promise identifiability or a global optimum.

Native source support, optional motion limits, phase/material smoothness and
solver acceptance remain mandatory. A trial that violates them returns the
original error; it is not clipped into the mesh or replaced by a fabricated
penalty. This is unconstrained local design in the declared coordinates, not
box-constrained machine planning or differentiation of event/adaptive times.
