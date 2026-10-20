# Derivative-free thermal scan fitting

`fs_dfo::neldermead::state::NelderMeadState` can drive
`ScanPowerStudy::evaluate_value` directly. This searches the original nonlinear
thermal model without requiring an objective derivative or an adjoint tape.
Valid phase-chart corners remain usable. The forward solver's material,
energy, spatial-source and cancellation checks are not weakened, and callback
errors are not converted into penalties.

## Run the physical example

From the repository root:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- 0.2 0.4 50 400
cargo test --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit
cargo test -p fs-dfo --lib neldermead
```

The example declares its synthetic model rather than inferring one: a unit cube
with density 1 kg/m3, conductivity 2 W/(m K), an adiabatic exterior and a
100 J/kg latent plateau at 310 K. The initial specific enthalpy is exactly
100 J/kg, at the solidus, and both selected heater powers are initially zero.
The heater moves between two declared positions, then switches off at t=2 s.
Measurements are actual volume means of nodal liquid **mass** fraction. Here
uniform reference density makes that mean a bulk mass fraction as well.

The default targets are 0.1 after one second and 0.3 after two seconds. When
all nodes remain on the plateau, their mean energy balance implies 10 W in the
first interval and 20 W in the second. The alternate command targets 0.2 and
0.4, implying 20/20 W in that regime. These are reference expectations, not
claims of experimental calibration or a global optimum.

The optional fourth positional value is the total physical evaluation budget,
including one reserved replay. The example first spends an eight-callback
search allowance, then explicitly extends both the solver and physical-study
budgets. It keeps the same search state; it does not restart the simplex.
A value evaluation runs the two original thermal endpoints, and the example
has zero retained-adjoint allowances to prohibit a hidden derivative path.
The exported schedule is evaluated again before reporting recovery. Output
contains actual unsnapped powers, endpoint measurements, objective mismatch,
stop reason and work counts. An exhausted or collapsed search that misses the
targets exits unsuccessfully; a small simplex alone is never called recovery.

## Solver surface

Construct `NelderMeadState::new(initial, signed_axis_steps, config)` with a
`NelderMeadConfig` containing a total callback allowance, simplex-shape cap,
optional objective target, two absolute spread tolerances, and optional finite
box bounds. Construction invokes no callback, even with a zero work allowance.
Every projected initial axis step must remain distinct from the starting point.

`advance(cx, additional_evaluations, callback)` accepts a fallible callback:

```rust
let report = simplex.advance(cx, allowance, &mut |point| {
    study.evaluate_value(cx, point).map(|r| r.physical.value)
})?;
```

Keep the study, callback, coordinates and physical inputs fixed across calls.
Use `extend_evaluation_budget` to increase a lifetime allowance. The solver and
study counters are independent: the application must leave enough study budget
for each admitted callback and any final replay. An objective-only evaluation
is not a check of a `ConstrainedScanStudy`'s additional inequalities.

The example uses nonnegative-linear powers with a 1 W reference, so its two
optimizer coordinates numerically equal watts. General solver box bounds are
in **optimizer coordinates**, not automatically SI units. For other power or
geometry scales, transform the bounds explicitly.

## Interruptions and hard work limits

Each callback attempt is reserved before invocation, including failures.
Initialization, expansion, contraction and every shrink vertex count against
the same ceiling. A returned error retains the candidate for a separately
charged retry; a nonfinite score is an error, not an infinite penalty.

A successful callback result is cached before the next cancellation check.
Resume consumes that result without repeating the expensive simulation. An
unfinished shrink retains its evaluated replacements separately and publishes
the new simplex only when every replacement has a finite score. The reported
best point may already come from an unfinished iteration, because it is a real
completed observation. No unevaluated point is reported as best.

`Paused`, `EvaluationLimit`, `TargetReached` and `SimplexSmall` are distinct.
The shape allowance bounds simplex entries, not peak allocated bytes; trial,
centroid, best and staged-shrink arrays require additional O(n^2) storage.
Vector traversals poll cancellation; sorting is bounded by the admitted
simplex dimension, without a hard real-time latency claim.

The original `fs_dfo::nelder_mead` convenience function uses the same worker.
It now respects `max_evals` even when the initial simplex does not fit. Its
callback-only API has no `Cx`, and retains a panic boundary for invalid input
or numerical refusal. Zero budget is a modeling error there; use the stateful
API to represent a zero-work pause. Negative-infinite target still means no
target stopping. Other targets, objectives and generated points must be finite.

## Limits

This is deterministic local Nelder-Mead, not a global optimizer or a KKT
certificate. Optional box projection changes candidate locations before their
actual evaluations and can collapse a simplex at a boundary. It does not add
nonlinear state constraints. No outcome is silently snapped to a boundary.
A derivative-free search can still stagnate on a flat saturated response.
Neither this optimizer nor the example smooths phase corners, changes the time
grid, relaxes physical acceptance, or treats a failed simulation as valid data.
