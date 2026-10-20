# Heater activation and exact off states

`ScanPowerStudy::new_joint_with_power_coordinates` selects either the original
`ScanPowerCoordinates::PositiveLog` or `ScanPowerCoordinates::NonnegativeLinear`.
The former remains the default of `new` and `new_joint`.

The linear option uses `P = reference_w * x`, with `x >= 0` and an explicit,
positive watt scale on every `ScanPowerControl`. Selected initial values may be
zero. Its derivative is `dP/dx = reference_w`, including at zero, so an inactive
heater can receive a nonzero activation gradient. Log coordinates never reach
zero; squared coordinates lose that activation gradient at zero.

Use the linear option through `ConstrainedScanStudy`. Each selected power adds
one automatic `-x <= 0` inequality to the native SQP, before all declared power,
geometry, energy, point-state and sampled-peak rows. These automatic rows count
against `max_kkt_dimension`. Their physical watts are reported separately in
`ScanConstraintEvaluation::nonnegative_power_values_w`. No additional user
constraint is needed just to admit nonnegative heating. A selected linear power
can also have an explicit upper bound of zero; the logarithmic map cannot.

Objective, state, peak and integrated-energy gradients all use the same selected
coordinate map. Geometry remains in the original signed-center/log-width
coordinates. Unselected power values, including final shutdowns, stay exactly
as authored. Full-storage and checkpointed adjoints retain their existing
physical solvers, event sides and replay checks.

A finite negative selected linear-power trial is unavailable to SQP before any
source rebuild or physics call. The native SQP can backtrack that known domain
failure. Other source, phase, material, numerical and cancellation errors still
propagate. Public evaluations of negative power return errors; sources are
never clipped to zero, and a positive value that underflows is refused. Native
SQP callback counts include unavailable probes; physical evaluation counts do
not include probes rejected before physics. All attempted physical evaluations
remain charged on failure. The unconstrained L-BFGS entry points reject a study
with selected nonnegative-linear powers rather than silently dropping its bounds.

A zero-power source still needs a mesh-resolved footprint for its power
sensitivity. This changes neither source support admission nor caloric smoothness
rules. A selected off state is exactly representable, but finite-tolerance SQP may
return a small positive power instead; reports and exports do not snap it to zero.

## Executable example

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --heater-off
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --heater-off 152 157 0 0
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --heater-off 149 148 1 1
cargo test --manifest-path examples/scan_power_fit/Cargo.toml -- zero_power
```

The synthetic unit-volume, unit-reference-density body starts at 150 J/kg on a
310 K latent plateau. Two powers act over [0,1) and [1,2), followed by a fixed
zero shutdown. The default targets are 154 and 152 J/kg, starting with both
powers off. With no cooling mechanism those targets cannot both be reached:
the analytical least-squares optimum is 3 W followed by 0 W, producing 153 J/kg
at both observations and objective 0.01. Targets 152/157 instead require 2/5 W;
targets 149/148 have the best admissible solution 0/0 W and objective 0.025.

The example solves the native enthalpy objective and SQP, exercises budget
continuation, uses eight endpoints with two-step replay sections, and replays
exported sources. It prints actual watts (without rounding them to off), original
predictions, objective, domain multipliers and KKT diagnostics. An unconverged
or infeasible candidate exits unsuccessfully. A converged constrained best fit
can retain nonzero target mismatch, which is reported rather than called recovery.

This is fixed-grid local inverse design, not binary heater scheduling, physical
switching-time optimization, experimental calibration or a global-optimum proof.
Native Rust execution must be distinguished from independent analytical checks.
