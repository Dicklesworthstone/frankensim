# Event-limited adaptive thermal histories

`ThermalHistoryRun`, `EnthalpyHistoryRun`, `RadiativeThermalHistoryRun` and
`RadiativeEnthalpyHistoryRun` expose `advance_adaptive_step`. Each compares one
original backward-Euler step with two half steps from the same accepted state.
Only the actual two-half-step trajectory is published, never an extrapolated
field. Both original half-step reports are returned.

## Temperature and phase controls

After constructing a run and its original solver configuration, select explicit
local error and work limits:

```rust
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
let adaptive = AdaptiveHistoryConfig {
    absolute_temperature_k: 0.002,
    relative_temperature: 0.0,
    min_step_s: 1e-7,
    max_step_s: 0.25,
    max_attempts: 24,
};
```

Call the temperature or radiative temperature run with
`advance_adaptive_step(cx, target_s, proposed_step_s, step_config, adaptive)`.
After success, use `suggested_step_s` for the next proposal. Repeat until the
accepted time reaches the target. For radiation wrappers inspect the clock
through `conduction().time_s()`; other runs expose `time_s()` directly.

Enthalpy runs also take a final `PhaseErrorTolerance` argument from
`bc::history::march::enthalpy::adaptive`. Its fields are
`absolute_enthalpy_j_kg`, `relative_enthalpy` and `absolute_liquid_fraction`.
Temperature, specific enthalpy and liquid fraction must all pass. A latent
plateau cannot hide time error just because its temperature remains unchanged.
Relative allowances scale changes from accepted history, not an arbitrary
enthalpy origin or the absolute Kelvin offset.

The original linear/nonlinear method, reference masses or per-element capacity,
conductivity assignments, phase charts and contact interfaces remain attached.
Radiation wrappers retain every sourced emissivity and reservoir history and
return separate convection and nonlinear radiation reports. All source,
boundary, reference-temperature and attached reservoir knots cap the proposed
interval. Held data use exact interval-end left limits, preserving short pulses.

## Work and publication

Each successful call consumes two lifetime physical-step slots. A maximum of
`3 * max_attempts` endpoint solves is allowed per call, each with its original
solver work budget. Coarse and rejected trials consume work but never advance
physical state. Cancellation, invalid histories, failed physical solves or
unmet temporal accuracy publish no partial pair. Budget extension retains the
accepted fields and time. No hidden controller state is required for resume.

Sum both half-step storage changes and time-integrated supplies; do not sum
unweighted heat rates. Radiative physical balances must be read from the
physical radiation reports, not relabeled frozen inner Robin reports.

The controller halves rejected intervals and suggests growth only when the
normalized discrepancy is at most 1/8. An event or final time may shorten an
interval below the configured floor; it is still tested and is not accepted
merely because it is small. Intervals without a representable interior midpoint
refuse; direct event stepping remains available for adjacent-float intervals.

This is an estimated local time-error test, not a global error bound, certified
phase-transition time or continuum accuracy claim. Material corners and sharp
transients can violate smooth asymptotic assumptions. Physical solver and
energy tolerances should be tighter than temporal allowances; those gates are
never widened to make temporal adaptation succeed. Native Rust execution must
be reported separately from independent equation checks.
