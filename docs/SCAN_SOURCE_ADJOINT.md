# Moving-source controls through implicit thermal trajectories

`MovingHeatSource::linearize_at` and `linearize_before` bind the original
normalized P1 heat source at an explicit time. `MovingHeatLinearization::apply`
returns the source-density Jacobian action; `pullback` consumes one objective
seed per source-density node and returns center, half-width and absorbed-watt
derivatives. Both differentiate the normalization. Translation and width
changes redistribute the declared absorbed watts; they do not change their
integral, even on a clipped footprint. A P1 nodal heat-load seed is a different
coordinate system and must not be supplied as a density seed.

`ScanningEnthalpyRun::advance_linearized_to` runs the original nonlinear primal
once and attaches the existing revalidated physical enthalpy tangent. With
radiation enabled, the tangent includes the full mean-temperature-dependent
radiation feedback, not a frozen Robin coefficient. The response retains the
original source observations, energy checks, phase fields and radiation report.
`advance_linearized_event_step` additionally respects schedule and selected
motion limits, using the identical pre-switch loads and derivative side.

```rust
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{
    ScanSensitivityBudget, ScanSensitivityError, ScanStepGradient,
};
use fs_conduction::transient::enthalpy::EnthalpyStepConfig;
use fs_conduction::LinearConfig;
use fs_exec::Cx;

pub fn final_temperature_gradient(
    cx: &Cx<'_>,
    run: &mut ScanningEnthalpyRun<'_, '_, '_>,
    endpoint_s: f64,
    step: EnthalpyStepConfig,
    budget: ScanSensitivityBudget,
    temperature_seed: &[f64],
    transpose: LinearConfig,
) -> Result<ScanStepGradient, ScanSensitivityError> {
    let endpoint = run.advance_linearized_to(cx, endpoint_s, step, budget)?;
    endpoint.pullback_temperature(cx, temperature_seed, transpose)
}
```

For an objective in nodal J/kg coordinates, call `pullback` instead. One checked
transpose solve supplies every source's seven instantaneous control derivatives.
The returned `storage.transport.previous_specific_enthalpy` is the seed for the
previous fixed-grid endpoint. Retain a bounded collection of bindings, visit it
backward, and pass this seed to the preceding binding's `pullback`. Sum the
per-step derivatives for a shared constant power or width. For ramp/table/path
parameters, additionally apply the chosen fixed-time interpolation chain rule;
that parameterization is not inferred by this API. The run does not build or
retain a growing tape automatically.

`ScanSensitivityBudget::max_source_vertices` bounds retained source count times
node count (up to three scalar arrays per binding, excluding allocator overhead).
`max_radiation_feedback_entries` is passed unchanged to the existing radiative
adjoint owner. Original spatial, nonlinear, energy, source-support and transpose
budgets remain independent. Failed advance/binding leaves all accepted fields,
time and the step count unchanged. A subsequent invalid pullback seed does not
undo an already accepted primal step.

These are derivatives of the fixed discrete endpoint and its original phase
charts. Phase/material corners or validity endpoints refuse a classical
linearization. On a latent plateau, temperature derivatives vanish while
enthalpy derivatives still describe heat storage. An off source must have an
admitted unit-power footprint to differentiate power; negative power directions
are infeasible at zero. Source-kernel support entry/exit is C1, not necessarily C2.

There is no derivative of event times, adaptive acceptance, motion-limit
selection, support-count admission, mesh geometry, optical absorption, phase-chart
parameters or material selection. Radiation feedback is included, but reservoir
and emissivity parameter gradients are not emitted by this source-control API.
The reported transpose residual is not a continuum or global time-error bound.
