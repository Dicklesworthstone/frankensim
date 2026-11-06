//! Simultaneous wall/film/gas steps on the original bounded window clock.

use super::*;
use crate::finite_gas_thermal::wall_film::implicit::{WallFilmImplicitControls, WorkMeter};

/// Uniform-source window controls and numerical allowances shared by EVERY
/// implicit attempt, including rejected phase-admission proposals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFilmImplicitWindowControls {
    pub window: WallFilmWindowSolve,
    /// Newton iteration and residual evaluation limits apply to the complete
    /// requested window. The other fields select the original numerical policy.
    pub solve: WallFilmImplicitControls,
    /// Shared gas inversion allowance. Each invocation also retains the smaller
    /// original `environment.solve.phase.gas.gas.maximum_iterations` bound.
    pub maximum_gas_iterations: usize,
}
impl WallFilmImplicitWindowControls {
    pub fn new(window: WallFilmWindowSolve) -> Self {
        Self {
            window,
            solve: WallFilmImplicitControls::default(),
            maximum_gas_iterations: 100_000,
        }
    }
}

/// Complete physical accounts and all attempted numerical work. No rejected
/// proposal contributes physical heat, transferred mass, events or elapsed time.
#[derive(Clone, Debug, PartialEq)]
pub struct WallFilmImplicitWindowReport {
    pub physical: WallFilmWindowReport,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    /// Known original gas iterations from successful finalizers, including
    /// proposals subsequently rejected by the returned-wall phase check.
    pub gas_iterations: usize,
    /// Gas work charged to the window. A failed finalizer has no partial receipt,
    /// so its reserved bounded loop remains charged; successful calls refund
    /// unused iterations. This can exceed known `gas_iterations` after retries.
    pub charged_gas_iterations: usize,
    /// Accepted intervals' absolute native phase plus gas-reconstruction mass
    /// defects. Checked alongside actual final volatile, total and carrier mass.
    pub sum_absolute_mass_residual_kg: f64,
}

impl WallFilmThermalEnvironment<'_> {
    /// Continue the simultaneous wall, zero-capable film and finite gas endpoint
    /// over a complete physical window with the original source allocation.
    ///
    /// Body and gas input joules are full-window totals, uniform in time. Only
    /// explicit original phase `StepSizeLimit` refusals subdivide an interval;
    /// source/numerical errors and exhausted shared work abort atomically even
    /// after earlier internal intervals succeeded. Accepted steps preserve both
    /// source-owned liquid calorics and native latent storage in the wall.
    ///
    /// This is an admission-driven first-order continuation. Refine the physical
    /// step for time accuracy; its frozen-endpoint dry-out times are estimates.
    pub fn advance_implicit_window(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: WallFilmImplicitWindowControls,
    ) -> Result<UniformEnthalpyStep<WallFilmImplicitWindowReport>, WallFilmWindowError> {
        poll(cx)?;
        let meter = WorkMeter::new(controls.solve, controls.maximum_gas_iterations)?;
        let (candidate, sum_absolute_mass_residual_kg) =
            self.advance_window_with(cx, input, controls.window, |environment, part| {
                environment
                    .advance_implicit_metered(cx, part, controls.solve, &meter)
                    .map(WindowStep::implicit)
                    .map_err(WallFilmWindowError::from)
            })?;
        Ok(UniformEnthalpyStep {
            state: candidate.state,
            external_heat_j: candidate.external_heat_j,
            energy_residual_tolerance_j: candidate.energy_residual_tolerance_j,
            report: WallFilmImplicitWindowReport {
                physical: candidate.report,
                newton_iterations: meter.newton_iterations(),
                residual_evaluations: meter.residual_evaluations(),
                gas_iterations: meter.gas_iterations(),
                charged_gas_iterations: meter.charged_gas_iterations(),
                sum_absolute_mass_residual_kg,
            },
        })
    }
}
