//! Shared constitutive and executor owners for the reduced deposition network.
//! The spatial P1 enthalpy solver remains stationary; this adapter uses its
//! existing material-chart owner, not a second phase or heat-capacity model.

/// Implicit, radiation-capable parcel transport using the existing nonlinear
/// enthalpy and conservative finite-contact solvers. Its aggregate Biot gate
/// refuses outside the isothermal regime; it is not a growing spatial P1 solve.
pub mod implicit;

use fs_exec::Cx;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;

use crate::transient::deposition::{
    DepositionError, DepositionSimulation, EnergyAudit, StepConfig, ThermalLaw,
};

impl ThermalLaw for EquilibriumEnthalpyPhaseCurve {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.state_at_specific_enthalpy(specific_h_j_kg)
            .map(|state| state.temperature_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        // The chart is piecewise linear in h; plateaus contribute zero. Unlike
        // a local derivative query this remains valid across all chart knots.
        self.knots()
            .windows(2)
            .map(|pair| {
                (pair[1].temperature_k - pair[0].temperature_k)
                    / (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
            })
            .fold(0.0, f64::max)
    }
}

/// Advance within an existing executor scope. Cancellation discards the whole
/// candidate, including newly activated material and event cursor updates.
pub fn advance_to<'a>(
    cx: &Cx<'_>,
    simulation: &DepositionSimulation<'a>,
    target_s: f64,
    config: StepConfig,
) -> Result<(DepositionSimulation<'a>, EnergyAudit), DepositionError> {
    simulation.advance_to(target_s, config, || {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)
    })
}
