//! Constitutive adapters for the shared conservative thermal timeline.
//!
//! The material remains the only owner of the enthalpy/phase relation. This
//! implementation makes the same immutable chart available to electrical,
//! thermal and mechanical consumers without a dependency on a FEM solver.

use crate::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

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
