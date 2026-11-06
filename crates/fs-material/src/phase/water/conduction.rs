//! Original IF97 liquid/water/steam calorics on the native implicit heat solver.
//!
//! The immutable pressure chart owns h, T, sensible Cp, latent inventory and
//! the continuous dh/T entropy primitive. A fused state query performs one
//! original enthalpy inversion and uses its analytic Cp; saturation retains
//! zero dT/dh, including the explicit plateau-side tangent at both endpoints.
//! No resampled enthalpy curve, finite apparent Cp or extra latent source is used.
//!
//! Conductivity is an explicitly declared scalar BULK k(T) for this fixed
//! pressure and reference geometry, including the coexistence interval. This
//! closure makes k independent of vapor inventory at Ts. The caller must justify
//! that effective law: it is not inferred from IF97, phase conductivities or a
//! phase fraction. Convection, boiling/condensation transfer coefficients,
//! evolving phase geometry and separate phase temperatures require other owners.
//! Geometry remains reference heat-exchange geometry, not a rigid sealed vessel.

use super::IsobaricWater;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Borrow the actual water caloric owner with a declared bulk conductivity.
/// The same borrowed chart must be registered on the thermal timeline.
#[derive(Clone, Debug)]
pub struct WaterConductionLaw<'a> {
    caloric: &'a IsobaricWater,
    conductivity: ConductivityCurve,
}

impl<'a> WaterConductionLaw<'a> {
    /// Bind k(T) over the complete original chart. This checks support, not
    /// measured validity or the physical suitability of the bulk closure.
    pub fn try_new(
        caloric: &'a IsobaricWater,
        conductivity: ConductivityCurve,
    ) -> Result<Self, DepositionError> {
        let [lower, upper] = conductivity.temperature_range_k();
        if lower > caloric.controls.minimum_temperature_k
            || upper < caloric.controls.maximum_temperature_k
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Self {
            caloric,
            conductivity,
        })
    }

    /// Original chart, with its pressure, source states, Cp bounds and entropy.
    pub const fn caloric(&self) -> &'a IsobaricWater {
        self.caloric
    }

    /// Unmodified caller-declared bulk curve, including every breakpoint.
    pub const fn conductivity(&self) -> &ConductivityCurve {
        &self.conductivity
    }

    /// One original inverse supplies T and dT/dh; the original k(T) segment
    /// supplies k and `(dk/dT)(dT/dh)` at that same state. Sensible Cp is never
    /// frozen at the initial temperature or replaced with an enthalpy chord.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let state = self
            .caloric
            .state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?
            .into_state();
        let temperature_k = state.temperature_k();
        let temperature_slope_kg_k_j = state
            .isobaric_heat_capacity_response()
            .temperature_enthalpy_derivative_kg_k_j();
        let (conductivity_w_m_k, dk_dt) = self
            .conductivity
            .conductivity_and_temperature_derivative(temperature_k)?;
        let conductivity_slope_kg_w_m_k_j = dk_dt * temperature_slope_kg_k_j;
        if !temperature_slope_kg_k_j.is_finite()
            || temperature_slope_kg_k_j < 0.0
            || !conductivity_slope_kg_w_m_k_j.is_finite()
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitConductionState {
            temperature_k,
            temperature_slope_kg_k_j,
            conductivity_w_m_k,
            conductivity_slope_kg_w_m_k_j,
        })
    }
}

impl ConductivityLaw for WaterConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.conductivity
            .conductivity_w_m_k(self.caloric.temperature_k(h)?)
    }
}

impl ImplicitConductionLaw for WaterConductionLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.enthalpy_bounds_j_kg()
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Self::state_at_specific_enthalpy(self, h)
    }

    /// Temperature and VAPOR mass fraction from one original caloric inverse.
    /// A latent step can change inventory while leaving temperature unchanged,
    /// so adaptive time control must not infer its phase error from T alone.
    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        let state = self
            .caloric
            .state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?
            .into_state();
        Ok(ImplicitConductionObservables {
            temperature_k: state.temperature_k(),
            phase_mass_fraction: Some(state.vapor_mass_fraction()),
        })
    }
}

#[cfg(test)]
mod tests;
