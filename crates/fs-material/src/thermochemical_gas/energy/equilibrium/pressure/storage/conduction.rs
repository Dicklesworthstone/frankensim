//! Reacting isobaric calorics on the original implicit heat operator.
//!
//! The caller declares bulk k(T) for these fixed element ratios and pressure.
//! Frozen transport values alone do not supply the reacting composition
//! derivative dk/dT. This adapter therefore uses one explicit original curve,
//! retaining all its knots, and differentiates through the reacting Cp.
//! The same borrowed caloric law must be registered on the thermal timeline.

use super::GasEquilibriumThermalStorage;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// One original reacting storage law and a declared bulk conductivity chart.
pub struct GasEquilibriumConductionLaw<'storage, 'source> {
    caloric: &'storage GasEquilibriumThermalStorage<'source>,
    conductivity: ConductivityCurve,
}
impl<'storage, 'source> GasEquilibriumConductionLaw<'storage, 'source> {
    /// Require the declared curve to cover the complete prepared caloric span.
    /// Applicability to the actual pressure and reacting mixture is explicit
    /// caller data; this constructor does not invent material evidence.
    pub fn try_new(
        caloric: &'storage GasEquilibriumThermalStorage<'source>,
        conductivity: ConductivityCurve,
    ) -> Result<Self, DepositionError> {
        let [lower, upper] = caloric.temperature_domain();
        let [k_lower, k_upper] = conductivity.temperature_range_k();
        if k_lower > lower.value() || k_upper < upper.value() {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Self { caloric, conductivity })
    }

    pub const fn caloric(&self) -> &'storage GasEquilibriumThermalStorage<'source> {
        self.caloric
    }
    pub const fn conductivity(&self) -> &ConductivityCurve { &self.conductivity }

    /// One original inverse supplies T and dT/dh; k and dk/dh use that state.
    pub fn state_at_specific_enthalpy(
        &self, h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let state = self.caloric.caloric_state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let temperature_k = state.temperature.value();
        let temperature_slope_kg_k_j =
            self.caloric.feed_mass_kg() / state.response.heat_capacity_p_j_k;
        let (conductivity_w_m_k, dk_dt) =
            self.conductivity.conductivity_and_temperature_derivative(temperature_k)?;
        let conductivity_slope_kg_w_m_k_j = dk_dt * temperature_slope_kg_k_j;
        if !temperature_slope_kg_k_j.is_finite() || temperature_slope_kg_k_j <= 0.0
            || !conductivity_slope_kg_w_m_k_j.is_finite()
            || (dk_dt != 0.0 && conductivity_slope_kg_w_m_k_j == 0.0) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitConductionState {
            temperature_k, temperature_slope_kg_k_j,
            conductivity_w_m_k, conductivity_slope_kg_w_m_k_j,
        })
    }
}
impl ConductivityLaw for GasEquilibriumConductionLaw<'_, '_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.conductivity.conductivity_w_m_k(self.caloric.temperature_k(h)?)
    }
}
impl ImplicitConductionLaw for GasEquilibriumConductionLaw<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw { self.caloric }
    fn conductivity_law(&self) -> &dyn ConductivityLaw { self }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] { self.caloric.specific_enthalpy_domain_j_kg() }
    fn state_at_specific_enthalpy(
        &self, h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Self::state_at_specific_enthalpy(self, h)
    }

    /// d(1/rho)/dh = (alpha/rho) dT/dh at the same fixed pressure and
    /// conserved inventory. Reaction-induced molecular-count changes are
    /// already included in alpha; no frozen ideal-gas 1/T shortcut is used.
    fn specific_volume_derivative_at_specific_enthalpy(
        &self, h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let state = self.caloric.caloric_state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let derivative = (state.response.thermal_expansion_per_k / state.response.density_kg_m3)
            * (self.caloric.feed_mass_kg() / state.response.heat_capacity_p_j_k);
        if !derivative.is_finite()
            || (state.response.thermal_expansion_per_k != 0.0 && derivative == 0.0) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(Some(derivative))
    }

    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.caloric.temperature_k(h)?,
            phase_mass_fraction: None,
        })
    }
}
