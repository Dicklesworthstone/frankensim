//! Rigid reacting storage on the original implicit thermal operator.
//!
//! Bulk k(T) is explicitly authored for the fixed conserved density and
//! elemental ratios. Its pressure/composition dependence along that path is
//! part of the supplied curve, not inferred from frozen transport snapshots.

use super::GasIsochoricEquilibriumThermalStorage;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{ImplicitConductionLaw, ImplicitConductionState};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// One original reacting internal-energy law and a declared bulk conductivity.
/// Register this same borrowed caloric law on the original thermal timeline.
pub struct GasIsochoricEquilibriumConductionLaw<'storage, 'source> {
    caloric: &'storage GasIsochoricEquilibriumThermalStorage<'source>,
    conductivity: ConductivityCurve,
}
impl<'storage, 'source> GasIsochoricEquilibriumConductionLaw<'storage, 'source> {
    /// Require the complete caloric temperature span to be supported by k(T).
    /// This admits numerical support, not empirical applicability of the curve.
    pub fn try_new(
        caloric: &'storage GasIsochoricEquilibriumThermalStorage<'source>,
        conductivity: ConductivityCurve,
    ) -> Result<Self, DepositionError> {
        let [lower, upper] = caloric.temperature_domain();
        let [k_lower, k_upper] = conductivity.temperature_range_k();
        if k_lower > lower.value() || k_upper < upper.value() {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Self {
            caloric,
            conductivity,
        })
    }

    /// The identical original reacting storage law to register on the timeline.
    pub const fn caloric(&self) -> &'storage GasIsochoricEquilibriumThermalStorage<'source> {
        self.caloric
    }
    /// Authored bulk conductivity along this fixed-density reacting branch.
    pub const fn conductivity(&self) -> &ConductivityCurve {
        &self.conductivity
    }

    /// One original UV query supplies T and dT/du; k and dk/du use that state.
    pub fn state_at_specific_internal_energy(
        &self,
        u: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let state = self
            .caloric
            .caloric_state_at_specific_internal_energy(u, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let temperature_k = state.temperature.value();
        let temperature_slope_kg_k_j =
            self.caloric.feed_mass_kg() / state.response.heat_capacity_v_j_k;
        let (conductivity_w_m_k, dk_dt) = self
            .conductivity
            .conductivity_and_temperature_derivative(temperature_k)?;
        let conductivity_slope_kg_w_m_k_j = dk_dt * temperature_slope_kg_k_j;
        if !temperature_slope_kg_k_j.is_finite()
            || temperature_slope_kg_k_j <= 0.0
            || !conductivity_slope_kg_w_m_k_j.is_finite()
            || (dk_dt != 0.0 && conductivity_slope_kg_w_m_k_j == 0.0)
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
impl ConductivityLaw for GasIsochoricEquilibriumConductionLaw<'_, '_> {
    fn conductivity_w_m_k(&self, u: f64) -> Result<f64, DepositionError> {
        self.conductivity
            .conductivity_w_m_k(self.caloric.temperature_k(u)?)
    }
}
impl ImplicitConductionLaw for GasIsochoricEquilibriumConductionLaw<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.specific_internal_energy_domain_j_kg()
    }
    fn state_at_specific_enthalpy(
        &self,
        u: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        self.state_at_specific_internal_energy(u)
    }
    fn specific_volume_derivative_at_specific_enthalpy(
        &self,
        u: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.caloric.bulk_density_kg_m3(u)?;
        // Sealed rigid density is constant even as molecular count and pressure
        // respond. The generic operator's historical h label carries u here.
        Ok(Some(0.0))
    }
}
