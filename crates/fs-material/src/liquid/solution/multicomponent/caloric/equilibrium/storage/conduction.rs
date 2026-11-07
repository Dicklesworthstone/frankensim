//! Equilibrium liquid redistribution on the native implicit heat operator.
//!
//! The original caloric branch owns temperature, phase composition, enthalpy
//! and the unregularized redistribution heat capacity. One inversion supplies
//! both temperature and its exact local inverse tangent. The physical energy
//! residual always uses the original branch, including its source Cp knots.
//!
//! Bulk conductivity is explicitly authored for this composition, pressure and
//! reference geometry. Phase conductivities alone do not determine a separated
//! mixture's morphology or interphase resistance. This adapter therefore never
//! averages them or treats the caloric phase weights as a conductivity law.
//! Cells share a local temperature and retain their total component inventory;
//! no phase-selective flow, interphase heat-transfer rate or momentum is added.

use super::MulticomponentLiquidEquilibriumThermalStorage;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Original fixed-feed equilibrium storage and a declared bulk conductivity.
/// The same borrowed storage must be registered on the thermal timeline.
#[derive(Clone, Debug)]
pub struct MulticomponentLiquidEquilibriumConductionLaw<'a> {
    caloric: &'a MulticomponentLiquidEquilibriumThermalStorage,
    conductivity: ConductivityCurve,
}

impl<'a> MulticomponentLiquidEquilibriumConductionLaw<'a> {
    /// Bind a bulk curve over the complete prepared caloric domain. Its
    /// applicability to this mixture and morphology remains the caller's
    /// physical input; this constructor does not infer material evidence.
    pub fn try_new(
        caloric: &'a MulticomponentLiquidEquilibriumThermalStorage,
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

    pub const fn caloric(&self) -> &'a MulticomponentLiquidEquilibriumThermalStorage {
        self.caloric
    }

    /// Original declared bulk conductivity, retaining all of its breakpoints.
    pub const fn conductivity(&self) -> &ConductivityCurve {
        &self.conductivity
    }

    /// One original equilibrium inverse supplies T and dT/dh = mass/Cp.
    /// Conductivity and dk/dh use that same state; no finite differencing or
    /// second phase solve is needed to obtain the implicit Jacobian entries.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let state = self
            .caloric
            .caloric_state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let temperature_k = state.temperature().value();
        let capacity = state.isobaric_heat_capacity_j_per_k();
        let temperature_slope_kg_k_j = self.caloric.mass_kg() / capacity;
        if temperature_slope_kg_k_j <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        conductivity_state(&self.conductivity, temperature_k, temperature_slope_kg_k_j)
    }
}

// Both original sensible branches and exact coexistence plateaus use this
// same chain rule. A plateau legitimately has dT/dh=0; a resolved branch's
// caller separately requires its strictly positive inverse heat capacity.
pub(super) fn conductivity_state(
    conductivity: &ConductivityCurve,
    temperature_k: f64,
    temperature_slope_kg_k_j: f64,
) -> Result<ImplicitConductionState, DepositionError> {
    let (conductivity_w_m_k, dk_dt) =
        conductivity.conductivity_and_temperature_derivative(temperature_k)?;
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

impl ConductivityLaw for MulticomponentLiquidEquilibriumConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.conductivity
            .conductivity_w_m_k(self.caloric.temperature_k(h)?)
    }
}

impl ImplicitConductionLaw for MulticomponentLiquidEquilibriumConductionLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.specific_enthalpy_domain_j_kg()
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Self::state_at_specific_enthalpy(self, h)
    }

    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.caloric.temperature_k(h)?,
            // Several liquid phase fractions cannot be represented by the
            // native scalar solid/liquid-fraction observable.
            phase_mass_fraction: None,
        })
    }
}
