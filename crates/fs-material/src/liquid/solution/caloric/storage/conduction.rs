//! Equilibrium binary-liquid calorics on the existing implicit heat solver.
//!
//! Storage retains the original pure-component Cp sources and regular-solution
//! equilibrium. Its local inverse tangent includes the changing equilibrium
//! phase compositions, without resampling enthalpy or fitting an apparent Cp.
//! At a miscibility corner the caloric owner explicitly supplies a homogeneous
//! one-sided Newton tangent. It supplies that same approximate tangent inside
//! an unresolved numerical binodal bracket; neither case claims a unique Cp.
//! The original enthalpy residual still decides acceptance.
//!
//! Conductivity is an explicitly declared scalar bulk curve for this fixed
//! composition, pressure and reference geometry. It is NOT inferred from the
//! individual phase conductivities: a phase fraction alone does not determine
//! mixture morphology or interphase thermal resistance. This stationary,
//! co-temperature thermal-cell closure does not provide a two-phase fluid,
//! diffusion, interphase heat transfer or a phase-fraction error observable.

use super::BinaryLiquidThermalStorage;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Original equilibrium storage plus caller-declared bulk heat conduction.
/// The borrowed storage must be the same law registered on the thermal owner.
#[derive(Clone, Debug)]
pub struct BinaryLiquidConductionLaw<'a> {
    caloric: &'a BinaryLiquidThermalStorage,
    conductivity: ConductivityCurve,
}

impl<'a> BinaryLiquidConductionLaw<'a> {
    /// Bind a bulk conductivity chart covering the complete prepared thermal
    /// domain. The caller is responsible for its physical validity for this
    /// mixture and geometry; this constructor does not infer source evidence.
    pub fn try_new(
        caloric: &'a BinaryLiquidThermalStorage,
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

    /// Original, unmodified storage object used by the thermal timeline.
    pub const fn caloric(&self) -> &'a BinaryLiquidThermalStorage {
        self.caloric
    }

    /// Original declared bulk curve, including every conductivity breakpoint.
    pub const fn conductivity(&self) -> &ConductivityCurve {
        &self.conductivity
    }

    /// One equilibrium inversion supplies temperature and its phase-aware
    /// tangent. Conductivity and dk/dh are evaluated at that same temperature.
    /// No additional equilibrium solve or numerical differentiation is needed
    /// to obtain the local heat capacity from the accepted caloric state.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let state = self
            .caloric
            .caloric_state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let capacity = state
            .isobaric_heat_capacity_response()
            .map_err(|_| DepositionError::MaterialDomain)?
            .heat_capacity_j_per_k();
        let temperature_k = state.temperature().value();
        let temperature_slope_kg_k_j = self.caloric.mass_kg / capacity;
        let (conductivity_w_m_k, dk_dt) = self
            .conductivity
            .conductivity_and_temperature_derivative(temperature_k)?;
        let conductivity_slope_kg_w_m_k_j = dk_dt * temperature_slope_kg_k_j;
        if !temperature_slope_kg_k_j.is_finite()
            || temperature_slope_kg_k_j <= 0.0
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

impl ConductivityLaw for BinaryLiquidConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.conductivity
            .conductivity_w_m_k(self.caloric.temperature_k(h)?)
    }
}

impl ImplicitConductionLaw for BinaryLiquidConductionLaw<'_> {
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
            phase_mass_fraction: None,
        })
    }
}
