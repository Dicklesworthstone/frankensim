//! Native implicit conduction across original sensible and latent liquid heat.
//!
//! The original joined calorics supply T and dT/dh in one inversion. The
//! genuine coexistence interval has exactly zero dT/dh and dk/dh. Its evolving
//! upper-ensemble mass fraction remains an explicit adaptive error observable,
//! even when temperature is constant. Both ensembles are liquid: this fraction
//! is neither a solid/liquid fraction nor vapor quality.
//!
//! Conductivity is a caller-authored bulk k(T) over the complete path. Its
//! constant value at coexistence is that explicit closure, not an average of
//! phase conductivities or a claimed morphology/interphase transfer model.

use super::super::conduction::conductivity_state;
use super::MulticomponentLiquidCoexistenceThermalStorage;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Original joined storage and an explicitly declared bulk conductivity curve.
#[derive(Clone, Debug)]
pub struct MulticomponentLiquidCoexistenceConductionLaw<'a> {
    caloric: &'a MulticomponentLiquidCoexistenceThermalStorage,
    conductivity: ConductivityCurve,
}

impl<'a> MulticomponentLiquidCoexistenceConductionLaw<'a> {
    /// The bulk curve must cover both full sensible branches and their common
    /// temperature. Physical applicability to the phase morphology, pressure
    /// and reference geometry remains an explicit input of its author.
    pub fn try_new(
        caloric: &'a MulticomponentLiquidCoexistenceThermalStorage,
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

    pub const fn caloric(&self) -> &'a MulticomponentLiquidCoexistenceThermalStorage {
        self.caloric
    }

    pub const fn conductivity(&self) -> &ConductivityCurve {
        &self.conductivity
    }

    /// One original heat-coordinate query supplies the branch tangent or
    /// exact latent zero. No apparent heat capacity or numerical differencing
    /// enters the original native heat residual/Jacobian.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let state = self
            .caloric
            .caloric_state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        conductivity_state(
            &self.conductivity,
            state.temperature().value(),
            state.temperature_slope_kg_k_j(),
        )
    }
}

impl ConductivityLaw for MulticomponentLiquidCoexistenceConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.conductivity
            .conductivity_w_m_k(self.caloric.temperature_k(h)?)
    }
}

impl ImplicitConductionLaw for MulticomponentLiquidCoexistenceConductionLaw<'_> {
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
        let state = self
            .caloric
            .caloric_state_at_specific_enthalpy(h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(ImplicitConductionObservables {
            temperature_k: state.temperature().value(),
            phase_mass_fraction: Some(state.upper_branch_mass_fraction()),
        })
    }
}
