//! Native implicit conduction from the original volumetric-capacity source.
//!
//! The fixed-reference-volume law stores e = integral C_vol(T) dT / rho0.
//! Its exact local tangent is dT/de = rho0/C_vol(T), so the same source k(T)
//! gives dk/de = (dk/dT)*rho0/C_vol(T). The existing chain/graph solvers own
//! residuals, energy increments, time integration, sources and acceptance.
//!
//! This extends the original law rather than introducing a second caloric
//! object. Its mass/volume admission, energy datum, source receipts, optional
//! sensible entropy and DeclaredSensibleStorage basis remain unchanged.
//! Scalar conductivity retains the caller's direction and fixed-geometry
//! interpretation. No density measurement, Cp/Cv conversion, phase fraction,
//! stream enthalpy or moving-volume work is inferred from volumetric capacity.

use super::{ThermalTransportError, timeline::FixedVolumeThermalLaw};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

impl FixedVolumeThermalLaw<'_> {
    /// Evaluate original storage, conductivity and exact caloric tangents once.
    ///
    /// Original affine C_vol is integrated and inverted analytically; it is
    /// never frozen at the old temperature or replaced by enthalpy chords.
    /// Interior conductivity knots use their right derivative and the upper
    /// endpoint its left derivative. Successful queries allocate nothing and
    /// perform no source selection or hashing. All original support and
    /// representability refusals remain active.
    pub fn conduction_state_at_specific_storage(
        &self,
        storage_j_kg: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let temperature_k = self
            .temperature_k(storage_j_kg)
            .map_err(|_| invalid("storage outside original fixed-volume law"))?;
        let ([capacity, conductivity], dk_dt) = self.path().conduction_values(temperature_k)?;
        let temperature_slope = self.reference_density_kg_m3() / capacity;
        let conductivity_slope = dk_dt * temperature_slope;
        if !temperature_slope.is_finite()
            || temperature_slope <= 0.0
            || !conductivity_slope.is_finite()
            || (dk_dt != 0.0 && conductivity_slope == 0.0)
        {
            return Err(invalid("volumetric storage tangents are not representable"));
        }
        Ok(ImplicitConductionState {
            temperature_k,
            temperature_slope_kg_k_j: temperature_slope,
            conductivity_w_m_k: conductivity,
            conductivity_slope_kg_w_m_k_j: conductivity_slope,
        })
    }
}

impl ConductivityLaw for FixedVolumeThermalLaw<'_> {
    fn conductivity_w_m_k(&self, storage_j_kg: f64) -> Result<f64, DepositionError> {
        let temperature = self.temperature_k(storage_j_kg)?;
        self.path()
            .at_temperature(temperature)
            .map(|state| state.conductivity_w_m_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

impl ImplicitConductionLaw for FixedVolumeThermalLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.specific_storage_span_j_kg()
    }

    fn state_at_specific_enthalpy(
        &self,
        storage_j_kg: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        self.conduction_state_at_specific_storage(storage_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn time_observables(
        &self,
        storage_j_kg: f64,
    ) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.temperature_k(storage_j_kg)?,
            phase_mass_fraction: None,
        })
    }
}
