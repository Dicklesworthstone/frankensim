//! Exact fixed-density Cv(T) storage and k(T) on the native thermal solver.
//!
//! The original caloric law stores internal energy u, with du/dT=Cv. This
//! adapter supplies dT/du=1/Cv and dk/du=(dk/dT)/Cv to the existing implicit
//! chain, graph, SDIRK2 and adaptive owners. Their historical h-labelled
//! coordinates retain the original InternalEnergy storage basis throughout.
//!
//! Conductivity uses the same material card, fixed density and full source
//! context as Cv. Every k breakpoint is retained independently of the Cv grid.
//! The original law still validates cell mass/volume and owns entropy; no
//! Cp substitution, pressure path, flow work or moving geometry is inferred.

use super::super::{ThermalTransportError, conductivity::PreparedScalarConductivity};
use super::IsochoricHeatCapacityLaw;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// One source-exact Cv owner and its supported scalar conductivity path.
/// Queries allocate nothing and use the shared prepared conductivity evaluator.
#[derive(Clone, Debug)]
pub struct IsochoricConductionLaw<'a> {
    caloric: &'a IsochoricHeatCapacityLaw,
    conductivity: PreparedScalarConductivity,
    identity: ContentHash,
}

impl<'a> IsochoricConductionLaw<'a> {
    /// Select positive scalar static conductivity along the original isochore.
    /// Its source must explicitly support the caloric density at every T.
    /// Fixed-pressure-only conductivity and other unsupported contexts refuse.
    pub fn try_new(
        card: &MaterialCard,
        caloric: &'a IsochoricHeatCapacityLaw,
        conductivity: &PropertyKey,
        selection: ClaimSelection,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        if card.content_hash() != caloric.material_card_identity() {
            return Err(invalid(
                "conduction must retain the original isochoric material card",
            ));
        }
        let [lower, upper] = caloric.query_endpoints();
        let conductivity = PreparedScalarConductivity::try_new(
            card,
            conductivity,
            caloric.temperature_axis(),
            lower,
            upper,
            selection,
        )?;
        let claim = card
            .claims()
            .claim(conductivity.source_line().knots[0].1.receipt.selected)
            .ok_or(invalid("selected conductivity source is missing"))?;
        if !claim.validity.bounds().contains_key(caloric.density_axis()) {
            return Err(invalid(
                "conductivity source must support the original isochoric density",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.isochoric-conduction.v1");
        hash.update(caloric.identity().as_bytes());
        conductivity.append_identity(&mut hash);
        Ok(Self {
            caloric,
            conductivity,
            identity: hash.finalize(),
        })
    }

    /// The same Cv law registered on the thermal timeline, including its
    /// original internal-energy storage basis and cell density admission.
    #[must_use]
    pub const fn caloric(&self) -> &'a IsochoricHeatCapacityLaw {
        self.caloric
    }

    /// Original conductivity source receipts and uncertainty, without resampling.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        self.conductivity.source_line()
    }

    pub(in crate::thermal::transport) fn integrated_conductivity_curve(
        &self,
    ) -> Result<crate::thermal::conduction::ConductivityCurve, DepositionError> {
        self.conductivity.integrated_curve()
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Original scalar conductivity at an already resolved caloric temperature.
    ///
    /// The SAME prepared source path and temperature support are retained.
    /// Relative-energy consumers can reuse their original inversion without a
    /// second absolute-energy round trip or a second caloric source evaluation.
    pub fn conductivity_at_temperature_w_m_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        self.conductivity
            .at_temperature(temperature_k)
            .map(|(k, _)| k)
    }

    /// Fused source-exact state and derivatives with respect to specific u.
    /// At k corners the shared evaluator selects the right derivative, except
    /// at the upper endpoint where it uses the left derivative.
    pub fn state_at_specific_internal_energy(
        &self,
        u: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        let temperature = self.caloric.temperature_k(u).map_err(|_| {
            ThermalTransportError::Invalid("internal energy outside original isochoric law")
        })?;
        let cv = self.caloric.specific_heat_cv_j_kg_k(temperature)?;
        self.conductivity.state(temperature, cv)
    }
}

impl ConductivityLaw for IsochoricConductionLaw<'_> {
    fn conductivity_w_m_k(&self, u: f64) -> Result<f64, DepositionError> {
        let temperature = self.caloric.temperature_k(u)?;
        self.conductivity
            .at_temperature(temperature)
            .map(|(conductivity, _)| conductivity)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

impl ImplicitConductionLaw for IsochoricConductionLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.specific_internal_energy_span_j_kg()
    }
    fn state_at_specific_enthalpy(
        &self,
        u: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        self.state_at_specific_internal_energy(u)
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn time_observables(&self, u: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.caloric.temperature_k(u)?,
            phase_mass_fraction: None,
        })
    }
}

#[cfg(test)]
#[path = "conduction/tests.rs"]
mod tests;
