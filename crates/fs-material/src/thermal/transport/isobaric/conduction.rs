//! Source-exact isobaric storage and conduction on the existing thermal owner.
//!
//! The borrowed caloric law integrates its original affine Cp(T) exactly. This
//! adapter compiles conductivity on the same material card, temperature span,
//! pressure and fixed composition. It supplies dT/dh=1/Cp and dk/dh=(dk/dT)/Cp
//! to the existing implicit chain/graph and their adaptive time integrators.
//! No resampled enthalpy chart, effective capacity or second solver is created.
//!
//! Conductivity is a positive scalar for the caller's declared heat-flow
//! direction. Fixed reference geometry, single phase and the original isobaric
//! interpretation still apply. Spatial resolution and temporal accuracy are
//! properties of the thermal solve, not of this constitutive adapter.

use super::super::{ThermalTransportError, conductivity::PreparedScalarConductivity};
use super::IsobaricHeatCapacityLaw;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

#[cfg(test)]
use super::super::CONDUCTIVITY;
#[cfg(test)]
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

/// One original caloric owner and its continuously supported conductivity path.
/// Construction resolves source evidence once; runtime queries allocate nothing.
#[derive(Clone, Debug)]
pub struct IsobaricConductionLaw<'a> {
    caloric: &'a IsobaricHeatCapacityLaw,
    conductivity: PreparedScalarConductivity,
    identity: ContentHash,
}

impl<'a> IsobaricConductionLaw<'a> {
    /// Select conductivity from the exact original card and fixed source context.
    /// Every conductivity breakpoint is retained independently of the Cp knots.
    /// Scalar, explicitly directional static conductivity is admitted; tensor,
    /// dynamic, pressure-unbounded and out-of-support sources refuse.
    pub fn try_new(
        card: &MaterialCard,
        caloric: &'a IsobaricHeatCapacityLaw,
        conductivity: &PropertyKey,
        selection: ClaimSelection,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        if card.content_hash() != caloric.material_card_identity() {
            return Err(invalid(
                "conduction must retain the original isobaric material card",
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
        if !claim
            .validity
            .bounds()
            .contains_key(caloric.pressure_axis())
        {
            return Err(invalid(
                "conductivity source must support the original isobaric pressure",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.isobaric-conduction.v1");
        hash.update(caloric.identity().as_bytes());
        conductivity.append_identity(&mut hash);
        Ok(Self {
            caloric,
            conductivity,
            identity: hash.finalize(),
        })
    }

    /// The exact caloric law registered on the thermal timeline.
    #[must_use]
    pub const fn caloric(&self) -> &'a IsobaricHeatCapacityLaw {
        self.caloric
    }

    /// Original material-source receipts and uncertainty, without resampling.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        self.conductivity.source_line()
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// One original inverse and analytic derivatives on the same temperature.
    /// Interior conductivity knots select their right derivative; the upper
    /// endpoint uses the left derivative. Cp remains continuous at its knots.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        let temperature = self.caloric.temperature_k(h).map_err(|_| {
            ThermalTransportError::Invalid("enthalpy outside original isobaric law")
        })?;
        let cp = self.caloric.specific_heat_cp_j_kg_k(temperature)?;
        self.conductivity.state(temperature, cp)
    }
}

impl ConductivityLaw for IsobaricConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = self.caloric.temperature_k(h)?;
        self.conductivity
            .at_temperature(temperature)
            .map(|(conductivity, _)| conductivity)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

impl ImplicitConductionLaw for IsobaricConductionLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.specific_enthalpy_span_j_kg()
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Self::state_at_specific_enthalpy(self, h).map_err(|_| DepositionError::MaterialDomain)
    }

    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.caloric.temperature_k(h)?,
            phase_mass_fraction: None,
        })
    }
}

#[cfg(test)]
mod tests;
