//! Original mixture conductivity on homogeneous multicomponent enthalpy.
//!
//! The existing scalar source-line owner compiles the already selected liquid
//! conductivity at the actual fixed composition and pressure. All k(T) knots,
//! receipts and source support survive; pure-fluid averaging is not used.
//! One original caloric inverse supplies T and Cp for exact dT/dh and dk/dh.

use super::MulticomponentLiquidThermalStorage;
use crate::thermal::transport::{ThermalTransportError, conductivity::PreparedScalarConductivity};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyLineSlice};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Original homogeneous storage and original pinned mixture conductivity.
#[derive(Clone, Debug)]
pub struct MulticomponentLiquidConductionLaw<'a> {
    caloric: &'a MulticomponentLiquidThermalStorage,
    conductivity: PreparedScalarConductivity,
}

impl<'a> MulticomponentLiquidConductionLaw<'a> {
    /// Admit a continuous conductivity path on the original mixture card.
    /// The retained liquid conductivity role and its selected source determine
    /// the key; this cannot silently select a replacement conductivity claim.
    /// The source must explicitly support the original static pressure.
    pub fn try_new(
        card: &MaterialCard,
        caloric: &'a MulticomponentLiquidThermalStorage,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let reference = &caloric.owner.liquid.reference;
        if card.content_hash() != reference.material().card_identity() {
            return Err(invalid(
                "multicomponent conduction requires the original mixture card",
            ));
        }
        let role = &reference.requirements().properties()[2];
        let selected = reference
            .material()
            .property(role.name())
            .ok_or(invalid("original liquid conductivity role is missing"))?
            .answer()
            .receipt
            .selected;
        let claim = card
            .claims()
            .claim(selected)
            .ok_or(invalid("original liquid conductivity source is missing"))?;
        if !claim
            .validity
            .bounds()
            .contains_key(&caloric.owner.pressure_axis)
        {
            return Err(invalid(
                "liquid conductivity must support the original static pressure",
            ));
        }
        let conductivity = PreparedScalarConductivity::try_new(
            card,
            &claim.key,
            &caloric.owner.liquid.temperature_axis,
            &caloric.source_endpoints[0],
            &caloric.source_endpoints[1],
            ClaimSelection::Pinned(selected),
        )?;
        Ok(Self {
            caloric,
            conductivity,
        })
    }

    pub const fn caloric(&self) -> &'a MulticomponentLiquidThermalStorage {
        self.caloric
    }
    /// Original conductivity evidence, independently of pure caloric Cp knots.
    pub const fn source_line(&self) -> &PropertyLineSlice {
        self.conductivity.source_line()
    }
    /// Original inversion plus exact analytic local tangents. At k corners the
    /// derivative is right-sided, except at the upper endpoint (left-sided).
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        let state = self
            .caloric
            .caloric_state_at_specific_enthalpy(h, || false)
            .map_err(|_| {
                ThermalTransportError::Invalid("enthalpy outside original liquid calorics")
            })?;
        self.conductivity.state(
            state.temperature().value(),
            state.isobaric_heat_capacity_j_per_k() / self.caloric.mass_kg,
        )
    }
}

impl ConductivityLaw for MulticomponentLiquidConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.conductivity
            .at_temperature(self.caloric.temperature_k(h)?)
            .map(|(k, _)| k)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

impl ImplicitConductionLaw for MulticomponentLiquidConductionLaw<'_> {
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
        Self::state_at_specific_enthalpy(self, h).map_err(|_| DepositionError::MaterialDomain)
    }
    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.caloric.temperature_k(h)?,
            phase_mass_fraction: None,
        })
    }
}
