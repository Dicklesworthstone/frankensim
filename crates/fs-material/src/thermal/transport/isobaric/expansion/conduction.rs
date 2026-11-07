//! Source-exact isobaric calorics, conductivity and free-expansion density.
//!
//! This composition borrows the original Cp/k owner and its already qualified
//! free-expansion path. Native storage remains m*h: Cp already includes the
//! maintained-pressure boundary work, so no separate p*dV heat term is added.

use super::IsobaricFreeExpansion;
use crate::thermal::transport::{
    ThermalTransportError, isobaric::conduction::IsobaricConductionLaw,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::PropertyLineSlice;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// The same source isobar with an explicit original density closure.
///
/// Register `caloric()` on the thermal timeline, and this object as its
/// conductivity/implicit law. The native planar material-density consumer may
/// derive half-thickness from conserved mass and its declared fixed footprint.
/// That reduction does not establish 3-D strain compatibility: a fixed footprint
/// must be physically applicable to the declared expansion path. Tangential
/// expansion, evolving area, stress and radial motion require their own owners.
#[derive(Clone, Debug)]
pub struct IsobaricExpansionConductionLaw<'a> {
    expansion: &'a IsobaricFreeExpansion<'a>,
    conduction: &'a IsobaricConductionLaw<'a>,
    identity: ContentHash,
}

impl<'a> IsobaricExpansionConductionLaw<'a> {
    /// Compose the exact SAME Cp object; equal values or a matching material
    /// name do not replace source identity, pressure, composition or datum.
    /// Preparation borrows existing tables and makes no new material query.
    pub fn try_new(
        expansion: &'a IsobaricFreeExpansion<'a>,
        conduction: &'a IsobaricConductionLaw<'a>,
    ) -> Result<Self, ThermalTransportError> {
        if !core::ptr::eq(expansion.caloric(), conduction.caloric()) {
            return Err(ThermalTransportError::Invalid(
                "expanding conduction requires the same original isobaric caloric object",
            ));
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-material.isobaric-expansion-conduction.v1");
        hash.update(expansion.identity().as_bytes());
        hash.update(conduction.identity().as_bytes());
        Ok(Self {
            expansion,
            conduction,
            identity: hash.finalize(),
        })
    }

    /// Exact caloric/density object registered on the native timeline.
    #[must_use]
    pub const fn caloric(&self) -> &'a IsobaricFreeExpansion<'a> {
        self.expansion
    }

    /// Original Cp/k preparation, retaining its original sources and knots.
    #[must_use]
    pub const fn original_conduction(&self) -> &'a IsobaricConductionLaw<'a> {
        self.conduction
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        self.conduction.source_line()
    }

    pub(in crate::thermal::transport) fn integrated_conductivity_curve(
        &self,
    ) -> Result<crate::thermal::conduction::ConductivityCurve, DepositionError> {
        self.conduction.integrated_conductivity_curve()
    }

    /// Original inverse, Cp(T), k(T) and their exact enthalpy derivatives.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        self.conduction.state_at_specific_enthalpy(h)
    }
}

impl ConductivityLaw for IsobaricExpansionConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.conduction.conductivity_w_m_k(h)
    }
}

impl ImplicitConductionLaw for IsobaricExpansionConductionLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.expansion
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.conduction.enthalpy_bounds_j_kg()
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Self::state_at_specific_enthalpy(self, h).map_err(|_| DepositionError::MaterialDomain)
    }

    fn specific_volume_derivative_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.expansion
            .specific_volume_derivative_at_specific_enthalpy(h)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        self.conduction.time_observables(h)
    }
}
