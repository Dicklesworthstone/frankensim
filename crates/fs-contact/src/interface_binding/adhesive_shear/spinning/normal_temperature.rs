//! Shear-source continuation through an admitted thermal normal-contact update.
//!
//! The normal owner supplies the new elastic state, work of adhesion and actual
//! solid disk. This operation refreshes the SAME selected shear claims at that
//! owner's interface temperature, preserving independent bulk temperatures.

use fs_tribo::{FrictionLaw, adhesive::AdhesiveFrictionLaw};

use super::super::{
    AdhesiveShearError, KINETIC_JUNCTION_SHEAR_PROPERTY, STATIC_JUNCTION_SHEAR_PROPERTY,
    SourceAdhesiveShear,
};
use crate::interface_binding::BoundNormalInterfaceState;
use crate::normal_patch::adhesive::{
    SourceJkrSphereContact,
    temperature::{SourceDmtThermalRefresh, SourceJkrThermalRefresh, SourceMaugisThermalRefresh},
};
use crate::normal_patch::{dmt::SourceDmtSphereContact, maugis::SourceMaugisSphereContact};

impl SourceAdhesiveShear {
    /// Continue original junction-strength claims on an admitted JKR successor.
    ///
    /// This differs from the frozen-normal `SourceAdhesiveSlidingPatch` thermal
    /// entry: the returned owner belongs to the new normal state. Its sliding
    /// patch MUST be made from the successor's actual normal response, so both
    /// shear force and torsional capacity follow the changed solid-contact area.
    /// The immutable normal refresh proves original model/source continuity.
    /// Missing support, a foreign predecessor or invalid strengths refuses.
    ///
    /// No thermal energy is inferred from the change in strength. A uniform
    /// junction shear law has no stored elastic tangential energy; actual slip
    /// and spin work remains the responsibility of its mechanical integrator.
    pub fn rebind_jkr_temperature(
        &self,
        original: &SourceJkrSphereContact,
        refreshed: &SourceJkrThermalRefresh,
    ) -> Result<Self, AdhesiveShearError> {
        if self.normal_identity != original.normal_state().identity()
            || refreshed.original_model_identity() != original.identity()
            || refreshed.before().normal_state_identity() != self.normal_identity
        {
            return Err(AdhesiveShearError::Invalid {
                reason: "thermal junction continuation requires its original JKR and shear owners",
            });
        }
        self.rebind_thermal_normal(refreshed.source().normal_state())
    }

    /// Continue the same junction claims on the original thermally updated DMT
    /// law. Only its actual solid disk can carry shear; its attractive gap has
    /// zero solid area regardless of the refreshed work of adhesion.
    pub fn rebind_dmt_temperature(
        &self,
        original: &SourceDmtSphereContact,
        refreshed: &SourceDmtThermalRefresh,
    ) -> Result<Self, AdhesiveShearError> {
        if self.normal_identity != original.normal_state().identity()
            || refreshed.original_model_identity() != original.identity()
            || refreshed.before().normal_state_identity() != self.normal_identity
        {
            return Err(AdhesiveShearError::Invalid {
                reason: "thermal junction continuation requires its original DMT and shear owners",
            });
        }
        self.rebind_thermal_normal(refreshed.source().normal_state())
    }

    /// Continue the original Maugis junction claims while retaining its complete
    /// contacting/noncontacting phase. The cohesive annulus does not become a
    /// shear-bearing disk when its temperature or attraction changes.
    pub fn rebind_maugis_temperature(
        &self,
        original: &SourceMaugisSphereContact,
        refreshed: &SourceMaugisThermalRefresh,
    ) -> Result<Self, AdhesiveShearError> {
        if self.normal_identity != original.normal_state().identity()
            || refreshed.original_model_identity() != original.identity()
            || refreshed.before().model_identity() != original.identity()
        {
            return Err(AdhesiveShearError::Invalid {
                reason: "thermal junction continuation requires its original Maugis and shear owners",
            });
        }
        self.rebind_thermal_normal(refreshed.source().normal_state())
    }

    // Public entries first establish the immutable original-model successor;
    // this shared leaf resolves only the original pinned junction claims.
    fn rebind_thermal_normal(
        &self,
        normal: &BoundNormalInterfaceState,
    ) -> Result<Self, AdhesiveShearError> {
        let point = normal.interface_state();
        let resolved = if self.resolved.query_point() == point.query_point()
            && self.resolved.axis_quantities() == point.axis_quantities()
        {
            self.resolved.clone()
        } else {
            self.resolved
                .at_query_point(normal.elastic().interface().card(), &point.query())?
        };
        if resolved
            .properties()
            .iter()
            .any(|property| !property.has_only_fixed_external_conditions(Some("T"), false))
        {
            return Err(AdhesiveShearError::Invalid {
                reason: "thermal junction strengths require owned T and fixed other source conditions",
            });
        }
        let strength = |name| {
            resolved
                .property(name)
                .map(|property| property.value_si())
                .ok_or(AdhesiveShearError::Invalid {
                    reason: "original thermal junction strength is absent or ambiguous",
                })
        };
        let law = AdhesiveFrictionLaw::new(
            FrictionLaw::Coulomb {
                static_mu: 0.0,
                kinetic_mu: 0.0,
            },
            strength(STATIC_JUNCTION_SHEAR_PROPERTY)?,
            strength(KINETIC_JUNCTION_SHEAR_PROPERTY)?,
        )?;
        Ok(Self {
            normal_identity: normal.identity(),
            resolved,
            law,
        })
    }
}
