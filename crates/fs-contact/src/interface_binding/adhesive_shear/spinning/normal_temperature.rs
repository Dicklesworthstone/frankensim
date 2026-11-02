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
use crate::normal_patch::adhesive::{SourceJkrSphereContact, temperature::SourceJkrThermalRefresh};

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
        let normal = refreshed.source().normal_state();
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
