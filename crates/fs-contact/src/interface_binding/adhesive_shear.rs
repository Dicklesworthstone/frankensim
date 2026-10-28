//! Constant junction shear from the original ordered interface's property claims.
//!
//! The explicit model is `F_static = tau_static A`, `F_sliding = tau_kinetic A`.
//! No Coulomb term is added: this avoids silently counting a measured junction
//! contribution twice. Neither strength is inferred from work of adhesion or
//! bulk yield. Normal geometry remains owned by the normal-contact solver.

use core::fmt;

use fs_blake3::ContentHash;
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_qty::Pressure;
use fs_tribo::adhesive::{AdhesiveContactState, AdhesiveFrictionLaw, AdhesiveFrictionResponse};
use fs_tribo::{FrictionLaw, TangentialSlip, TriboError};

use super::BoundNormalInterfaceState;
use crate::normal_patch::adhesive::JkrResponse;

/// Interface property for the static junction shear strength, Pa.
pub const STATIC_JUNCTION_SHEAR_PROPERTY: &str = "adhesive_static_shear_strength";
/// Interface property for the positive-speed junction shear strength, Pa.
pub const KINETIC_JUNCTION_SHEAR_PROPERTY: &str = "adhesive_kinetic_shear_strength";

/// A constant-shear-strength law with retained source receipts and normal identity.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceAdhesiveShear {
    normal_identity: ContentHash,
    resolved: ResolvedInterfaceStatePoint,
    law: AdhesiveFrictionLaw,
}

/// Explicit refusal of missing data, an unowned source coordinate or identity drift.
#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveShearError {
    /// Original material query or physical scalar admission failed.
    Material(MaterialStatePointError),
    /// The constitutive leaf refused the strengths, geometry or velocity.
    Tribo(TriboError),
    /// Source/model compatibility was not established.
    Invalid {
        /// Stable refusal reason.
        reason: &'static str,
    },
}

impl fmt::Display for AdhesiveShearError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Material(error) => write!(f, "adhesive shear source: {error}"),
            Self::Tribo(error) => write!(f, "adhesive shear law: {error}"),
            Self::Invalid { reason } => write!(f, "adhesive shear: {reason}"),
        }
    }
}

impl std::error::Error for AdhesiveShearError {}

impl From<MaterialStatePointError> for AdhesiveShearError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}

impl From<TriboError> for AdhesiveShearError {
    fn from(error: TriboError) -> Self {
        Self::Tribo(error)
    }
}

impl SourceAdhesiveShear {
    /// Resolve both strengths at the EXACT ordered card and query of the normal
    /// owner. Temperature and other explicitly fixed external conditions may be
    /// supported; unowned rate, strain, gap and history restrictions refuse.
    /// Calling this constructor explicitly selects the constant-junction model.
    pub fn resolve(
        normal: &BoundNormalInterfaceState,
        selection: MaterialPropertySelection,
    ) -> Result<Self, AdhesiveShearError> {
        let requirements = [
            ScalarPropertyRequirement::try_new(
                STATIC_JUNCTION_SHEAR_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::NonNegative,
            )?,
            ScalarPropertyRequirement::try_new(
                KINETIC_JUNCTION_SHEAR_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::NonNegative,
            )?,
        ];
        let resolved = resolve_interface_state_point(
            normal.elastic().interface().card(),
            &normal.interface_state().query(),
            &requirements,
            selection,
        )?;
        if resolved
            .properties()
            .iter()
            .any(|property| !property.is_static_at_fixed_external_conditions())
        {
            return Err(AdhesiveShearError::Invalid {
                reason: "junction strengths require fixed external source conditions",
            });
        }
        let strength = |name| {
            resolved
                .property(name)
                .map(|property| property.value_si())
                .ok_or(AdhesiveShearError::Invalid {
                    reason: "junction strength is absent or ambiguous",
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

    /// Original selected property values, conditions, uncertainty and usage receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedInterfaceStatePoint {
        &self.resolved
    }

    /// Exact normal source identity required at evaluation.
    #[must_use]
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_identity
    }

    /// Constitutive leaf selected by this explicit constant-junction model.
    #[must_use]
    pub const fn law(&self) -> &AdhesiveFrictionLaw {
        &self.law
    }

    /// Evaluate supplied contact geometry under the unchanged normal owner.
    /// Provenance is borrowed from that owner without upgrading its authority.
    pub fn evaluate<'a>(
        &self,
        normal: &'a BoundNormalInterfaceState,
        contact: AdhesiveContactState,
        slip: TangentialSlip,
    ) -> Result<AdhesiveFrictionResponse<'a>, AdhesiveShearError> {
        if normal.identity() != self.normal_identity {
            return Err(AdhesiveShearError::Invalid {
                reason: "normal source, ordered interface or state point changed",
            });
        }
        Ok(self
            .law
            .evaluate(normal.elastic().interface().interface(), contact, slip)?)
    }

    /// Use actual JKR area and repulsive force, including a tensile net reaction.
    /// A response from a different ordered source state refuses before evaluation.
    /// Frictional power is returned for the thermal owner to allocate once.
    pub fn evaluate_jkr<'a>(
        &self,
        normal: &'a BoundNormalInterfaceState,
        response: &JkrResponse,
        slip: TangentialSlip,
    ) -> Result<AdhesiveFrictionResponse<'a>, AdhesiveShearError> {
        if response.normal_state_identity() != self.normal_identity {
            return Err(AdhesiveShearError::Invalid {
                reason: "JKR response belongs to a different normal source state",
            });
        }
        self.evaluate(normal, response.friction_contact_state()?, slip)
    }
}
