//! Thermal junction strength on an explicitly held normal-contact state.
//!
//! Only the two original shear claims follow the common face temperature. The
//! original normal query, elastic coefficients, adhesion energy, branch, area
//! and radius remain with the frozen normal owner. This split supplies no
//! normal-law thermal feedback, thermal expansion or temperature-driven peeling.

use core::fmt;

use fs_material::state_point::{MaterialStatePointError, ResolvedInterfaceStatePoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};
use fs_tribo::adhesive::spinning::contact::{
    AdhesiveWrenchContact, AdhesiveWrenchError, AdhesiveWrenchParameters,
};

use super::super::{KINETIC_JUNCTION_SHEAR_PROPERTY, STATIC_JUNCTION_SHEAR_PROPERTY};
use super::SourceAdhesiveSlidingPatch;

/// Original frozen normal patch and its selected shear claims at one trial T.
///
/// This is an immutable material/geometry proposal, not an accepted thermal
/// state. The enclosing owner must close its actual face temperature, native
/// recoil and exported sliding-plus-spin heat together before publication.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceThermalAdhesivePatch<'a> {
    original: SourceAdhesiveSlidingPatch<'a>,
    source: ResolvedInterfaceStatePoint,
    law: AdhesiveWrenchContact,
}

impl<'a> SourceThermalAdhesivePatch<'a> {
    /// Unchanged ordered sources and actual normal-owner solid-contact disk.
    #[must_use]
    pub const fn original(&self) -> SourceAdhesiveSlidingPatch<'a> {
        self.original
    }

    /// Actual thermal query, original pinned claims and complete usage receipts.
    #[must_use]
    pub const fn friction_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }

    /// Original uniform-disk law at the refreshed static and kinetic strengths.
    #[must_use]
    pub const fn law(&self) -> &AdhesiveWrenchContact {
        &self.law
    }
}

/// A thermal source or native constitutive refusal without any accepted change.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceThermalAdhesivePatchError {
    /// The original pinned query or scalar admissibility failed.
    Material(MaterialStatePointError),
    /// The original uniform-junction law refused the refreshed strengths.
    Wrench(AdhesiveWrenchError),
    /// Temperature ownership or actual solid-contact geometry was not admitted.
    Invalid(&'static str),
}

impl From<MaterialStatePointError> for SourceThermalAdhesivePatchError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}

impl From<AdhesiveWrenchError> for SourceThermalAdhesivePatchError {
    fn from(error: AdhesiveWrenchError) -> Self {
        Self::Wrench(error)
    }
}

impl fmt::Display for SourceThermalAdhesivePatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Material(error) => write!(f, "thermal adhesive source: {error}"),
            Self::Wrench(error) => write!(f, "thermal adhesive junction: {error}"),
            Self::Invalid(reason) => write!(f, "thermal adhesive patch: {reason}"),
        }
    }
}

impl std::error::Error for SourceThermalAdhesivePatchError {}

impl<'a> SourceAdhesiveSlidingPatch<'a> {
    /// Refresh only the original junction strengths at a positive face T [K].
    ///
    /// `temperature_axis` must be an existing static absolute-temperature axis,
    /// or the exact legacy `T` convention. Its original quantity convention and
    /// every other query coordinate stay fixed. Selected claims never switch,
    /// even if another claim becomes eligible; out-of-domain queries refuse.
    /// No cooling, weakening or heating trial mutates either source owner.
    ///
    /// Normal geometry and coefficients are explicitly held at their original
    /// source state, including under zero or tensile resultant load. A zero-area
    /// detached state refuses this solid-face entry; it cannot acquire a thermal
    /// contact, friction force or torque from attraction across a gap. Rebinding
    /// the normal law, contact area, or its temperature is a separate operation.
    pub fn at_temperature(
        &self,
        temperature_axis: &str,
        temperature_k: f64,
    ) -> Result<SourceThermalAdhesivePatch<'a>, SourceThermalAdhesivePatchError> {
        let radius = self.patch().parameters().contact_radius_m;
        if radius <= 0.0 || self.contact_area_m2() <= 0.0 {
            return Err(SourceThermalAdhesivePatchError::Invalid(
                "thermal adhesive junction requires an actual attached solid disk",
            ));
        }
        let original = self.shear_source().resolved();
        let quantity = original.axis_quantities().get(temperature_axis).copied();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let legacy = temperature_axis == "T"
            && (quantity.is_none()
                || quantity == Some(QuantitySpec::dimensional(Temperature::DIMS)));
        if !temperature_k.is_finite()
            || temperature_k <= 0.0
            || !(legacy || quantity == Some(absolute))
            || !original
                .state_coordinate(temperature_axis)
                .is_some_and(|v| v.is_finite() && v > 0.0)
            || original.properties().iter().any(|property| {
                !property.has_only_fixed_external_conditions(Some(temperature_axis), false)
            })
        {
            return Err(SourceThermalAdhesivePatchError::Invalid(
                "junction refresh requires one original absolute-temperature axis",
            ));
        }
        let source = if original.state_coordinate(temperature_axis) == Some(temperature_k) {
            // Retain the original receipt policy as well as its exact values
            // when no source coordinate changed, as the thermal Coulomb owner does.
            original.clone()
        } else {
            let point = original
                .query()
                .with(temperature_axis, temperature_k)
                .map_err(|source| MaterialStatePointError::Query {
                    property: format!("adhesive junction temperature:{temperature_axis}"),
                    source,
                })?;
            original.at_query_point(self.normal_state().elastic().interface().card(), &point)?
        };
        let strength = |name| {
            source.property(name).map(|p| p.value_si()).ok_or(
                SourceThermalAdhesivePatchError::Invalid(
                    "original junction strength is absent or ambiguous",
                ),
            )
        };
        let law = AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
            contact_radius_m: radius,
            static_shear_strength_pa: strength(STATIC_JUNCTION_SHEAR_PROPERTY)?,
            kinetic_shear_strength_pa: strength(KINETIC_JUNCTION_SHEAR_PROPERTY)?,
        })?;
        Ok(SourceThermalAdhesivePatch {
            original: *self,
            source,
            law,
        })
    }
}
