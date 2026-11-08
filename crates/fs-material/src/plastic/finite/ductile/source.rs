//! Source admission for finite ductile attenuation, without numerical refitting.
//!
//! The work ordinate is explicitly equivalent KIRCHHOFF plastic work per
//! REFERENCE area after initiation: `Gp = integral g R(alpha) h d(alpha)`
//! with fixed reference width h and the selected native linear/Voce resistance.
//! It is distinct from total fracture dissipation, cohesive work, current-area
//! work and the existing small-strain property. Equal units do not alias them.
//! Resolution retains the default exponential shape. An optional explicit
//! builder authors a different normalized work profile while preserving the
//! same capacity receipt; the receipt alone does not validate that shape.
//!
//! The original finite model must be a member of the same immutable material
//! card as the work claim and initiation locus. Complete fixed typed contexts
//! must agree. Trials use the prepared law, never query sources, and borrow
//! it from this owner so its source receipts remain alive with the history.

use super::{FiniteDuctileError, FiniteDuctileJ2};
use crate::MaterialError;
use crate::ductile::evolution::source::{DuctileJ2SourceSettings, PLASTIC_WORK_PER_AREA_DIMS};
use crate::ductile::evolution::{DuctileJ2Settings, DuctileWorkEvolution};
use crate::ductile::source::ResolvedDuctileLocus;
use crate::plastic::finite::source::ResolvedFiniteJ2;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, is_fixed_external_coordinate,
    resolve_material_state_point,
};
use fs_matdb::{ConstitutiveModelCard, MaterialCard, QueryPoint};
use fs_qty::QuantitySpec;

/// This exact finite reference-work meaning cannot be supplied by other Gc/Gp keys.
pub const FINITE_POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY: &str =
    "finite_ductile_post_initiation_equivalent_kirchhoff_plastic_work_per_reference_area";

/// One explicit finite-work claim; other physics may add its own requirements.
pub fn finite_ductile_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    Ok(vec![ScalarPropertyRequirement::try_new(
        FINITE_POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY,
        PLASTIC_WORK_PER_AREA_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?])
}

/// Original source owners and the prepared finite ductile law.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedFiniteDuctileJ2 {
    material: ResolvedMaterialStatePoint,
    plastic: ResolvedFiniteJ2,
    law: FiniteDuctileJ2<ResolvedDuctileLocus>,
}

impl ResolvedFiniteDuctileJ2 {
    /// Adopt already selected sources. The complete original card proves actual
    /// model membership; matching numerical coefficients or names is insufficient.
    pub fn try_from_finite_j2(
        card: &MaterialCard,
        material: ResolvedMaterialStatePoint,
        plastic: ResolvedFiniteJ2,
        initiation: ResolvedDuctileLocus,
        settings: DuctileJ2SourceSettings,
    ) -> Result<Self, FiniteDuctileSourceError> {
        let point = material.query();
        if card.content_hash() != material.card_identity()
            || card.id() != material.material()
            || !card
                .models()
                .iter()
                .any(|model| model == plastic.model_card())
            || plastic.context() != &point
            || !initiation.source().matches_material(&material)
            || point.axes().iter().any(|(axis, value)| {
                !is_fixed_external_coordinate(
                    axis,
                    point.axis_quantities().get(axis).copied(),
                    Some(*value),
                )
            })
        {
            return Err(FiniteDuctileSourceError::Invalid(
                "finite plasticity, reference-work and initiation require the same original card, model membership and complete fixed typed context",
            ));
        }
        let property = material
            .property(FINITE_POST_INITIATION_EQUIVALENT_PLASTIC_WORK_PROPERTY)
            .ok_or(FiniteDuctileSourceError::Invalid(
                "missing or ambiguous finite Kirchhoff plastic-work per reference-area calibration",
            ))?;
        let requirement = property.requirement();
        if requirement.quantity() != QuantitySpec::dimensional(PLASTIC_WORK_PER_AREA_DIMS)
            || requirement.hardness_test().is_some()
            || requirement.elastic_component().is_some()
            || requirement.strain_component().is_some()
            || requirement.stress_component().is_some()
            || !property.is_static_at_fixed_external_conditions()
            || !property.value_si().is_finite()
            || property.value_si() <= 0.0
        {
            return Err(FiniteDuctileSourceError::Invalid(
                "finite reference-work requires its original positive scalar J/m2 claim at fixed external conditions",
            ));
        }
        let law = FiniteDuctileJ2::new(
            plastic.law().clone(),
            initiation,
            DuctileJ2Settings {
                characteristic_length_m: settings.characteristic_length_m,
                post_initiation_plastic_work_j_m2: property.value_si(),
                maximum_damage_increment: settings.maximum_damage_increment,
                energy_tolerances: settings.energy_tolerances,
            },
        )?;
        Ok(Self {
            material,
            plastic,
            law,
        })
    }

    /// Explicitly author a normalized post-initiation shape before borrowing
    /// histories. The original equivalent Kirchhoff-work CAPACITY, initiation,
    /// native finite model, typed context and source receipts remain unchanged.
    /// A scalar capacity claim does not establish measured support for a new
    /// shape, total fracture energy, crack separation, or material removal.
    pub fn with_evolution(mut self, evolution: DuctileWorkEvolution) -> Self {
        self.law = self.law.with_evolution(evolution);
        self
    }

    /// Exact work claim, immutable card identity and original property-use receipt.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Original finite model, including its Kirchhoff convention, Voce and domain.
    pub const fn plastic(&self) -> &ResolvedFiniteJ2 {
        &self.plastic
    }
    /// Borrowed prepared law with original initiation support and source receipts.
    pub const fn law(&self) -> &FiniteDuctileJ2<ResolvedDuctileLocus> {
        &self.law
    }
}

/// Resolve an explicitly selected original finite model and reference-work claim.
/// No small-strain coefficient, current-area work or inferred source is promoted.
pub fn resolve_finite_ductile_j2(
    card: &MaterialCard,
    model: &ConstitutiveModelCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    initiation: ResolvedDuctileLocus,
    settings: DuctileJ2SourceSettings,
) -> Result<ResolvedFiniteDuctileJ2, FiniteDuctileSourceError> {
    let plastic = ResolvedFiniteJ2::try_from_model_card_at(model, point)?;
    let material =
        resolve_material_state_point(card, point, &finite_ductile_requirements()?, selection)?;
    ResolvedFiniteDuctileJ2::try_from_finite_j2(card, material, plastic, initiation, settings)
}

/// Original source, finite constitutive or evolution refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteDuctileSourceError {
    MaterialState(MaterialStatePointError),
    Plastic(MaterialError),
    Native(FiniteDuctileError),
    Invalid(&'static str),
}
impl core::fmt::Display for FiniteDuctileSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ductile source refused: {self:?}")
    }
}
impl std::error::Error for FiniteDuctileSourceError {}
impl From<MaterialStatePointError> for FiniteDuctileSourceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::MaterialState(error)
    }
}
impl From<MaterialError> for FiniteDuctileSourceError {
    fn from(error: MaterialError) -> Self {
        Self::Plastic(error)
    }
}
impl From<FiniteDuctileError> for FiniteDuctileSourceError {
    fn from(error: FiniteDuctileError) -> Self {
        Self::Native(error)
    }
}
