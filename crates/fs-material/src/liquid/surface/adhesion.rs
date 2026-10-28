//! Thermodynamic adhesion data feeding the existing native equilibrium wall.
//!
//! The model-predicted angle retains the original work and liquid/gas tension
//! receipts. It is never inserted into the material database as a measured
//! contact-angle claim. The caller supplies the same physical vapor environment.

use super::{
    ResolvedLiquidSurfaceState, ResolvedLiquidWettingState, SURFACE_TENSION_DIMS,
    SURFACE_TENSION_PROPERTY, require_frozen_coefficient,
};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ResolvedScalarProperty, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point,
};
use fs_blake3::ContentHash;
use fs_matdb::{InterfaceSystemCard, QueryPoint};

mod law;
pub use law::{YoungDupreError, YoungDupreLaw};

/// Reversible equilibrium solid/liquid separation work in the declared vapor
/// environment, J/m2. This is not dry fracture toughness or dissipative peel work.
pub const EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY: &str = "equilibrium-work-of-adhesion";

/// Whether an equilibrium boundary angle is observed or explicitly modeled.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EquilibriumAngleOrigin {
    /// Original equilibrium-contact-angle source claim.
    SourcedAngle,
    /// Young--Dupre prediction from original adhesion and tension claims.
    YoungDupre,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct YoungDupreSource {
    law: YoungDupreLaw,
    tension: ResolvedScalarProperty,
    liquid_card: ContentHash,
}
impl YoungDupreSource {
    pub(super) fn require_liquid(
        &self,
        liquid: &ResolvedLiquidSurfaceState,
    ) -> Result<(), MaterialStatePointError> {
        if liquid.material().card_identity() != self.liquid_card
            || liquid.material().property(SURFACE_TENSION_PROPERTY) != Some(&self.tension)
            || liquid.surface_tension_n_m() != Some(self.law.surface_tension_n_m())
        {
            return Err(invalid(
                "Young-Dupre boundary requires its original liquid/gas tension source and state",
            ));
        }
        Ok(())
    }
}

impl ResolvedLiquidWettingState {
    /// Explicitly select ideal partial-wetting equilibrium from reversible work
    /// and the actual liquid/gas tension. The supplied interface bundle retains
    /// all its original properties and receipts; no angle claim is fabricated.
    pub fn try_from_adhesion_work(
        card: &InterfaceSystemCard,
        resolved: ResolvedInterfaceStatePoint,
        liquid: &ResolvedLiquidSurfaceState,
    ) -> Result<Self, MaterialStatePointError> {
        if resolved.card_identity() != card.content_hash() {
            return Err(invalid(
                "adhesion wetting requires the original ordered interface card",
            ));
        }
        let work = resolved
            .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
            .ok_or_else(|| invalid("missing or ambiguous reversible equilibrium adhesion work"))?;
        require_frozen_coefficient(
            work,
            EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
            SURFACE_TENSION_DIMS,
        )?;
        let sigma = liquid.surface_tension_n_m().ok_or_else(|| {
            invalid("Young-Dupre wetting requires active sourced liquid/gas tension")
        })?;
        let tension = liquid
            .material()
            .property(SURFACE_TENSION_PROPERTY)
            .ok_or_else(|| invalid("missing original liquid/gas surface-tension receipt"))?;
        require_frozen_coefficient(tension, SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS)?;
        let law = YoungDupreLaw::new(work.value_si(), sigma).map_err(|_| {
            invalid(
                "Young-Dupre data have no representable finite-angle partial-wetting equilibrium",
            )
        })?;
        let result = Self {
            resolved,
            surfaces: [card.surface_a().clone(), card.surface_b().clone()],
            context: card.context().clone(),
            contact_angle_radians: law.contact_angle_radians(),
            young_dupre: Some(YoungDupreSource {
                law,
                tension: tension.clone(),
                liquid_card: liquid.material().card_identity(),
            }),
        };
        result.require_liquid_state(liquid)?;
        Ok(result)
    }

    /// Provenance of the angle value without confusing a prediction with a claim.
    #[must_use]
    pub const fn angle_origin(&self) -> EquilibriumAngleOrigin {
        if self.young_dupre.is_some() {
            EquilibriumAngleOrigin::YoungDupre
        } else {
            EquilibriumAngleOrigin::SourcedAngle
        }
    }

    /// Selected thermodynamic model, when this boundary derives from adhesion.
    #[must_use]
    pub fn young_dupre_law(&self) -> Option<YoungDupreLaw> {
        self.young_dupre.as_ref().map(|source| source.law)
    }

    /// Original liquid/gas tension receipt behind a Young--Dupre prediction.
    #[must_use]
    pub fn surface_tension_source(&self) -> Option<&ResolvedScalarProperty> {
        self.young_dupre.as_ref().map(|source| &source.tension)
    }

    /// Refresh the original interface AND a supplied already-refreshed liquid
    /// at the same complete point. Original card and selected tension claim
    /// cannot change. The old boundary remains untouched on any refusal.
    /// This constructs a new frozen endpoint, not a thermal energy evolution.
    pub fn at_query_point_with_liquid(
        &self,
        card: &InterfaceSystemCard,
        point: &QueryPoint,
        liquid: &ResolvedLiquidSurfaceState,
    ) -> Result<Self, MaterialStatePointError> {
        let Some(previous) = &self.young_dupre else {
            let refreshed = self.at_query_point(card, point)?;
            refreshed.require_liquid_state(liquid)?;
            return Ok(refreshed);
        };
        let tension = liquid
            .material()
            .property(SURFACE_TENSION_PROPERTY)
            .ok_or_else(|| invalid("Young-Dupre refresh lost the original tension claim"))?;
        if liquid.material().card_identity() != previous.liquid_card
            || tension.answer().receipt.selected != previous.tension.answer().receipt.selected
        {
            return Err(invalid(
                "Young-Dupre refresh cannot substitute liquid/gas tension sources",
            ));
        }
        let resolved = self.resolved.at_query_point(card, point)?;
        Self::try_from_adhesion_work(card, resolved, liquid)
    }
}

/// Resolve reversible adhesion at the liquid's exact physical query, then
/// return the EXISTING equilibrium-boundary type consumed by native LBM.
/// Smooth, homogeneous, rigid partial-wetting equilibrium is an explicit model
/// choice; no hysteresis, surfactant kinetics, film spreading or peel strength
/// is inferred. Both liquid and substrate sources retain their original context.
pub fn resolve_young_dupre_wetting(
    card: &InterfaceSystemCard,
    liquid: &ResolvedLiquidSurfaceState,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidWettingState, MaterialStatePointError> {
    let requirement = ScalarPropertyRequirement::try_new(
        EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
        SURFACE_TENSION_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?;
    let resolved =
        resolve_interface_state_point(card, &liquid.material().query(), &[requirement], selection)?;
    ResolvedLiquidWettingState::try_from_adhesion_work(card, resolved, liquid)
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
