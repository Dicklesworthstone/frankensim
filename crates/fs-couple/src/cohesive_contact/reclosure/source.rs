//! Explicit residual friction data on the original ordered interface.

use fs_matdb::{InterfaceSystemCard, MaterialStateId, QueryPoint, SurfaceSpec};
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::FrictionLaw;

/// Dimensionless residual static coefficient after complete cohesive failure.
pub const RESIDUAL_STATIC_FRICTION_PROPERTY: &str = "cohesive_residual_static_friction";
/// Dimensionless residual kinetic coefficient after complete cohesive failure.
pub const RESIDUAL_KINETIC_FRICTION_PROPERTY: &str = "cohesive_residual_kinetic_friction";

/// Original resolved source retained separately from irreversible bond history.
/// No residual coefficient is inferred from bond strength or fracture energy.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidualCoulombSource {
    resolved: ResolvedInterfaceStatePoint,
    surfaces: [SurfaceSpec; 2],
    law: FrictionLaw,
}
impl ResidualCoulombSource {
    /// Preserve an existing source bundle, including all extra properties and
    /// receipts. Coefficients must be constant at the fixed external query.
    pub fn try_from_resolved(
        card: &InterfaceSystemCard,
        resolved: ResolvedInterfaceStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        if resolved.card_identity() != card.content_hash() {
            return Err(invalid(
                "residual friction requires its original interface card",
            ));
        }
        let mut values = [0.0; 2];
        for (value, name) in values.iter_mut().zip([
            RESIDUAL_STATIC_FRICTION_PROPERTY,
            RESIDUAL_KINETIC_FRICTION_PROPERTY,
        ]) {
            let property = resolved.property(name).ok_or_else(|| {
                MaterialStatePointError::InvalidRequirement {
                    property: name.into(),
                    reason: "missing or ambiguous residual friction coefficient",
                }
            })?;
            if property.requirement().quantity() != QuantitySpec::dimensional(Dims::NONE)
                || !property.is_static_at_fixed_external_conditions()
            {
                return Err(invalid(
                    "residual friction requires dimensionless fixed-source coefficients",
                ));
            }
            *value = property.value_si();
        }
        let law = FrictionLaw::Coulomb {
            static_mu: values[0],
            kinetic_mu: values[1],
        };
        law.validate()
            .map_err(|_| invalid("residual friction requires 0 <= kinetic <= static"))?;
        Ok(Self {
            resolved,
            surfaces: [card.surface_a().clone(), card.surface_b().clone()],
            law,
        })
    }
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedInterfaceStatePoint {
        &self.resolved
    }
    #[must_use]
    pub const fn surfaces(&self) -> &[SurfaceSpec; 2] {
        &self.surfaces
    }
    #[must_use]
    pub const fn law(&self) -> &FrictionLaw {
        &self.law
    }
    pub fn require_surface_materials(
        &self,
        receiving: [&MaterialStateId; 2],
    ) -> Result<(), MaterialStatePointError> {
        if [&self.surfaces[0].material, &self.surfaces[1].material] != receiving {
            return Err(invalid(
                "residual friction receiving surface assignment mismatch",
            ));
        }
        Ok(())
    }
    /// This first handoff retains the original interface card and fixed query.
    /// A changed fracture-surface system needs a separate admitted transfer.
    pub fn require_cohesive_source(
        &self,
        cohesive: &ResolvedBilinearCohesive,
    ) -> Result<(), MaterialStatePointError> {
        if self.resolved.card_identity() != cohesive.resolved().card_identity()
            || self.resolved.query_point() != cohesive.resolved().query_point()
            || self.resolved.axis_quantities() != cohesive.resolved().axis_quantities()
            || &self.surfaces != cohesive.surfaces()
        {
            return Err(invalid(
                "residual friction must retain the original cohesive interface and fixed query",
            ));
        }
        Ok(())
    }
}

/// Select only the explicitly named residual data through the existing material
/// source resolver; the complete original card identity remains attached.
pub fn resolve_residual_coulomb(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<ResidualCoulombSource, MaterialStatePointError> {
    let requirements = [
        RESIDUAL_STATIC_FRICTION_PROPERTY,
        RESIDUAL_KINETIC_FRICTION_PROPERTY,
    ]
    .into_iter()
    .map(|name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
    })
    .collect::<Result<Vec<_>, _>>()?;
    let resolved = resolve_interface_state_point(card, point, &requirements, selection)?;
    ResidualCoulombSource::try_from_resolved(card, resolved)
}
fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
