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

    /// Refresh the original residual-surface coefficients at one prescribed
    /// absolute temperature. Every selected claim, extra property, surface and
    /// held coordinate stays bound to the same immutable interface card.
    /// This changes no cohesive or friction history and supplies no heat.
    pub fn at_temperature(
        &self,
        card: &InterfaceSystemCard,
        temperature_axis: &str,
        temperature_k: f64,
    ) -> Result<Self, MaterialStatePointError> {
        self.require_temperature_axis(temperature_axis)?;
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid(
                "residual friction temperature requires positive kelvin",
            ));
        }
        let point = self
            .resolved
            .query()
            .with(temperature_axis, temperature_k)
            .map_err(|source| MaterialStatePointError::Query {
                property: format!("residual friction temperature:{temperature_axis}"),
                source,
            })?;
        Self::try_from_resolved(card, self.resolved.at_query_point(card, &point)?)
    }

    /// Admit a temperature successor without replacing the original source or
    /// selecting different friction data. Numerical equality of mu is not a
    /// substitute for the original claim and ordered-surface identity.
    pub fn require_temperature_successor(
        &self,
        next: &Self,
        temperature_axis: &str,
    ) -> Result<(), MaterialStatePointError> {
        self.require_temperature_axis(temperature_axis)?;
        next.require_temperature_axis(temperature_axis)?;
        if self.resolved.card_identity() != next.resolved.card_identity()
            || self.surfaces != next.surfaces
            || self.resolved.axis_quantities() != next.resolved.axis_quantities()
        {
            return Err(invalid(
                "residual friction temperature successor changed source",
            ));
        }
        let a = self.resolved.query_point();
        let b = next.resolved.query_point();
        if a.len() != b.len()
            || a.iter().zip(b).any(|((name, value), (other, new))| {
                name != other || (name != temperature_axis && value != new)
            })
        {
            return Err(invalid(
                "residual friction temperature successor changed a held coordinate",
            ));
        }
        let a = self.resolved.properties();
        let b = next.resolved.properties();
        if a.len() != b.len()
            || a.iter().zip(b).any(|(old, new)| {
                old.requirement() != new.requirement()
                    || old.answer().receipt.selected != new.answer().receipt.selected
            })
        {
            return Err(invalid(
                "residual friction temperature successor changed selected claims",
            ));
        }
        Ok(())
    }

    fn require_temperature_axis(&self, axis: &str) -> Result<(), MaterialStatePointError> {
        use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
        let quantity = self.resolved.axis_quantities().get(axis).copied();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let legacy = axis == "T"
            && (quantity.is_none()
                || quantity == Some(QuantitySpec::dimensional(fs_qty::Temperature::DIMS)));
        if !(legacy || quantity == Some(absolute))
            || !self
                .resolved
                .state_coordinate(axis)
                .is_some_and(|v| v.is_finite() && v > 0.0)
            || self
                .resolved
                .properties()
                .iter()
                .any(|p| !p.has_only_fixed_external_conditions(Some(axis), false))
        {
            return Err(invalid(
                "residual friction refresh requires one original absolute-temperature axis",
            ));
        }
        Ok(())
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
