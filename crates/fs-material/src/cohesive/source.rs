//! Source-defined cohesive failure at fixed external conditions.
//!
//! Interface definitions and their source receipts remain immutable. The
//! existing cohesive point owns separation and irreversible damage. The explicit
//! `at_temperature` path retains the same envelope and claims; its successor can
//! be consumed only with the point owner's temperature-history transfer.

/// Explicit normal/shear source data for the frozen-onset BK extension.
pub mod mixed_mode;
pub use mixed_mode::*;

/// Original ordered-interface source binding for physical bond viscosity.
pub mod rate;
pub use rate::*;

use super::BilinearCohesive;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, interface_stiffness_per_area_dims,
    resolve_interface_state_point,
};
use fs_matdb::{InterfaceSystemCard, MaterialStateId, QueryPoint, SurfaceSpec};
use fs_qty::{Dims, Pressure, QuantitySpec};

/// Intact normal/shear stiffness per original bonded area, Pa/m.
pub const COHESIVE_STIFFNESS_PROPERTY: &str = "cohesive_stiffness";
/// Common peak normal/shear traction, Pa.
pub const COHESIVE_STRENGTH_PROPERTY: &str = "cohesive_strength";
/// Complete fracture work per original bonded area, J/m².
pub const COHESIVE_FRACTURE_ENERGY_PROPERTY: &str = "cohesive_fracture_energy";
/// Undamaged compression penalty per original bonded area, Pa/m.
pub const COHESIVE_COMPRESSION_STIFFNESS_PROPERTY: &str = "cohesive_compression_stiffness";
/// Dimensionless static friction coefficient of the residual fracture surfaces.
pub const RESIDUAL_STATIC_FRICTION_PROPERTY: &str = "cohesive_residual_static_friction";
/// Dimensionless kinetic friction coefficient of the residual fracture surfaces.
pub const RESIDUAL_KINETIC_FRICTION_PROPERTY: &str = "cohesive_residual_kinetic_friction";
/// Energy per area, J/m². This is not bulk energy density or toughness in Pa√m.
pub const FRACTURE_ENERGY_PER_AREA_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);

/// A bilinear cohesive law tied to the complete original ordered interface.
///
/// The standard constructor selects common-mode data; the explicitly named
/// mixed-mode constructor selects independent normal/shear data and BK exponent.
/// This does not infer mixed-mode calibration, adhesion from bulk chemistry,
/// physical bonding from geometric contact, or fracture heat partition.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedBilinearCohesive {
    resolved: ResolvedInterfaceStatePoint,
    surfaces: [SurfaceSpec; 2],
    law: BilinearCohesive,
}

impl ResolvedBilinearCohesive {
    /// Consume an existing multiphysics interface bundle without selecting its
    /// sources again. Extra properties and their receipts stay in the bundle.
    ///
    /// Coefficients may depend on explicitly fixed temperature, composition,
    /// moisture or external pressure. The canonical static-source admission
    /// rejects separation, damage, rate, frequency and other unowned axes,
    /// including applicability-only bounds on otherwise constant coefficients.
    /// The ordinary history update keeps these conditions fixed; an explicit
    /// temperature successor requires its history-transfer owner.
    /// The numerical separation limit is explicit and must include failure.
    pub fn try_from_resolved(
        card: &InterfaceSystemCard,
        resolved: ResolvedInterfaceStatePoint,
        maximum_separation_m: f64,
    ) -> Result<Self, MaterialStatePointError> {
        if resolved.card_identity() != card.content_hash() {
            return Err(invalid(
                "cohesive law requires its original ordered interface card",
            ));
        }
        let mut values = [0.0; 4];
        for (value, (name, dims)) in values.iter_mut().zip(properties()) {
            let property = resolved.property(name).ok_or_else(|| {
                MaterialStatePointError::InvalidRequirement {
                    property: name.into(),
                    reason: "missing or ambiguous cohesive interface coefficient",
                }
            })?;
            if property.requirement().quantity() != QuantitySpec::dimensional(dims) {
                return Err(invalid(
                    "cohesive coefficient has incompatible dimensions or quantity semantics",
                ));
            }
            if !property.is_static_at_fixed_external_conditions() {
                return Err(invalid(
                    "cohesive coefficients require fixed external conditions; evolving separation, damage and rates cannot be frozen",
                ));
            }
            *value = property.value_si();
            if !value.is_finite() || *value <= 0.0 {
                return Err(MaterialStatePointError::OutsideConsumerDomain {
                    property: name.into(),
                    value: *value,
                    admissibility: ScalarAdmissibility::StrictlyPositive,
                });
            }
        }
        let [stiffness, strength, fracture_energy, compression] = values;
        let law = BilinearCohesive::new(
            stiffness, strength, fracture_energy, compression, maximum_separation_m,
        ).map_err(|_| invalid("cohesive source coefficients require representable onset < failure <= separation limit"))?;
        Ok(Self {
            resolved,
            surfaces: [card.surface_a().clone(), card.surface_b().clone()],
            law,
        })
    }

    /// Existing immutable law; borrowing it preserves the lifetime of its source.
    #[must_use]
    pub const fn law(&self) -> &BilinearCohesive {
        &self.law
    }

    /// Complete original interface state, selected claims and usage receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedInterfaceStatePoint {
        &self.resolved
    }

    /// Ordered surface identities and texture frames from the original card.
    #[must_use]
    pub const fn surfaces(&self) -> &[SurfaceSpec; 2] {
        &self.surfaces
    }

    /// Match the material states actually assigned to the receiving surfaces.
    /// Texture orientation and co-location are separate geometry obligations.
    pub fn require_surface_materials(
        &self,
        first: &MaterialStateId,
        second: &MaterialStateId,
    ) -> Result<(), MaterialStatePointError> {
        if &self.surfaces[0].material != first || &self.surfaces[1].material != second {
            return Err(invalid(
                "cohesive interface surface material assignment mismatch",
            ));
        }
        Ok(())
    }

    /// Resolve a new temperature from the SAME card and every original selected
    /// claim. Other coordinates, quantity conventions, model/rate bounds and
    /// extra interface properties stay bound. No source is silently replaced.
    /// This returns immutable coefficients, not an evolved cohesive history.
    /// The receiving point must use `CohesivePoint::transfer_temperature`.
    pub fn at_temperature(
        &self,
        card: &InterfaceSystemCard,
        temperature_axis: &str,
        temperature_k: f64,
    ) -> Result<Self, MaterialStatePointError> {
        self.require_temperature_axis(temperature_axis)?;
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid(
                "cohesive temperature must be finite positive kelvin",
            ));
        }
        let point = self
            .resolved
            .query()
            .with(temperature_axis, temperature_k)
            .map_err(|source| MaterialStatePointError::Query {
                property: format!("cohesive temperature:{temperature_axis}"),
                source,
            })?;
        let resolved = self.resolved.at_query_point(card, &point)?;
        let envelope = if self.law.mixed_mode_parameters().is_some() {
            CohesiveEnvelope::FrozenOnsetMixedMode
        } else {
            CohesiveEnvelope::CommonMode
        };
        match self.law.viscosity() {
            Some(viscosity) => Self::try_from_viscous_resolved(
                card,
                resolved,
                envelope,
                self.law.maximum_separation_m(),
                viscosity.maximum_rate_m_s(),
            ),
            None => match envelope {
                CohesiveEnvelope::CommonMode => {
                    Self::try_from_resolved(card, resolved, self.law.maximum_separation_m())
                }
                CohesiveEnvelope::FrozenOnsetMixedMode => Self::try_from_mixed_mode_resolved(
                    card,
                    resolved,
                    self.law.maximum_separation_m(),
                ),
            },
        }
    }

    /// Admit an already resolved temperature successor for native history
    /// transfer. Required even if its numerical coefficients happen to agree.
    /// The physical source, receiving order, selected claims and all nonthermal
    /// coordinates remain identical. The receiving history retains its original
    /// common-mode or frozen-physical-initiation-mode interpretation.
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
            || self.law.maximum_separation_m() != next.law.maximum_separation_m()
            || self.law.mixed_mode_parameters().is_some()
                != next.law.mixed_mode_parameters().is_some()
            || self.law.viscosity().map(|v| v.maximum_rate_m_s())
                != next.law.viscosity().map(|v| v.maximum_rate_m_s())
        {
            return Err(invalid(
                "cohesive temperature successor changed source or model bounds",
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
                "cohesive temperature successor changed a held coordinate",
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
                "cohesive temperature successor changed original selected claims",
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
                .any(|property| !property.has_only_fixed_external_conditions(Some(axis), false))
        {
            return Err(invalid(
                "cohesive refresh requires one original absolute-temperature axis",
            ));
        }
        Ok(())
    }
}

/// Resolve the four explicit cohesive properties using the shared source owner.
/// No missing coefficient is inferred from either bulk material. The maximum
/// separation remains a caller-supplied applicability/numerical control in m.
pub fn resolve_bilinear_cohesive(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    maximum_separation_m: f64,
) -> Result<ResolvedBilinearCohesive, MaterialStatePointError> {
    let requirements: Result<Vec<_>, _> = properties()
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        })
        .collect();
    let resolved = resolve_interface_state_point(card, point, &requirements?, selection)?;
    ResolvedBilinearCohesive::try_from_resolved(card, resolved, maximum_separation_m)
}

fn properties() -> [(&'static str, Dims); 4] {
    let stiffness = interface_stiffness_per_area_dims();
    [
        (COHESIVE_STIFFNESS_PROPERTY, stiffness),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
        ),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness),
    ]
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
