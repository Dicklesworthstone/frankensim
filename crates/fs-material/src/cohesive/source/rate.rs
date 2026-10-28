//! Physical interface viscosity selected with its original intrinsic bond law.
use super::*;
use crate::cohesive::CohesiveViscosity;

/// Viscous traction per normal separation rate, Pa s/m.
pub const COHESIVE_NORMAL_VISCOSITY_PROPERTY: &str = "cohesive_normal_viscosity";
/// Isotropic viscous traction per shear separation rate, Pa s/m.
pub const COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY: &str = "cohesive_tangential_viscosity";
/// Interface viscosity has one more inverse length than bulk dynamic viscosity.
pub const COHESIVE_VISCOSITY_DIMS: Dims = Dims([-2, 1, -1, 0, 0, 0]);

/// Explicit intrinsic fracture envelope; no fallback is selected from data gaps.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CohesiveEnvelope {
    CommonMode,
    FrozenOnsetMixedMode,
}
impl ResolvedBilinearCohesive {
    /// Retain all original receipts and the ordered interface. Eta is fixed at
    /// the same external conditions as the intrinsic law. Evolving separation,
    /// damage, temperature and rate-conditioned source coefficients cannot be
    /// silently frozen. The explicit rate bound belongs to this local model.
    pub fn try_from_viscous_resolved(
        card: &InterfaceSystemCard,
        resolved: ResolvedInterfaceStatePoint,
        envelope: CohesiveEnvelope,
        maximum_separation_m: f64,
        maximum_rate_m_s: f64,
    ) -> Result<Self, MaterialStatePointError> {
        let mut eta = [0.0; 2];
        for (value, name) in eta.iter_mut().zip([
            COHESIVE_NORMAL_VISCOSITY_PROPERTY,
            COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
        ]) {
            let property = resolved.property(name).ok_or_else(|| {
                MaterialStatePointError::InvalidRequirement {
                    property: name.into(),
                    reason: "missing or ambiguous physical cohesive viscosity",
                }
            })?;
            if property.requirement().quantity()
                != QuantitySpec::dimensional(COHESIVE_VISCOSITY_DIMS)
            {
                return Err(invalid("cohesive viscosity requires Pa s/m, not bulk Pa s"));
            }
            if !property.is_static_at_fixed_external_conditions() {
                return Err(invalid(
                    "cohesive viscosity requires fixed external source conditions",
                ));
            }
            *value = property.value_si();
        }
        let viscosity = CohesiveViscosity::new(eta[0], eta[1], maximum_rate_m_s)
            .map_err(|_| invalid("invalid physical cohesive viscosity or declared rate bound"))?;
        let mut bound = match envelope {
            CohesiveEnvelope::CommonMode => {
                Self::try_from_resolved(card, resolved, maximum_separation_m)
            }
            CohesiveEnvelope::FrozenOnsetMixedMode => {
                Self::try_from_mixed_mode_resolved(card, resolved, maximum_separation_m)
            }
        }?;
        bound.law = bound.law.with_viscosity(viscosity);
        Ok(bound)
    }
}

/// Resolve intrinsic fracture and both viscosities in one original source query.
/// Fracture inputs must exclude the viscous work represented by these eta values;
/// total measured finite-rate fracture work is not automatically intrinsic Gc.
pub fn resolve_viscous_bilinear_cohesive(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    envelope: CohesiveEnvelope,
    maximum_separation_m: f64,
    maximum_rate_m_s: f64,
) -> Result<ResolvedBilinearCohesive, MaterialStatePointError> {
    let mut names = match envelope {
        CohesiveEnvelope::CommonMode => properties().to_vec(),
        CohesiveEnvelope::FrozenOnsetMixedMode => mixed_mode::mixed_properties().to_vec(),
    };
    names.extend([
        (COHESIVE_NORMAL_VISCOSITY_PROPERTY, COHESIVE_VISCOSITY_DIMS),
        (
            COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
            COHESIVE_VISCOSITY_DIMS,
        ),
    ]);
    let requirements: Result<Vec<_>, _> = names
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        })
        .collect();
    let resolved = resolve_interface_state_point(card, point, &requirements?, selection)?;
    ResolvedBilinearCohesive::try_from_viscous_resolved(
        card,
        resolved,
        envelope,
        maximum_separation_m,
        maximum_rate_m_s,
    )
}
