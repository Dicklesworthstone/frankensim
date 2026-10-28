//! Source binding for the explicitly selected frozen-onset mixed-mode law.
use super::*;
use crate::cohesive::mixed_mode::MixedModeParameters;

/// Pure opening stiffness, Pa/m.
pub const COHESIVE_NORMAL_STIFFNESS_PROPERTY: &str = "cohesive_normal_stiffness";
/// Stiffness in each shear direction, Pa/m.
pub const COHESIVE_TANGENTIAL_STIFFNESS_PROPERTY: &str = "cohesive_tangential_stiffness";
/// Pure opening strength, Pa.
pub const COHESIVE_NORMAL_STRENGTH_PROPERTY: &str = "cohesive_normal_strength";
/// Pure shear strength, Pa.
pub const COHESIVE_TANGENTIAL_STRENGTH_PROPERTY: &str = "cohesive_tangential_strength";
/// Pure opening fracture work, J/m2.
pub const COHESIVE_MODE_I_FRACTURE_ENERGY_PROPERTY: &str = "cohesive_mode_i_fracture_energy";
/// Pure shear fracture work, J/m2.
pub const COHESIVE_MODE_II_FRACTURE_ENERGY_PROPERTY: &str = "cohesive_mode_ii_fracture_energy";
/// Dimensionless calibrated BK interaction exponent.
pub const COHESIVE_BK_EXPONENT_PROPERTY: &str = "cohesive_bk_exponent";

impl ResolvedBilinearCohesive {
    /// Select mode-dependent fracture from explicit original interface claims.
    /// Missing normal/shear data do not fall back to common-mode or bulk data.
    /// The existing solid, matching/nonmatching P1, and normal thermal-joint
    /// owners consume this same immutable law and preserve its source bundle.
    pub fn try_from_mixed_mode_resolved(
        card: &InterfaceSystemCard,
        resolved: ResolvedInterfaceStatePoint,
        maximum_separation_m: f64,
    ) -> Result<Self, MaterialStatePointError> {
        if resolved.card_identity() != card.content_hash() {
            return Err(invalid(
                "mixed cohesive law requires original ordered interface card",
            ));
        }
        let mut values = [0.0; 8];
        for (value, (name, dims)) in values.iter_mut().zip(mixed_properties()) {
            let property = resolved.property(name).ok_or_else(|| {
                MaterialStatePointError::InvalidRequirement {
                    property: name.into(),
                    reason: "missing or ambiguous mixed cohesive coefficient",
                }
            })?;
            if property.requirement().quantity() != QuantitySpec::dimensional(dims) {
                return Err(invalid(
                    "mixed cohesive coefficient has incompatible quantity",
                ));
            }
            if !property.is_static_at_fixed_external_conditions() {
                return Err(invalid(
                    "mixed cohesive coefficients require fixed external conditions",
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
        let [kn, kt, tn, tt, gi, gii, eta, compression] = values;
        let law = BilinearCohesive::new_mixed_mode(
            MixedModeParameters {
                normal_stiffness_pa_m: kn,
                tangential_stiffness_pa_m: kt,
                normal_strength_pa: tn,
                tangential_strength_pa: tt,
                mode_i_fracture_energy_j_m2: gi,
                mode_ii_fracture_energy_j_m2: gii,
                bk_exponent: eta,
            },
            compression,
            maximum_separation_m,
        )
        .map_err(|_| {
            invalid("mixed cohesive coefficients violate positive-softening or separation bounds")
        })?;
        Ok(Self {
            resolved,
            surfaces: [card.surface_a().clone(), card.surface_b().clone()],
            law,
        })
    }
}

/// Resolve all eight mode-specific coefficients through the original query owner.
/// Calling this function explicitly selects the frozen-onset BK extension.
pub fn resolve_mixed_mode_bilinear_cohesive(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    maximum_separation_m: f64,
) -> Result<ResolvedBilinearCohesive, MaterialStatePointError> {
    let requirements: Result<Vec<_>, _> = mixed_properties()
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        })
        .collect();
    let resolved = resolve_interface_state_point(card, point, &requirements?, selection)?;
    ResolvedBilinearCohesive::try_from_mixed_mode_resolved(card, resolved, maximum_separation_m)
}

pub(super) fn mixed_properties() -> [(&'static str, Dims); 8] {
    let stiffness = interface_stiffness_per_area_dims();
    [
        (COHESIVE_NORMAL_STIFFNESS_PROPERTY, stiffness),
        (COHESIVE_TANGENTIAL_STIFFNESS_PROPERTY, stiffness),
        (COHESIVE_NORMAL_STRENGTH_PROPERTY, Pressure::DIMS),
        (COHESIVE_TANGENTIAL_STRENGTH_PROPERTY, Pressure::DIMS),
        (
            COHESIVE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
        ),
        (
            COHESIVE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
        ),
        (COHESIVE_BK_EXPONENT_PROPERTY, Dims::NONE),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness),
    ]
}
