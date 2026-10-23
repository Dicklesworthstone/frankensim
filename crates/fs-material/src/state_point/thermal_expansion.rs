//! Exact orthotropic expansion along an admitted, fixed-context source path.
//!
//! The three coefficients are linear strain per kelvin in one fixed principal
//! frame. They are not endpoint strains or a general six-component tensor law.
//! A consumer must explicitly bind that principal frame to its geometry frame.

use super::*;

/// Principal linear thermal-expansion coefficients [1/K], in axis order 1, 2, 3.
/// The material axes must remain fixed over the complete temperature path.
pub const ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES: [&str; 3] = [
    "linear_thermal_expansion_coefficient_1",
    "linear_thermal_expansion_coefficient_2",
    "linear_thermal_expansion_coefficient_3",
];

/// Three source-integrated free linear strains in a fixed material frame.
///
/// All endpoint values and every path segment use the exact claims selected at
/// the reference state. Negative coefficients and cooling are supported. This
/// is the additive small-strain integral, not logarithmic or finite strain.
#[derive(Clone, Debug, PartialEq)]
pub struct IntegratedOrthotropicThermalExpansion {
    reference: ResolvedMaterialStatePoint,
    current: ResolvedMaterialStatePoint,
    selected_claims: [ClaimId; 3],
    free_principal_strains: [f64; 3],
    identity: ContentHash,
}

impl IntegratedOrthotropicThermalExpansion {
    /// Complete selected expansion sources at the reference temperature.
    #[must_use]
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }

    /// The same selected expansion sources at the current temperature.
    #[must_use]
    pub const fn current(&self) -> &ResolvedMaterialStatePoint {
        &self.current
    }

    /// Exact integrated claim in principal-axis order 1, 2, 3.
    #[must_use]
    pub const fn selected_claims(&self) -> &[ClaimId; 3] {
        &self.selected_claims
    }

    /// Signed `integral(alpha_i(T) dT)` in principal-axis order 1, 2, 3.
    #[must_use]
    pub const fn free_principal_strains(&self) -> &[f64; 3] {
        &self.free_principal_strains
    }

    /// Identity binding the immutable card, exact path, selected laws and result.
    /// Orientation is deliberately supplied and identified by the consumer.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Integrate three independently sourced orthotropic coefficients exactly.
///
/// The complete non-temperature query context and its quantity schema must
/// remain fixed. The existing source owner admits every crossed support cell,
/// including fixed-context paths through two-coordinate surfaces. Unsupported
/// gaps, missing axes, ambiguous selection and coefficient schema mismatches
/// refuse the complete result. No scalar fallback or missing-axis completion
/// occurs. Principal axes are fixed along this path; rotating principal frames
/// and general off-diagonal coefficient laws require a different model.
pub fn integrate_orthotropic_thermal_expansion(
    card: &MaterialCard,
    reference_point: &QueryPoint,
    current_point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<IntegratedOrthotropicThermalExpansion, MaterialStatePointError> {
    let (reference_temperature_k, current_temperature_k) =
        matching_temperature_path(reference_point, current_point)?;
    let requirements = ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                INVERSE_TEMPERATURE_DIMS,
                ScalarAdmissibility::Finite,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reference = resolve_material_state_point(card, reference_point, &requirements, selection)?;
    let current = reference.at_query_point(card, current_point)?;
    let selected_claims = ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES.map(|name| {
        reference
            .property(name)
            .expect("all three canonical expansion requirements were resolved")
            .answer()
            .receipt
            .selected
    });
    let mut free_principal_strains = [0.0; 3];
    for (strain, selected) in free_principal_strains.iter_mut().zip(selected_claims) {
        *strain = integrate_selected_thermal_expansion(
            card,
            reference_point,
            current_point,
            selected,
            reference_temperature_k,
            current_temperature_k,
        )?;
    }
    let mut identity =
        DomainHasher::new("org.frankensim.fs-material.integrated-orthotropic-thermal-expansion.v1");
    identity.update(card.content_hash().as_bytes());
    identity.update(reference.identity().as_bytes());
    identity.update(current.identity().as_bytes());
    for (selected, strain) in selected_claims.iter().zip(free_principal_strains) {
        identity.update(selected.0.as_bytes());
        identity.update(&strain.to_bits().to_le_bytes());
    }
    Ok(IntegratedOrthotropicThermalExpansion {
        reference,
        current,
        selected_claims,
        free_principal_strains,
        identity: identity.finalize(),
    })
}
