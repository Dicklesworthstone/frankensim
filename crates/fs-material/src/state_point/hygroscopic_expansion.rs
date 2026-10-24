//! Source-integrated small-strain swelling on a declared dry-mass moisture axis.
//! Coefficients are derivatives with respect to kg water/kg dry solid, measured
//! relative to one explicit stress-free reference length. They are neither
//! thermal coefficients nor logarithmic/finite-strain swelling stretches.

use super::*;
use fs_qty::semantic::MoistureBasis;

/// Principal `d(epsilon_free)/dM` in fixed material directions 1, 2, 3.
pub const ORTHOTROPIC_HYGROSCOPIC_EXPANSION_PROPERTIES: [&str; 3] = [
    "linear_hygroscopic_expansion_coefficient_1",
    "linear_hygroscopic_expansion_coefficient_2",
    "linear_hygroscopic_expansion_coefficient_3",
];
/// Dry-basis moisture content of the source law's stress-free reference length.
pub const HYGROSCOPIC_REFERENCE_MOISTURE_PROPERTY: &str =
    "hygroscopic_strain_reference_moisture_content";

/// Exact integrated engineering strains and complete source state at both ends.
#[derive(Clone, Debug, PartialEq)]
pub struct IntegratedOrthotropicHygroscopicExpansion {
    reference: ResolvedMaterialStatePoint,
    current: ResolvedMaterialStatePoint,
    free_principal_strains: [f64; 3],
    identity: ContentHash,
}
impl IntegratedOrthotropicHygroscopicExpansion {
    /// Original datum and selected derivative claims.
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }
    /// Same claims re-admitted at the actual current moisture and environment.
    pub const fn current(&self) -> &ResolvedMaterialStatePoint {
        &self.current
    }
    /// Signed additive strain integral, in fixed principal material directions.
    pub const fn free_principal_strains(&self) -> &[f64; 3] {
        &self.free_principal_strains
    }
    /// Binds card, datum, query context, selected source path and result.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Integrate all three hygroscopic coefficients on one continuously supported
/// moisture path. Negative coefficients and drying are admitted. The source's
/// reference-moisture property must equal the reference query; no arbitrary
/// change of the denominator of engineering strain is allowed.
///
/// Both queries explicitly type `moisture_content` as dry-mass moisture and
/// `temperature` as absolute temperature. All other coordinates and the entire
/// quantity schema are fixed. Source selection is pinned at the reference;
/// exact piecewise-linear integration uses the existing material-line owner,
/// including fixed-context cuts through surfaces and all intervening support.
/// Hysteresis, stress-assisted diffusion and sorption heat are separate laws.
pub fn integrate_orthotropic_hygroscopic_expansion(
    card: &MaterialCard,
    reference_point: &QueryPoint,
    current_point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<IntegratedOrthotropicHygroscopicExpansion, MaterialStatePointError> {
    let invalid = |quantity| MaterialStatePointError::InvalidDerived { quantity };
    let dry = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass),
        ValueForm::Static,
    ));
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if reference_point.axis_quantities() != current_point.axis_quantities()
        || reference_point.axes().len() != current_point.axes().len()
        || reference_point.axis_quantities().get("moisture_content") != Some(&dry)
        || reference_point.axis_quantities().get("temperature") != Some(&absolute)
        || reference_point.axes().contains_key("T")
        || !reference_point
            .axes()
            .get("temperature")
            .is_some_and(|t| *t > 0.0)
    {
        return Err(invalid(
            "hygroscopic path needs identical schemas, absolute temperature and dry-mass moisture",
        ));
    }
    for (axis, value) in reference_point.axes() {
        let next = current_point
            .axes()
            .get(axis)
            .ok_or_else(|| invalid("hygroscopic path coordinate set"))?;
        if axis != "moisture_content" && value.to_bits() != next.to_bits() {
            return Err(invalid(
                "non-moisture coordinate changed along hygroscopic path",
            ));
        }
    }
    let from = reference_point.axes()["moisture_content"];
    let to = current_point.axes()["moisture_content"];
    if from < 0.0 || to < 0.0 {
        return Err(invalid("negative dry-basis moisture"));
    }
    let mut requirements = ORTHOTROPIC_HYGROSCOPIC_EXPANSION_PROPERTIES
        .into_iter()
        .map(|key| ScalarPropertyRequirement::try_new(key, Dims::NONE, ScalarAdmissibility::Finite))
        .collect::<Result<Vec<_>, _>>()?;
    requirements.push(ScalarPropertyRequirement::try_with_quantity(
        HYGROSCOPIC_REFERENCE_MOISTURE_PROPERTY,
        dry,
        ScalarAdmissibility::NonNegative,
    )?);
    let reference = resolve_material_state_point(card, reference_point, &requirements, selection)?;
    let current = reference.at_query_point(card, current_point)?;
    if reference
        .properties()
        .iter()
        .any(|p| !p.is_static_at_fixed_external_conditions())
    {
        return Err(invalid(
            "hygroscopic strain needs static environmental source laws",
        ));
    }
    let (lower, upper) = if from <= to {
        (reference_point, current_point)
    } else {
        (current_point, reference_point)
    };
    let mut free_principal_strains = [0.0; 3];
    for (index, key) in ORTHOTROPIC_HYGROSCOPIC_EXPANSION_PROPERTIES
        .into_iter()
        .chain(core::iter::once(HYGROSCOPIC_REFERENCE_MOISTURE_PROPERTY))
        .enumerate()
    {
        let property = reference
            .property(key)
            .expect("complete hygroscopic requirements");
        let selected = property.answer().receipt.selected;
        let claim = card
            .claims()
            .claim(selected)
            .ok_or_else(|| invalid("selected hygroscopic source"))?;
        if claim
            .value
            .axes()
            .any(|(axis, dims)| axis == "moisture_content" && dims != Dims::NONE)
        {
            return Err(invalid("dry-basis source coordinate must be dimensionless"));
        }
        let path = card
            .claims()
            .query_line_typed(
                &claim.key,
                "moisture_content",
                lower,
                upper,
                fs_matdb::ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: key.into(),
                source: Box::new(source),
            })?;
        let knots = path
            .knots
            .into_iter()
            .map(|(x, a)| (x, a.evidence.value.value))
            .collect::<Vec<_>>();
        if index == 3 {
            if knots
                .iter()
                .any(|(_, value)| value.to_bits() != from.to_bits())
            {
                return Err(invalid(
                    "hygroscopic reference moisture differs from the source length datum",
                ));
            }
        } else {
            free_principal_strains[index] = integrate_piecewise_linear(&knots, from, to)?;
            if !free_principal_strains[index].is_finite() {
                return Err(invalid("integrated hygroscopic strain is unrepresentable"));
            }
        }
    }
    let mut identity = DomainHasher::new("org.frankensim.fs-material.hygroscopic-expansion.v1");
    identity.update(card.content_hash().as_bytes());
    identity.update(reference.identity().as_bytes());
    identity.update(current.identity().as_bytes());
    for value in free_principal_strains {
        identity.update(&value.to_bits().to_le_bytes());
    }
    Ok(IntegratedOrthotropicHygroscopicExpansion {
        reference,
        current,
        free_principal_strains,
        identity: identity.finalize(),
    })
}
