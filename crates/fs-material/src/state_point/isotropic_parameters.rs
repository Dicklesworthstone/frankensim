//! Explicit isotropic elastic parameter changes of variables.
//!
//! An independent pair is a caller-selected constitutive assumption. It is
//! never inferred from a material name, tensor component, or available fields.
//! The original scalar resolution remains intact; derived constants are not
//! inserted into the source card or represented as measured property claims.

use super::*;

/// Independent source parameters for one small-strain isotropic law.
///
/// Requirements may use arbitrary source names and retain their complete keys
/// and admissibility rules. This scalar adapter requires dimension-only
/// density, moduli and Poisson ratio: no existing semantic kind denotes an
/// elastic modulus. Stress, pressure, hardness and tensor-component keys cannot
/// silently stand in for elastic constants. Query axes retain their exact
/// typed temperature, moisture, pressure, process and other source conditions.
#[derive(Clone, Debug, PartialEq)]
pub enum IsotropicElasticParameters {
    /// Young's modulus [Pa] and Poisson ratio.
    YoungPoisson {
        /// Selected Young's-modulus requirement.
        young_modulus: ScalarPropertyRequirement,
        /// Selected Poisson-ratio requirement.
        poisson_ratio: ScalarPropertyRequirement,
    },
    /// Young's modulus [Pa] and shear modulus [Pa].
    YoungShear {
        /// Selected Young's-modulus requirement.
        young_modulus: ScalarPropertyRequirement,
        /// Selected shear-modulus requirement.
        shear_modulus: ScalarPropertyRequirement,
    },
    /// Bulk modulus [Pa] and shear modulus [Pa].
    BulkShear {
        /// Selected bulk-modulus requirement.
        bulk_modulus: ScalarPropertyRequirement,
        /// Selected shear-modulus requirement.
        shear_modulus: ScalarPropertyRequirement,
    },
}

impl IsotropicElasticParameters {
    fn members(&self) -> (&ScalarPropertyRequirement, &ScalarPropertyRequirement, u8) {
        match self {
            Self::YoungPoisson {
                young_modulus,
                poisson_ratio,
            } => (young_modulus, poisson_ratio, 0),
            Self::YoungShear {
                young_modulus,
                shear_modulus,
            } => (young_modulus, shear_modulus, 1),
            Self::BulkShear {
                bulk_modulus,
                shear_modulus,
            } => (bulk_modulus, shear_modulus, 2),
        }
    }
}

/// Algebraic derivation from the retained independent source parameters.
///
/// These are nominal isotropic-model consequences, not new measurements.
/// Source marginal uncertainty stays in the original usage receipts. No joint
/// covariance is inferred and no uncertainty interval is propagated by this
/// adapter. It neither validates the material's isotropy nor reconciles
/// additional, unselected measurements.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotropicElasticDerivation {
    density_requirement: ScalarPropertyRequirement,
    parameters: IsotropicElasticParameters,
    source_identity: ContentHash,
    shear_modulus_pa: f64,
    bulk_modulus_pa: f64,
    identity: ContentHash,
}

impl IsotropicElasticDerivation {
    /// Original selected density key and consumer admissibility.
    #[must_use]
    pub const fn density_requirement(&self) -> &ScalarPropertyRequirement {
        &self.density_requirement
    }

    /// Independent source roles and their exact requirements.
    #[must_use]
    pub const fn parameters(&self) -> &IsotropicElasticParameters {
        &self.parameters
    }

    /// Unmodified identity of the original evidence-bearing source bundle.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Nominal shear modulus [Pa] implied by the selected isotropic model.
    #[must_use]
    pub const fn shear_modulus_pa(&self) -> f64 {
        self.shear_modulus_pa
    }

    /// Nominal bulk modulus [Pa] implied by the selected isotropic model.
    #[must_use]
    pub const fn bulk_modulus_pa(&self) -> f64 {
        self.bulk_modulus_pa
    }

    /// Binds source receipts, their parameter roles and the derived constants.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

impl IsotropicElasticStatePoint {
    /// Explicit parameter derivation, or None for the original canonical
    /// density/Young/Poisson resolver.
    #[must_use]
    pub const fn derivation(&self) -> Option<&IsotropicElasticDerivation> {
        self.derivation.as_ref()
    }

    /// Identity of the physical elastic state, including its parameter
    /// derivation when present. Canonical states retain their original identity.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.derivation.as_ref().map_or_else(
            || self.resolved.identity(),
            IsotropicElasticDerivation::identity,
        )
    }
}

/// Resolve an explicitly selected independent isotropic parameter pair.
///
/// All three source properties are evaluated atomically at the complete query
/// point using the existing claim-selection and validity machinery. Caller
/// admissibility is preserved and the constitutive positivity/stability
/// conditions are additionally enforced. Missing constants are derived only
/// from this chosen pair; no alternative pair or source is silently selected.
///
/// The resulting existing elastic-state type feeds ordinary solid mechanics.
/// Its `resolved()` bundle contains exactly the original source properties;
/// `derivation()` describes the algebra and `identity()` binds both. Re-query
/// this function at each changed physical state; a frozen value does not imply
/// temperature-, rate-, moisture- or history-independent behavior.
pub fn resolve_isotropic_elastic_state_point_from_parameters(
    card: &MaterialCard,
    point: &QueryPoint,
    density: &ScalarPropertyRequirement,
    parameters: &IsotropicElasticParameters,
    selection: MaterialPropertySelection,
) -> Result<IsotropicElasticStatePoint, MaterialStatePointError> {
    let (first, second, tag) = parameters.members();
    require_scalar_parameter(density, Density::DIMS)?;
    require_scalar_parameter(first, Pressure::DIMS)?;
    require_scalar_parameter(second, if tag == 0 { Dims::NONE } else { Pressure::DIMS })?;
    let requirements = [density.clone(), first.clone(), second.clone()];
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let value = |requirement: &ScalarPropertyRequirement| {
        resolved
            .property(requirement.name())
            .expect("selected independent source requirement was resolved")
            .value_si()
    };
    let rho = value(density);
    let a = value(first);
    let b = value(second);
    if !(rho.is_finite() && rho > 0.0 && a.is_finite() && a > 0.0 && b.is_finite()) {
        return Err(invalid_parameters());
    }
    let (young, poisson, shear, bulk) = match tag {
        0 => {
            if !(b > -1.0 && b < 0.5) {
                return Err(invalid_parameters());
            }
            (a, b, (0.5 * a) / (1.0 + b), (a / 3.0) / (1.0 - 2.0 * b))
        }
        1 => {
            if b <= 0.0 {
                return Err(invalid_parameters());
            }
            let ratio = a / b;
            let poisson = 0.5 * ratio - 1.0;
            (a, poisson, b, (a / 3.0) / (3.0 - ratio))
        }
        _ => {
            if b <= 0.0 {
                return Err(invalid_parameters());
            }
            // Normalize by the larger modulus before forming 9KG/(3K+G).
            // No K*G product or 3K sum can overflow before the final result.
            let (young, poisson) = if a >= b {
                let ratio = b / a;
                (
                    b * (9.0 / (3.0 + ratio)),
                    (3.0 - 2.0 * ratio) / (2.0 * (3.0 + ratio)),
                )
            } else {
                let ratio = a / b;
                (
                    a * (9.0 / (3.0 * ratio + 1.0)),
                    (3.0 * ratio - 2.0) / (2.0 * (3.0 * ratio + 1.0)),
                )
            };
            (young, poisson, b, a)
        }
    };
    if [young, shear, bulk]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || !(poisson.is_finite() && poisson > -1.0 && poisson < 0.5)
    {
        return Err(invalid_parameters());
    }
    // The existing solid consumer reconstructs these Lame coefficients.
    // Refuse a numerically unrepresentable law before its infallible ingress.
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let twice_mu = 2.0 * (young / (2.0 * (1.0 + poisson)));
    let diagonal = lambda + twice_mu;
    if !lambda.is_finite()
        || !twice_mu.is_finite()
        || twice_mu <= 0.0
        || !diagonal.is_finite()
        || diagonal <= 0.0
    {
        return Err(invalid_parameters());
    }
    // Check the represented normal block with the solid consumer's Cholesky
    // arithmetic. Positive mathematical moduli can lose definiteness at the
    // floating-point limit; the remaining three shear pivots are twice_mu.
    let l10 = lambda / diagonal.sqrt();
    let pivot1 = l10.mul_add(-l10, diagonal);
    if !(pivot1.is_finite() && pivot1 > 0.0) {
        return Err(invalid_parameters());
    }
    let l21 = l10.mul_add(-l10, lambda) / pivot1.sqrt();
    let pivot2 = l21.mul_add(-l21, pivot1);
    if !(l21.is_finite() && pivot2.is_finite() && pivot2 > 0.0) {
        return Err(invalid_parameters());
    }
    let source_identity = resolved.identity();
    let mut identity = DomainHasher::new("org.frankensim.fs-material.isotropic-parameters.v1");
    identity.update(source_identity.as_bytes());
    identity.update(&[tag]);
    // Source bundle order is canonical; these names additionally bind each
    // physical role, so swapping parameter meanings cannot preserve identity.
    for requirement in &requirements {
        identity.update(&(requirement.name().len() as u64).to_le_bytes());
        identity.update(requirement.name().as_bytes());
    }
    for value in [rho, young, poisson, shear, bulk] {
        identity.update(&value.to_bits().to_le_bytes());
    }
    let derivation = IsotropicElasticDerivation {
        density_requirement: density.clone(),
        parameters: parameters.clone(),
        source_identity,
        shear_modulus_pa: shear,
        bulk_modulus_pa: bulk,
        identity: identity.finalize(),
    };
    Ok(IsotropicElasticStatePoint {
        resolved,
        density_kg_m3: rho,
        young_modulus_pa: young,
        poisson_ratio: poisson,
        derivation: Some(derivation),
    })
}

fn invalid_parameters() -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived {
        quantity: "positive finite isotropic elastic moduli with -1 < nu < 0.5",
    }
}

fn require_scalar_parameter(
    requirement: &ScalarPropertyRequirement,
    dims: Dims,
) -> Result<(), MaterialStatePointError> {
    if requirement.quantity() != QuantitySpec::dimensional(dims)
        || requirement.hardness_test().is_some()
        || requirement.elastic_component().is_some()
        || requirement.strain_component().is_some()
        || requirement.stress_component().is_some()
    {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: requirement.name().to_owned(),
            reason: "isotropic parameter needs a dimension-only scalar of the correct units, not a stress, pressure, hardness or tensor component",
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SmallStrainLaw, elastic::IsotropicElastic};
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        InterpolationPolicy, PropertyClaim, PropertyValue, Provenance, UncertaintyModel,
    };
    use fs_qty::semantic::MoistureBasis;

    fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
        ScalarPropertyRequirement::try_with_key(
            &PropertyKey::new(name, dims),
            ScalarAdmissibility::Finite,
        )
        .unwrap()
    }

    fn parameters(tag: u8) -> IsotropicElasticParameters {
        match tag {
            0 => IsotropicElasticParameters::YoungPoisson {
                young_modulus: requirement("E-source", Pressure::DIMS),
                poisson_ratio: requirement("nu-source", Dims::NONE),
            },
            1 => IsotropicElasticParameters::YoungShear {
                young_modulus: requirement("E-source", Pressure::DIMS),
                shear_modulus: requirement("G-source", Pressure::DIMS),
            },
            _ => IsotropicElasticParameters::BulkShear {
                bulk_modulus: requirement("K-source", Pressure::DIMS),
                shear_modulus: requirement("G-source", Pressure::DIMS),
            },
        }
    }

    fn temperature() -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ))
    }
    fn moisture() -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Moisture(MoistureBasis::DryMass),
            ValueForm::Static,
        ))
    }
    fn point() -> QueryPoint {
        QueryPoint::new()
            .with_quantity("temperature", temperature(), 293.15)
            .unwrap()
            .with_quantity("moisture", moisture(), 0.12)
            .unwrap()
    }
    fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
        PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", temperature(), 273.15, 323.15)
                .with_quantity("moisture", moisture(), 0.10, 0.15),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("synthetic independent parameter fixture: {name}"),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        }
    }
    fn card(values: [f64; 5]) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for ((name, dims), value) in [
            ("rho-source", Density::DIMS),
            ("E-source", Pressure::DIMS),
            ("nu-source", Dims::NONE),
            ("G-source", Pressure::DIMS),
            ("K-source", Pressure::DIMS),
        ]
        .into_iter()
        .zip(values)
        {
            claims.insert_claim(claim(name, dims, value)).unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "solid".into(),
                process: "numerical reference".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }
    fn resolve(
        card: &MaterialCard,
        tag: u8,
    ) -> Result<IsotropicElasticStatePoint, MaterialStatePointError> {
        resolve_isotropic_elastic_state_point_from_parameters(
            card,
            &point(),
            &requirement("rho-source", Density::DIMS),
            &parameters(tag),
            MaterialPropertySelection::SingleClaimOnly,
        )
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= expected.abs().max(1.0) * 2e-14,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_three_parameterizations_produce_the_same_physical_stress() {
        let material = card([2500.0, 60e9, 0.25, 24e9, 40e9]);
        let mut identities = std::collections::BTreeSet::new();
        for tag in 0..3 {
            let state = resolve(&material, tag).unwrap();
            close(state.density_kg_m3(), 2500.0);
            close(state.young_modulus_pa(), 60e9);
            close(state.poisson_ratio(), 0.25);
            let derived = state.derivation().unwrap();
            close(derived.shear_modulus_pa(), 24e9);
            close(derived.bulk_modulus_pa(), 40e9);
            assert_eq!(derived.source_identity(), state.resolved().identity());
            assert_eq!(derived.identity(), state.identity());
            assert!(identities.insert(state.identity()));
            assert_eq!(state.resolved().properties().len(), 3);
            if tag != 0 {
                assert!(state.resolved().property("nu-source").is_none());
            }
            let law = IsotropicElastic::new(state.young_modulus_pa(), state.poisson_ratio(), 1e-3)
                .unwrap();
            // Independent Lame values: lambda=24 GPa, mu=24 GPa.
            let stress = law.stress(&[1e-4, 0.0, 0.0, 2e-5, 0.0, 0.0], &());
            for (actual, expected) in stress
                .into_iter()
                .zip([7.2e6, 2.4e6, 2.4e6, 0.96e6, 0.0, 0.0])
            {
                close(actual, expected);
            }
            for property in state.resolved().properties() {
                material
                    .claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.axis_quantities,
                    *point().axis_quantities()
                );
            }
        }
    }

    #[test]
    fn g0_canonical_resolver_retains_its_original_identity() {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            (DENSITY_PROPERTY, Density::DIMS, 2500.0),
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 60e9),
            (POISSON_RATIO_PROPERTY, Dims::NONE, 0.25),
        ] {
            claims.insert_claim(claim(name, dims, value)).unwrap();
        }
        let material = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "solid".into(),
                process: "canonical numerical reference".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap();
        let state = resolve_isotropic_elastic_state_point(
            &material,
            &point(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(state.derivation().is_none());
        assert_eq!(state.identity(), state.resolved().identity());
    }

    #[test]
    fn g0_preserves_full_source_conditions_and_rejects_wrong_quantity_roles() {
        let material = card([2500.0, 60e9, 0.25, 24e9, 40e9]);
        for point in [
            QueryPoint::new()
                .with_quantity("temperature", temperature(), 293.15)
                .unwrap(),
            point().with("temperature", 350.0).unwrap(),
            point()
                .with_quantity(
                    "moisture",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Moisture(MoistureBasis::WetMass),
                        ValueForm::Static,
                    )),
                    0.12,
                )
                .unwrap(),
        ] {
            assert!(
                resolve_isotropic_elastic_state_point_from_parameters(
                    &material,
                    &point,
                    &requirement("rho-source", Density::DIMS),
                    &parameters(1),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .is_err()
            );
        }
        let wrong = IsotropicElasticParameters::YoungShear {
            young_modulus: ScalarPropertyRequirement::try_with_quantity(
                "E-source",
                QuantitySpec::semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Static)),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
            shear_modulus: requirement("G-source", Pressure::DIMS),
        };
        assert!(matches!(
            resolve_isotropic_elastic_state_point_from_parameters(
                &material,
                &point(),
                &requirement("rho-source", Density::DIMS),
                &wrong,
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::InvalidRequirement { .. })
        ));
    }

    #[test]
    fn g0_refuses_unstable_or_unrepresentable_independent_pairs() {
        for (values, tag) in [
            ([2500.0, 60e9, 0.5, 24e9, 40e9], 0),
            ([2500.0, 60e9, -1.0, 24e9, 40e9], 0),
            ([2500.0, 80e9, 0.25, 24e9, 40e9], 1),
            ([2500.0, 60e9, 0.25, 0.0, 40e9], 1),
            ([2500.0, 60e9, 0.25, 24e9, -40e9], 2),
            ([0.0, 60e9, 0.25, 24e9, 40e9], 2),
            ([2500.0, 60e9, 0.25, 1.0, 1e300], 2),
        ] {
            assert!(matches!(
                resolve(&card(values), tag),
                Err(MaterialStatePointError::InvalidDerived { .. })
            ));
        }
    }

    #[test]
    fn g3_explicit_pins_change_the_derived_law_without_fusing_conflicting_claims() {
        let original = card([2500.0, 60e9, 0.25, 24e9, 40e9]);
        let mut claims = original.claims().clone();
        let other = claims
            .insert_claim(claim("G-source", Pressure::DIMS, 30e9))
            .unwrap();
        let material = MaterialCard::assemble(original.id().clone(), claims, Vec::new()).unwrap();
        assert!(resolve(&material, 1).is_err());
        let pins = |shear| {
            MaterialPropertySelection::PinnedByProperty(vec![
                (
                    "rho-source".into(),
                    material.claims().claims_for("rho-source")[0].0,
                ),
                (
                    "E-source".into(),
                    material.claims().claims_for("E-source")[0].0,
                ),
                ("G-source".into(), shear),
            ])
        };
        let mut states = Vec::new();
        for pin in [material.claims().claims_for("G-source")[0].0, other] {
            states.push(
                resolve_isotropic_elastic_state_point_from_parameters(
                    &material,
                    &point(),
                    &requirement("rho-source", Density::DIMS),
                    &parameters(1),
                    pins(pin),
                )
                .unwrap(),
            );
        }
        close(states[0].poisson_ratio(), 0.25);
        close(states[1].poisson_ratio(), 0.0);
        assert_ne!(states[0].identity(), states[1].identity());
        assert_ne!(
            states[0].resolved().identity(),
            states[1].resolved().identity()
        );
    }

    #[test]
    fn g3_modulus_rescaling_preserves_poisson_and_stress_scaling() {
        for scale in [1e-12, 1.0, 1e12] {
            let material = card([2500.0, 60e9 * scale, 0.25, 24e9 * scale, 40e9 * scale]);
            for tag in 0..3 {
                let state = resolve(&material, tag).unwrap();
                close(state.poisson_ratio(), 0.25);
                close(state.young_modulus_pa() / scale, 60e9);
                close(state.derivation().unwrap().bulk_modulus_pa() / scale, 40e9);
            }
        }
    }
}
