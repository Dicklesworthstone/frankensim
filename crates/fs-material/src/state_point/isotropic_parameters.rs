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
    /// Original source requirements, in physical-role order. Callers can add
    /// these to one larger multiphysics query without inventing derived claims.
    #[must_use]
    pub fn requirements(&self) -> [&ScalarPropertyRequirement; 2] {
        let (first, second, _) = self.members();
        [first, second]
    }

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

/// Stiffness projected from one complete material bundle and an explicit
/// isotropic parameter pair. Density is deliberately not required: a thermal
/// stiffness update must not replace a closed body's reference mass inventory.
///
/// Extra thermal, optical, electrical or other properties remain in `resolved`
/// and refresh atomically with the moduli. Derived elastic constants are model
/// consequences, not inserted source claims or propagated uncertainty bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedIsotropicElasticModuli {
    resolved: ResolvedMaterialStatePoint,
    parameters: IsotropicElasticParameters,
    young_modulus_pa: f64,
    poisson_ratio: f64,
    shear_modulus_pa: f64,
    bulk_modulus_pa: f64,
    identity: ContentHash,
}

impl ResolvedIsotropicElasticModuli {
    /// Adopt the original bundle without selecting sources again or dropping
    /// other physics. Both the original and this projection's admissibility
    /// rules apply. The pair is never inferred from the available properties.
    ///
    /// # Errors
    /// Refuses missing/wrong scalar keys, repeated roles, inadmissible source
    /// values, and an unstable or numerically unrepresentable isotropic law.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        parameters: &IsotropicElasticParameters,
    ) -> Result<Self, MaterialStatePointError> {
        let [young, poisson, shear, bulk] = elastic_moduli(&resolved, parameters)?;
        let (_, _, tag) = parameters.members();
        let mut identity = DomainHasher::new("org.frankensim.fs-material.isotropic-moduli.v1");
        identity.update(resolved.identity().as_bytes());
        identity.update(&[tag]);
        for requirement in parameters.requirements() {
            identity.update(&(requirement.name().len() as u64).to_le_bytes());
            identity.update(requirement.name().as_bytes());
            requirement.admissibility().encode(&mut identity);
        }
        for value in [young, poisson, shear, bulk] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            resolved,
            parameters: parameters.clone(),
            young_modulus_pa: young,
            poisson_ratio: poisson,
            shear_modulus_pa: shear,
            bulk_modulus_pa: bulk,
            identity: identity.finalize(),
        })
    }

    /// Entire original material bundle, including all other physical properties.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Explicit original source roles and consumer domains.
    #[must_use]
    pub const fn parameters(&self) -> &IsotropicElasticParameters {
        &self.parameters
    }

    /// Young's modulus [Pa] under the selected isotropic model.
    #[must_use]
    pub const fn young_modulus_pa(&self) -> f64 {
        self.young_modulus_pa
    }

    /// Poisson ratio under the selected isotropic model.
    #[must_use]
    pub const fn poisson_ratio(&self) -> f64 {
        self.poisson_ratio
    }

    /// Shear modulus [Pa] under the selected isotropic model.
    #[must_use]
    pub const fn shear_modulus_pa(&self) -> f64 {
        self.shear_modulus_pa
    }

    /// Bulk modulus [Pa] under the selected isotropic model.
    #[must_use]
    pub const fn bulk_modulus_pa(&self) -> f64 {
        self.bulk_modulus_pa
    }

    /// Binds the full material state, parameter roles, domains and response.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Refresh every retained property through the shared material owner, then
    /// rebuild the same elastic projection. No source switch or partial update
    /// is published. This endpoint query does not evolve strain/history or mass.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(self.resolved.at_query_point(card, point)?, &self.parameters)
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
    /// Adopt canonical density/Young/Poisson data from an existing complete
    /// bundle. Other physics and original property receipts remain intact.
    /// Physical bounds are checked even if the bundle requested only finiteness.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        let density = ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        let parameters = IsotropicElasticParameters::YoungPoisson {
            young_modulus: ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            poisson_ratio: ScalarPropertyRequirement::try_new(
                POISSON_RATIO_PROPERTY,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: -1.0,
                    upper: 0.5,
                },
            )?,
        };
        let rho = scalar_parameter_value(&resolved, &density)?;
        let [young, poisson, _, _] = elastic_moduli(&resolved, &parameters)?;
        Ok(Self {
            resolved,
            density_kg_m3: rho,
            young_modulus_pa: young,
            poisson_ratio: poisson,
            derivation: None,
        })
    }

    /// Adopt a complete bundle with explicit density and elastic source roles.
    /// The projection reuses the same algebra as density-independent stiffness;
    /// it never re-queries or synthesizes source properties.
    pub fn try_from_resolved_with_parameters(
        resolved: ResolvedMaterialStatePoint,
        density: &ScalarPropertyRequirement,
        parameters: &IsotropicElasticParameters,
    ) -> Result<Self, MaterialStatePointError> {
        require_scalar_parameter(density, Density::DIMS)?;
        let (first, second, tag) = parameters.members();
        if density.name() == first.name() || density.name() == second.name() {
            return Err(MaterialStatePointError::DuplicateRequirement {
                property: density.name().to_owned(),
            });
        }
        let rho = scalar_parameter_value(&resolved, density)?;
        if rho <= 0.0 {
            return Err(invalid_parameters());
        }
        let [young, poisson, shear, bulk] = elastic_moduli(&resolved, parameters)?;
        let requirements = [density, first, second];
        let source_identity = resolved.identity();
        // Existing resolvers used these exact requirements in the source
        // bundle. Preserve their v1 identity; stronger projection-only domains
        // additionally participate in a versioned preimage.
        let stronger_domains = requirements.iter().any(|requirement| {
            resolved
                .property_by_key(&PropertyKey::new(requirement.name(), requirement.dims()))
                .is_some_and(|property| property.requirement() != *requirement)
        });
        let domain = if stronger_domains {
            "org.frankensim.fs-material.isotropic-parameters.v2"
        } else {
            "org.frankensim.fs-material.isotropic-parameters.v1"
        };
        let mut identity = DomainHasher::new(domain);
        identity.update(source_identity.as_bytes());
        identity.update(&[tag]);
        for requirement in requirements {
            identity.update(&(requirement.name().len() as u64).to_le_bytes());
            identity.update(requirement.name().as_bytes());
            if stronger_domains {
                requirement.admissibility().encode(&mut identity);
            }
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
        Ok(Self {
            resolved,
            density_kg_m3: rho,
            young_modulus_pa: young,
            poisson_ratio: poisson,
            derivation: Some(derivation),
        })
    }

    /// Refresh the entire original bundle at the new physical point and rebuild
    /// its original parameterization with pinned source claims. This does not
    /// perform conservative mass transfer or constitutive-history integration.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        let resolved = self.resolved.at_query_point(card, point)?;
        match &self.derivation {
            Some(derivation) => Self::try_from_resolved_with_parameters(
                resolved,
                &derivation.density_requirement,
                &derivation.parameters,
            ),
            None => Self::try_from_resolved(resolved),
        }
    }

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

impl IsotropicSolidStatePoint {
    /// Adopt canonical density/Young/Poisson/yield data while retaining the
    /// complete bundle, including any coupled physical properties and receipts.
    pub fn try_from_resolved(resolved: ResolvedMaterialStatePoint)
        -> Result<Self, MaterialStatePointError>
    {
        let elastic = IsotropicElasticStatePoint::try_from_resolved(resolved)?;
        let yield_requirement = ScalarPropertyRequirement::try_new(
            YIELD_STRESS_PROPERTY, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive)?;
        Self::from_elastic(elastic, &yield_requirement)
    }

    /// Project an explicitly chosen isotropic pair and positive yield source
    /// from one complete material state. Derived E/nu values are never inserted
    /// into the original source properties or represented as measured claims.
    pub fn try_from_resolved_with_parameters(
        resolved: ResolvedMaterialStatePoint,
        density: &ScalarPropertyRequirement,
        parameters: &IsotropicElasticParameters,
        yield_requirement: &ScalarPropertyRequirement,
    ) -> Result<Self, MaterialStatePointError> {
        if std::iter::once(density).chain(parameters.requirements())
            .any(|requirement| requirement.name() == yield_requirement.name()) {
            return Err(MaterialStatePointError::DuplicateRequirement {
                property: yield_requirement.name().to_owned(),
            });
        }
        let elastic = IsotropicElasticStatePoint::try_from_resolved_with_parameters(
            resolved, density, parameters)?;
        Self::from_elastic(elastic, yield_requirement)
    }

    fn from_elastic(elastic: IsotropicElasticStatePoint,
        yield_requirement: &ScalarPropertyRequirement) -> Result<Self, MaterialStatePointError>
    {
        require_scalar_parameter(yield_requirement, Pressure::DIMS)?;
        let yield_stress_pa = scalar_parameter_value(elastic.resolved(), yield_requirement)?;
        if !yield_stress_pa.is_finite() || yield_stress_pa <= 0.0 {
            return Err(MaterialStatePointError::InvalidDerived { quantity: "positive_yield_stress_pa" });
        }
        let identity = if elastic.derivation().is_none() {
            elastic.resolved().identity()
        } else {
            let mut identity = DomainHasher::new("org.frankensim.fs-material.isotropic-solid-parameters.v1");
            identity.update(elastic.identity().as_bytes());
            identity.update(&(yield_requirement.name().len() as u64).to_le_bytes());
            identity.update(yield_requirement.name().as_bytes());
            yield_requirement.admissibility().encode(&mut identity);
            identity.update(&yield_stress_pa.to_bits().to_le_bytes());
            identity.finalize()
        };
        Ok(Self { elastic, yield_requirement: yield_requirement.clone(), yield_stress_pa, identity })
    }

    /// Refresh every original elastic, yield and coupled property atomically,
    /// with the same selected claims and physical source roles. This endpoint
    /// operation does not evolve contact geometry or constitutive history.
    pub fn at_query_point(&self, card: &MaterialCard, point: &QueryPoint)
        -> Result<Self, MaterialStatePointError>
    {
        Self::from_elastic(self.elastic.at_query_point(card, point)?, &self.yield_requirement)
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
/// `derivation()` describes the algebra and `identity()` binds both. Use
/// `at_query_point` for state changes that must retain these source selections;
/// a frozen value does not imply temperature-, rate-, moisture- or
/// history-independent behavior.
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
    IsotropicElasticStatePoint::try_from_resolved_with_parameters(resolved, density, parameters)
}

fn scalar_parameter_value(
    resolved: &ResolvedMaterialStatePoint,
    requirement: &ScalarPropertyRequirement,
) -> Result<f64, MaterialStatePointError> {
    let property = resolved
        .property_by_key(&PropertyKey::new(requirement.name(), requirement.dims()))
        .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
            property: requirement.name().to_owned(),
            reason: "elastic source role requires its exact dimension-only scalar key",
        })?;
    let value = property.value_si();
    if !requirement.admissibility().admits(value) {
        return Err(MaterialStatePointError::OutsideConsumerDomain {
            property: requirement.name().to_owned(),
            value,
            admissibility: requirement.admissibility(),
        });
    }
    Ok(value)
}

fn elastic_moduli(
    resolved: &ResolvedMaterialStatePoint,
    parameters: &IsotropicElasticParameters,
) -> Result<[f64; 4], MaterialStatePointError> {
    let (first, second, tag) = parameters.members();
    require_scalar_parameter(first, Pressure::DIMS)?;
    require_scalar_parameter(second, if tag == 0 { Dims::NONE } else { Pressure::DIMS })?;
    if first.name() == second.name() {
        return Err(MaterialStatePointError::DuplicateRequirement {
            property: first.name().to_owned(),
        });
    }
    let a = scalar_parameter_value(resolved, first)?;
    let b = scalar_parameter_value(resolved, second)?;
    if !(a.is_finite() && a > 0.0 && b.is_finite()) {
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
    Ok([young, poisson, shear, bulk])
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
    fn g3_contact_projection_retains_yield_loss_and_explicit_parameter_roles() {
        let base=card([2500.0,60e9,0.25,24e9,40e9]);
        let mut claims=base.claims().clone();
        claims.insert_claim(claim("yield-source",Pressure::DIMS,100e6)).unwrap();
        claims.insert_claim(claim("loss-source",fs_qty::DynViscosity::DIMS,10.0)).unwrap();
        let material=MaterialCard::assemble(base.id().clone(),claims,Vec::new()).unwrap();
        let density=requirement("rho-source",Density::DIMS);
        let yield_requirement=requirement("yield-source",Pressure::DIMS);
        let loss=requirement("loss-source",fs_qty::DynViscosity::DIMS);
        for tag in 0..3 {
            let parameters=parameters(tag);
            let mut requirements=vec![density.clone(),loss.clone()];
            requirements.extend(parameters.requirements().into_iter().cloned());
            let original=resolve_material_state_point(&material,&point(),&requirements,
                MaterialPropertySelection::SingleClaimOnly).unwrap();
            let expanded=original.with_additional_requirements(&material,
                &[yield_requirement.clone()],MaterialPropertySelection::SingleClaimOnly).unwrap();
            assert!(original.properties().iter().all(|p|expanded.properties().contains(p)));
            let solid=IsotropicSolidStatePoint::try_from_resolved_with_parameters(
                expanded.clone(),&density,&parameters,&yield_requirement).unwrap();
            assert_eq!(solid.resolved(),&expanded);
            close(solid.young_modulus_pa(),60e9);
            close(solid.poisson_ratio(),0.25);
            close(solid.yield_stress_pa(),100e6);
            assert_eq!(solid.elastic_state().derivation().unwrap().parameters(),&parameters);
            if tag==2 {
                assert!(solid.resolved().property("E-source").is_none());
                assert!(solid.resolved().property("nu-source").is_none());
            }
            let hot=solid.at_query_point(&material,&point().with("temperature",313.15).unwrap()).unwrap();
            assert_eq!(hot.resolved().properties().len(),5);
            assert_eq!(hot.yield_requirement(),&yield_requirement);
            assert_eq!(hot.elastic_state().derivation().unwrap().parameters(),&parameters);
            assert!(solid.at_query_point(&material,&point().with("temperature",330.0).unwrap()).is_err());
        }
        let mut requirements=vec![density.clone(),yield_requirement.clone(),loss];
        for (name,dims) in [("E-source",Pressure::DIMS),("nu-source",Dims::NONE),
            ("G-source",Pressure::DIMS),("K-source",Pressure::DIMS)] {
            requirements.push(requirement(name,dims));
        }
        let all=resolve_material_state_point(&material,&point(),&requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        let a=IsotropicSolidStatePoint::try_from_resolved_with_parameters(
            all.clone(),&density,&parameters(1),&yield_requirement).unwrap();
        let b=IsotropicSolidStatePoint::try_from_resolved_with_parameters(
            all,&density,&parameters(2),&yield_requirement).unwrap();
        assert_ne!(a.identity(),b.identity(),"the same bundle cannot erase different selected roles");
        let aa=a.reduced_modulus_against(&a).unwrap();
        let ab=a.reduced_modulus_against(&b).unwrap();
        close(aa.value_pa,ab.value_pa);
        assert_ne!(aa.identity,ab.identity);
    }

    #[test]
    fn g3_appending_yield_retains_original_pinned_receipts_with_competing_elastic_claims() {
        let mut claims=ClaimSet::new();
        let mut requirements=Vec::new();
        let mut pins=Vec::new();
        for (name,dims,value) in [("rho-source",Density::DIMS,2500.0),
            ("E-source",Pressure::DIMS,60e9),("G-source",Pressure::DIMS,24e9)] {
            let id=claims.insert_claim(claim(name,dims,value)).unwrap();
            requirements.push(requirement(name,dims));
            pins.push((name.to_owned(),id));
        }
        claims.insert_claim(claim("E-source",Pressure::DIMS,48e9)).unwrap();
        claims.insert_claim(claim("yield-source",Pressure::DIMS,100e6)).unwrap();
        let original_card=card([2500.0,60e9,0.25,24e9,40e9]);
        let material=MaterialCard::assemble(original_card.id().clone(),claims,Vec::new()).unwrap();
        assert!(resolve_material_state_point(&material,&point(),&requirements,
            MaterialPropertySelection::SingleClaimOnly).is_err());
        let original=resolve_material_state_point(&material,&point(),&requirements,
            MaterialPropertySelection::PinnedByProperty(pins)).unwrap();
        let yield_requirement=requirement("yield-source",Pressure::DIMS);
        let expanded=original.with_additional_requirements(&material,
            &[yield_requirement.clone()],MaterialPropertySelection::SingleClaimOnly).unwrap();
        assert!(original.properties().iter().all(|p|expanded.properties().contains(p)));
        assert_eq!(expanded.axis_quantities(),original.axis_quantities());
        let solid=IsotropicSolidStatePoint::try_from_resolved_with_parameters(expanded,
            &requirements[0],&parameters(1),&yield_requirement).unwrap();
        close(solid.young_modulus_pa(),60e9);
        assert!(original.with_additional_requirements(&material,&[requirements[1].clone()],
            MaterialPropertySelection::SingleClaimOnly).is_err());
        assert!(original.with_additional_requirements(&original_card,&[yield_requirement],
            MaterialPropertySelection::SingleClaimOnly).is_err());
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

    #[test]
    fn g1_stiffness_projection_needs_no_density_and_preserves_other_properties() {
        // Density deliberately cannot describe a body. A stiffness-only
        // projection must neither query it nor invent a replacement mass.
        let material = card([0.0, 60e9, 0.25, 24e9, 40e9]);
        let requirements = [
            requirement("E-source", Pressure::DIMS),
            requirement("nu-source", Dims::NONE),
            requirement("G-source", Pressure::DIMS),
            requirement("K-source", Pressure::DIMS),
        ];
        let source = resolve_material_state_point(
            &material,
            &point(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut identities = std::collections::BTreeSet::new();
        for tag in 0..3 {
            let moduli =
                ResolvedIsotropicElasticModuli::try_from_resolved(source.clone(), &parameters(tag))
                    .unwrap();
            assert_eq!(moduli.resolved(), &source);
            assert_eq!(moduli.parameters(), &parameters(tag));
            assert!(identities.insert(moduli.identity()));
            close(moduli.young_modulus_pa(), 60e9);
            close(moduli.poisson_ratio(), 0.25);
            close(moduli.shear_modulus_pa(), 24e9);
            close(moduli.bulk_modulus_pa(), 40e9);
            let law =
                IsotropicElastic::new(moduli.young_modulus_pa(), moduli.poisson_ratio(), 1e-3)
                    .unwrap();
            // Hydrostatic and tensor-shear loads independently check K and G.
            let stress = law.stress(&[1e-4, 1e-4, 1e-4, 2e-5, 0.0, 0.0], &());
            for (actual, expected) in stress.into_iter().zip([12e6, 12e6, 12e6, 0.96e6, 0.0, 0.0]) {
                close(actual, expected);
            }
        }
        assert!(
            IsotropicElasticStatePoint::try_from_resolved_with_parameters(
                source,
                &requirement("rho-source", Density::DIMS),
                &parameters(2),
            )
            .is_err()
        );
    }

    fn refreshing_bundle() -> (MaterialCard, ResolvedMaterialStatePoint) {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, dims, axis, axis_dims, knots) in [
            (
                "rho-source",
                Density::DIMS,
                "temperature",
                fs_qty::Temperature::DIMS,
                vec![(273.15, 2500.0), (323.15, 2000.0)],
            ),
            (
                "E-source",
                Pressure::DIMS,
                "temperature",
                fs_qty::Temperature::DIMS,
                vec![(273.15, 60e9), (323.15, 30e9)],
            ),
            (
                "G-source",
                Pressure::DIMS,
                "moisture",
                Dims::NONE,
                vec![(0.10, 24e9), (0.15, 12e9)],
            ),
            (
                THERMAL_CONDUCTIVITY_PROPERTY,
                THERMAL_CONDUCTIVITY_DIMS,
                "moisture",
                Dims::NONE,
                vec![(0.10, 2.0), (0.15, 4.0)],
            ),
        ] {
            let mut row = claim(name, dims, 0.0);
            row.value = PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: axis_dims,
                knots,
                dims,
            };
            row.interpolation = InterpolationPolicy::LinearInside;
            if name == THERMAL_CONDUCTIVITY_PROPERTY {
                row.validity =
                    row.validity
                        .with_quantity("temperature", temperature(), 273.15, 313.15);
            }
            claims.insert_claim(row).unwrap();
            requirements.push(requirement(name, dims));
        }
        // A conflicting source becomes applicable later. Refresh must keep the
        // originally selected G curve, not re-run an ambiguity policy or switch.
        let mut competitor = claim("G-source", Pressure::DIMS, 18e9);
        competitor.validity =
            competitor
                .validity
                .with_quantity("temperature", temperature(), 283.15, 323.15);
        claims.insert_claim(competitor).unwrap();
        let material = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic multiphysics".into(),
                phase: "solid".into(),
                process: "projection refresh fixture".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap();
        let cold = point()
            .with("temperature", 273.15)
            .unwrap()
            .with("moisture", 0.10)
            .unwrap();
        let source = resolve_material_state_point(
            &material,
            &cold,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        (material, source)
    }

    #[test]
    fn g3_complete_elastic_refresh_keeps_sources_and_other_physics_atomic() {
        let (material, source) = refreshing_bundle();
        let state = IsotropicElasticStatePoint::try_from_resolved_with_parameters(
            source.clone(),
            &requirement("rho-source", Density::DIMS),
            &parameters(1),
        )
        .unwrap();
        let moduli =
            ResolvedIsotropicElasticModuli::try_from_resolved(source.clone(), &parameters(1))
                .unwrap();
        let hot = source
            .query()
            .with("temperature", 303.15)
            .unwrap()
            .with("moisture", 0.13)
            .unwrap();
        let updated = state.at_query_point(&material, &hot).unwrap();
        let stiffness = moduli.at_query_point(&material, &hot).unwrap();
        close(updated.density_kg_m3(), 2200.0);
        close(updated.young_modulus_pa(), 42e9);
        close(updated.poisson_ratio(), 0.25);
        close(stiffness.shear_modulus_pa(), 16.8e9);
        close(stiffness.bulk_modulus_pa(), 28e9);
        close(
            stiffness
                .resolved()
                .property(THERMAL_CONDUCTIVITY_PROPERTY)
                .unwrap()
                .value_si(),
            3.2,
        );
        assert_eq!(updated.resolved(), stiffness.resolved());
        assert_eq!(updated.resolved().query(), hot);
        assert_eq!(updated.derivation().unwrap().parameters(), &parameters(1));
        assert_ne!(updated.identity(), state.identity());
        assert_ne!(stiffness.identity(), moduli.identity());
        for (before, after) in source
            .properties()
            .iter()
            .zip(updated.resolved().properties())
        {
            assert_eq!(before.requirement(), after.requirement());
            assert_eq!(
                before.answer().receipt.selected,
                after.answer().receipt.selected
            );
            material
                .claims()
                .verify_receipt(&after.answer().receipt)
                .unwrap();
        }
        assert!(
            resolve_isotropic_elastic_state_point_from_parameters(
                &material,
                &hot,
                &requirement("rho-source", Density::DIMS),
                &parameters(1),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err()
        );
        // Only the auxiliary heat property leaves support here. Dropping it
        // when adopting the elastic view would incorrectly admit this update.
        let unsupported = hot.clone().with("temperature", 320.0).unwrap();
        assert!(matches!(state.at_query_point(&material, &unsupported),
            Err(MaterialStatePointError::Query { property, .. }) if property == THERMAL_CONDUCTIVITY_PROPERTY));
        assert!(moduli.at_query_point(&material, &unsupported).is_err());
        assert_eq!(state.resolved(), &source);
        assert_eq!(moduli.resolved(), &source);
        assert_eq!(updated, state.at_query_point(&material, &hot).unwrap());
        assert_eq!(stiffness, moduli.at_query_point(&material, &hot).unwrap());
    }

    #[test]
    fn g0_projection_rechecks_roles_domains_and_derived_stability() {
        let material = card([2500.0, 60e9, 0.25, 24e9, 40e9]);
        let requirements = [
            requirement("E-source", Pressure::DIMS),
            requirement("G-source", Pressure::DIMS),
        ];
        let source = resolve_material_state_point(
            &material,
            &point(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let original = parameters(1);
        let mut bounded = original.clone();
        if let IsotropicElasticParameters::YoungShear { shear_modulus, .. } = &mut bounded {
            *shear_modulus = ScalarPropertyRequirement::try_new(
                "G-source",
                Pressure::DIMS,
                ScalarAdmissibility::OpenInterval {
                    lower: 20e9,
                    upper: 30e9,
                },
            )
            .unwrap();
        }
        let a =
            ResolvedIsotropicElasticModuli::try_from_resolved(source.clone(), &original).unwrap();
        let b =
            ResolvedIsotropicElasticModuli::try_from_resolved(source.clone(), &bounded).unwrap();
        assert_ne!(a.identity(), b.identity());
        if let IsotropicElasticParameters::YoungShear { shear_modulus, .. } = &mut bounded {
            *shear_modulus = ScalarPropertyRequirement::try_new(
                "G-source",
                Pressure::DIMS,
                ScalarAdmissibility::OpenInterval {
                    lower: 25e9,
                    upper: 30e9,
                },
            )
            .unwrap();
        }
        assert!(matches!(
            ResolvedIsotropicElasticModuli::try_from_resolved(source.clone(), &bounded),
            Err(MaterialStatePointError::OutsideConsumerDomain { .. })
        ));
        let repeated = IsotropicElasticParameters::YoungShear {
            young_modulus: requirements[0].clone(),
            shear_modulus: requirements[0].clone(),
        };
        assert!(matches!(
            ResolvedIsotropicElasticModuli::try_from_resolved(source.clone(), &repeated),
            Err(MaterialStatePointError::DuplicateRequirement { .. })
        ));
        assert!(ResolvedIsotropicElasticModuli::try_from_resolved(source, &parameters(2)).is_err());
        let unstable = card([2500.0, 80e9, 0.25, 24e9, 40e9]);
        let source = resolve_material_state_point(
            &unstable,
            &point(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(matches!(
            ResolvedIsotropicElasticModuli::try_from_resolved(source, &original),
            Err(MaterialStatePointError::InvalidDerived { .. })
        ));
    }
}
