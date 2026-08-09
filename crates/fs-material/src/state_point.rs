//! Evidence-bearing material-property resolution at one physical state point.
//!
//! This module is the shared bridge from immutable [`fs_matdb::MaterialCard`]
//! data to executable physics.  A consumer declares exactly which scalar
//! properties it needs, their dimensions, and their admissible numerical
//! domains.  Resolution evaluates every property at the same caller-supplied
//! condition point, retains every usage receipt, and publishes one bundle
//! identity.  No material name implies a property and no missing datum is
//! replaced by a representative preset.

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimId, ClaimSet, InterfaceSystemCard, MatDbError, MaterialAnswer, MaterialCard,
    MaterialStateId, PropertyUsageReceiptError, QueryPoint, SelectionPolicy,
};
use fs_qty::{Density, Dims, Pressure};

use crate::elastic::OrthotropicElastic;

/// Identity domain for a complete resolved material state-point bundle.
pub const MATERIAL_STATE_POINT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-material.state-point.v1";
/// Identity domain for a complete resolved ordered-interface property bundle.
pub const INTERFACE_STATE_POINT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.interface-state-point.v1";
/// Maximum scalar requirements admitted by one bounded resolution.
pub const MAX_MATERIAL_STATE_PROPERTIES: usize = 64;
/// Maximum UTF-8 bytes in one required property name.
pub const MAX_MATERIAL_PROPERTY_NAME_BYTES: usize = 128;

/// Canonical property key for bulk density in kg/m3.
pub const DENSITY_PROPERTY: &str = "density";
/// Canonical property key for isotropic Young's modulus in pascals.
pub const YOUNG_MODULUS_PROPERTY: &str = "young_modulus";
/// Canonical property key for isotropic Poisson ratio.
pub const POISSON_RATIO_PROPERTY: &str = "poisson_ratio";
/// Canonical property key for uniaxial yield stress in pascals.
pub const YIELD_STRESS_PROPERTY: &str = "yield_stress";
/// Canonical property key for isotropic linear thermal expansion [1/K].
pub const LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY: &str =
    "linear_thermal_expansion_coefficient";
/// Dimensions of an inverse thermodynamic-temperature interval [1/K].
pub const INVERSE_TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, -1, 0, 0]);
/// Canonical orthotropic Young's-modulus keys along material axes 1, 2, 3.
pub const ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES: [&str; 3] =
    ["young_modulus_1", "young_modulus_2", "young_modulus_3"];
/// Canonical orthotropic major Poisson-ratio keys `(nu12, nu13, nu23)`.
pub const ORTHOTROPIC_POISSON_RATIO_PROPERTIES: [&str; 3] =
    ["poisson_ratio_12", "poisson_ratio_13", "poisson_ratio_23"];
/// Canonical orthotropic shear-modulus keys `(G12, G23, G31)`.
pub const ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES: [&str; 3] =
    ["shear_modulus_12", "shear_modulus_23", "shear_modulus_31"];
/// Fixed visible wavelengths used by the current evidence-bearing complex-IOR
/// bridge [vacuum nm]. The physics value is the sampled constitutive response;
/// this grid is a bounded transport convention, not a material preset.
pub const VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM: [f64; 9] = [
    380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0,
];
/// Canonical material-card keys for the real part of visible complex index.
pub const VISIBLE_COMPLEX_IOR_ETA_PROPERTIES: [&str; 9] = [
    "optical_eta_380nm",
    "optical_eta_430nm",
    "optical_eta_480nm",
    "optical_eta_530nm",
    "optical_eta_580nm",
    "optical_eta_630nm",
    "optical_eta_680nm",
    "optical_eta_730nm",
    "optical_eta_780nm",
];
/// Canonical material-card keys for the nonnegative extinction coefficient.
pub const VISIBLE_COMPLEX_IOR_K_PROPERTIES: [&str; 9] = [
    "optical_k_380nm",
    "optical_k_430nm",
    "optical_k_480nm",
    "optical_k_530nm",
    "optical_k_580nm",
    "optical_k_630nm",
    "optical_k_680nm",
    "optical_k_730nm",
    "optical_k_780nm",
];
/// Canonical material-card key for the dimensionless Cauchy `A` coefficient.
pub const VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY: &str = "optical_cauchy_a";
/// Canonical material-card key for the Cauchy `B` coefficient in SI m2.
pub const VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY: &str = "optical_cauchy_b_m2";
/// Canonical material-card key for the Cauchy `C` coefficient in SI m4.
pub const VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY: &str = "optical_cauchy_c_m4";
/// Canonical material-card keys for linear-RGB reference transmittance.
pub const VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES: [&str; 3] = [
    "optical_transmittance_linear_r",
    "optical_transmittance_linear_g",
    "optical_transmittance_linear_b",
];
/// Canonical material-card key for the Beer-Lambert reference distance [m].
pub const VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY: &str =
    "optical_transmittance_reference_distance_m";

const LENGTH_SQUARED_DIMS: Dims = Dims([2, 0, 0, 0, 0, 0]);
const LENGTH_FOURTH_DIMS: Dims = Dims([4, 0, 0, 0, 0, 0]);

/// Numerical domain a resolved scalar must satisfy before a solver may use it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ScalarAdmissibility {
    /// Any finite scalar is admitted.
    Finite,
    /// The scalar must be strictly greater than zero.
    StrictlyPositive,
    /// The scalar must be greater than or equal to zero.
    NonNegative,
    /// The scalar must lie strictly between two finite ordered endpoints.
    OpenInterval {
        /// Exclusive lower endpoint.
        lower: f64,
        /// Exclusive upper endpoint.
        upper: f64,
    },
}

impl ScalarAdmissibility {
    fn validate(self) -> bool {
        match self {
            Self::Finite | Self::StrictlyPositive | Self::NonNegative => true,
            Self::OpenInterval { lower, upper } => {
                lower.is_finite() && upper.is_finite() && lower < upper
            }
        }
    }

    fn admits(self, value: f64) -> bool {
        if !value.is_finite() {
            return false;
        }
        match self {
            Self::Finite => true,
            Self::StrictlyPositive => value > 0.0,
            Self::NonNegative => value >= 0.0,
            Self::OpenInterval { lower, upper } => value > lower && value < upper,
        }
    }

    fn encode(self, hasher: &mut DomainHasher) {
        match self {
            Self::Finite => hasher.update(&[0]),
            Self::StrictlyPositive => hasher.update(&[1]),
            Self::NonNegative => hasher.update(&[2]),
            Self::OpenInterval { lower, upper } => {
                hasher.update(&[3]);
                hasher.update(&lower.to_bits().to_le_bytes());
                hasher.update(&upper.to_bits().to_le_bytes());
            }
        }
    }
}

/// One solver-declared scalar property requirement.
#[derive(Clone, Debug, PartialEq)]
pub struct ScalarPropertyRequirement {
    name: String,
    dims: Dims,
    admissibility: ScalarAdmissibility,
}

impl ScalarPropertyRequirement {
    /// Construct one bounded property requirement.
    pub fn try_new(
        name: impl Into<String>,
        dims: Dims,
        admissibility: ScalarAdmissibility,
    ) -> Result<Self, MaterialStatePointError> {
        let name = name.into();
        if name.trim().is_empty() || name.len() > MAX_MATERIAL_PROPERTY_NAME_BYTES {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: name,
                reason: "property name must be nonblank and within its byte cap",
            });
        }
        if !admissibility.validate() {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: name,
                reason: "scalar admissibility domain is malformed",
            });
        }
        Ok(Self {
            name,
            dims,
            admissibility,
        })
    }

    /// Exact property key queried from the material card.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Dimensions required by the consuming physics operator.
    #[must_use]
    pub const fn dims(&self) -> Dims {
        self.dims
    }

    /// Numerical domain required by the consuming physics operator.
    #[must_use]
    pub const fn admissibility(&self) -> ScalarAdmissibility {
        self.admissibility
    }
}

/// Explicit claim-selection decision applied to a property query.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MaterialPropertySelection {
    /// Exactly one in-domain claim must exist.
    SingleClaimOnly,
    /// Prefer an observation-backed claim, but still refuse a surviving tie.
    PreferObservationBacked,
    /// Use one caller-pinned immutable claim identity for every named property.
    /// The plan must cover the requirement set exactly; omissions and foreign
    /// property names refuse before any result is published.
    PinnedByProperty(Vec<(String, ClaimId)>),
}

/// One property value and the complete matdb evidence/usage receipt behind it.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedScalarProperty {
    requirement: ScalarPropertyRequirement,
    answer: MaterialAnswer,
}

impl ResolvedScalarProperty {
    /// Requirement that this answer satisfied.
    #[must_use]
    pub const fn requirement(&self) -> &ScalarPropertyRequirement {
        &self.requirement
    }

    /// Resolved SI scalar.
    #[must_use]
    pub fn value_si(&self) -> f64 {
        self.answer.evidence.value.value
    }

    /// Complete evidence and usage receipt from fs-matdb.
    #[must_use]
    pub const fn answer(&self) -> &MaterialAnswer {
        &self.answer
    }
}

/// Canonically ordered property bundle for one named material and state point.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMaterialStatePoint {
    material: MaterialStateId,
    card_identity: ContentHash,
    query_point: Vec<(String, f64)>,
    properties: Vec<ResolvedScalarProperty>,
    identity: ContentHash,
}

/// Canonically ordered property bundle for one immutable ordered interface.
///
/// Interface properties are deliberately separate from bulk material
/// properties: friction, adhesion, conductance, and wear depend on both
/// surfaces, their texture frames, medium, environment, and history. The
/// complete [`InterfaceSystemCard`] identity binds those conditions.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedInterfaceStatePoint {
    card_identity: ContentHash,
    surface_materials: [MaterialStateId; 2],
    query_point: Vec<(String, f64)>,
    properties: Vec<ResolvedScalarProperty>,
    identity: ContentHash,
}

impl ResolvedInterfaceStatePoint {
    /// Immutable ordered-interface card identity.
    #[must_use]
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    /// Bulk material-state identities in ordered surface roles.
    #[must_use]
    pub const fn surface_materials(&self) -> &[MaterialStateId; 2] {
        &self.surface_materials
    }

    /// Canonically ordered physical condition coordinates.
    #[must_use]
    pub fn query_point(&self) -> &[(String, f64)] {
        &self.query_point
    }

    /// Look up one state coordinate by exact axis name.
    #[must_use]
    pub fn state_coordinate(&self, axis: &str) -> Option<f64> {
        self.query_point
            .binary_search_by(|(name, _)| name.as_str().cmp(axis))
            .ok()
            .map(|index| self.query_point[index].1)
    }

    /// Canonically ordered resolved interface properties.
    #[must_use]
    pub fn properties(&self) -> &[ResolvedScalarProperty] {
        &self.properties
    }

    /// Fetch one resolved property by its exact card key.
    #[must_use]
    pub fn property(&self, name: &str) -> Option<&ResolvedScalarProperty> {
        self.properties
            .binary_search_by(|property| property.requirement.name.as_str().cmp(name))
            .ok()
            .map(|index| &self.properties[index])
    }

    /// Identity binding the complete card, state point, requirements, values,
    /// and property-use receipts.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

impl ResolvedMaterialStatePoint {
    /// Named chemistry/phase/process/revision resolved by this bundle.
    #[must_use]
    pub const fn material(&self) -> &MaterialStateId {
        &self.material
    }

    /// Immutable material-card identity.
    #[must_use]
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    /// Canonically ordered physical condition coordinates.
    #[must_use]
    pub fn query_point(&self) -> &[(String, f64)] {
        &self.query_point
    }

    /// Canonically ordered resolved properties.
    #[must_use]
    pub fn properties(&self) -> &[ResolvedScalarProperty] {
        &self.properties
    }

    /// Identity binding the card, state point, requirements, values, and receipts.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Fetch one resolved property by its canonical key.
    #[must_use]
    pub fn property(&self, name: &str) -> Option<&ResolvedScalarProperty> {
        self.properties
            .binary_search_by(|property| property.requirement.name.as_str().cmp(name))
            .ok()
            .map(|index| &self.properties[index])
    }
}

/// Resolve a bounded set of scalar properties from one immutable material card.
///
/// Every property is evaluated at the same `point`. Requirements are
/// canonicalized by name, so caller ordering cannot move the bundle identity.
/// Missing, ambiguous, out-of-domain, dimensionally wrong, or numerically
/// inadmissible properties refuse the complete bundle; no partial state is
/// published.
pub fn resolve_material_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    requirements: &[ScalarPropertyRequirement],
    selection: MaterialPropertySelection,
) -> Result<ResolvedMaterialStatePoint, MaterialStatePointError> {
    let (query_point, properties) =
        resolve_scalar_property_set(card.claims(), point, requirements, &selection)?;
    let card_identity = card.content_hash();
    let identity = resolved_identity(
        MATERIAL_STATE_POINT_IDENTITY_DOMAIN,
        card_identity,
        &query_point,
        &properties,
    );
    Ok(ResolvedMaterialStatePoint {
        material: card.id().clone(),
        card_identity,
        query_point,
        properties,
        identity,
    })
}

/// Resolve a bounded set of scalar properties from one immutable ordered
/// interface card at one complete physical state point.
///
/// This applies the same atomic, dimension-checked, receipt-preserving query
/// contract as [`resolve_material_state_point`], while binding the result to
/// both surface roles, texture frames, medium, environment, and history via
/// the interface-card content identity.
pub fn resolve_interface_state_point(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    requirements: &[ScalarPropertyRequirement],
    selection: MaterialPropertySelection,
) -> Result<ResolvedInterfaceStatePoint, MaterialStatePointError> {
    let (query_point, properties) =
        resolve_scalar_property_set(card.claims(), point, requirements, &selection)?;
    let card_identity = card.content_hash();
    let identity = resolved_identity(
        INTERFACE_STATE_POINT_IDENTITY_DOMAIN,
        card_identity,
        &query_point,
        &properties,
    );
    Ok(ResolvedInterfaceStatePoint {
        card_identity,
        surface_materials: [
            card.surface_a().material.clone(),
            card.surface_b().material.clone(),
        ],
        query_point,
        properties,
        identity,
    })
}

fn resolve_scalar_property_set(
    claims: &ClaimSet,
    point: &QueryPoint,
    requirements: &[ScalarPropertyRequirement],
    selection: &MaterialPropertySelection,
) -> Result<(Vec<(String, f64)>, Vec<ResolvedScalarProperty>), MaterialStatePointError> {
    if requirements.is_empty() || requirements.len() > MAX_MATERIAL_STATE_PROPERTIES {
        return Err(MaterialStatePointError::RequirementCount {
            observed: requirements.len(),
            maximum: MAX_MATERIAL_STATE_PROPERTIES,
        });
    }
    let mut requirements = requirements.to_vec();
    requirements.sort_by(|left, right| left.name.cmp(&right.name));
    let mut names = BTreeSet::new();
    for requirement in &requirements {
        if !names.insert(requirement.name.clone()) {
            return Err(MaterialStatePointError::DuplicateRequirement {
                property: requirement.name.clone(),
            });
        }
    }

    let pins = match selection {
        MaterialPropertySelection::SingleClaimOnly
        | MaterialPropertySelection::PreferObservationBacked => None,
        MaterialPropertySelection::PinnedByProperty(offered) => {
            let mut pins = offered.clone();
            pins.sort_by(|left, right| left.0.cmp(&right.0));
            if pins.len() != requirements.len()
                || pins.windows(2).any(|pair| pair[0].0 == pair[1].0)
                || pins
                    .iter()
                    .zip(&requirements)
                    .any(|((name, _), requirement)| name != &requirement.name)
            {
                return Err(MaterialStatePointError::InvalidSelectionPlan);
            }
            Some(pins)
        }
    };

    let mut properties = Vec::with_capacity(requirements.len());
    for (index, requirement) in requirements.into_iter().enumerate() {
        let answer = match selection {
            MaterialPropertySelection::SingleClaimOnly => {
                claims.query(&requirement.name, point, SelectionPolicy::SingleClaimOnly)
            }
            MaterialPropertySelection::PreferObservationBacked => claims.query(
                &requirement.name,
                point,
                SelectionPolicy::PreferObservationBacked,
            ),
            MaterialPropertySelection::PinnedByProperty(_) => claims.query_pinned(
                &requirement.name,
                point,
                pins.as_ref().expect("pin plan was admitted")[index].1,
            ),
        }
        .map_err(|source| MaterialStatePointError::Query {
            property: requirement.name.clone(),
            source,
        })?;
        let found = answer.evidence.value.dims;
        if found != requirement.dims {
            return Err(MaterialStatePointError::DimensionMismatch {
                property: requirement.name,
                expected: requirement.dims,
                found,
            });
        }
        let value = answer.evidence.value.value;
        if !requirement.admissibility.admits(value) {
            return Err(MaterialStatePointError::OutsideConsumerDomain {
                property: requirement.name,
                value,
                admissibility: requirement.admissibility,
            });
        }
        answer
            .receipt
            .try_content_hash()
            .map_err(|source| MaterialStatePointError::Receipt {
                property: requirement.name.clone(),
                source,
            })?;
        properties.push(ResolvedScalarProperty {
            requirement,
            answer,
        });
    }

    let query_point = point
        .axes()
        .iter()
        .map(|(axis, value)| (axis.clone(), *value))
        .collect::<Vec<_>>();
    Ok((query_point, properties))
}

/// The three scalar properties needed by isotropic linear elasticity.
///
/// Yield is deliberately absent: free-vibration eigenmodes need density and
/// the tangent stiffness, while contact/plastic admission additionally needs
/// a yield surface. Requiring unrelated data would prevent valid acoustic
/// calculations; silently inventing it would be worse.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotropicElasticStatePoint {
    resolved: ResolvedMaterialStatePoint,
    density_kg_m3: f64,
    young_modulus_pa: f64,
    poisson_ratio: f64,
}

/// Evidence-bearing instantaneous isotropic thermal-expansion coefficient.
///
/// This is a state-point value, not a total strain. A thermomechanical driver
/// integrates it over the actual temperature path before constructing a solid
/// operator's stress-free strain state. Negative coefficients are admissible;
/// material names never select a sign or value.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotropicThermalExpansionStatePoint {
    resolved: ResolvedMaterialStatePoint,
    linear_coefficient_per_k: f64,
}

impl IsotropicThermalExpansionStatePoint {
    /// Complete evidence-bearing property bundle.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Instantaneous isotropic linear expansion coefficient [1/K].
    #[must_use]
    pub const fn linear_coefficient_per_k(&self) -> f64 {
        self.linear_coefficient_per_k
    }
}

/// Resolve the instantaneous isotropic thermal-expansion coefficient.
///
/// The material card owns temperature dependence and validity. This resolver
/// neither extrapolates nor replaces a missing coefficient with a preset.
pub fn resolve_isotropic_thermal_expansion_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<IsotropicThermalExpansionStatePoint, MaterialStatePointError> {
    let requirements = [ScalarPropertyRequirement::try_new(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        ScalarAdmissibility::Finite,
    )?];
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let linear_coefficient_per_k = resolved
        .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
        .expect("canonical isotropic thermal-expansion requirement was resolved")
        .value_si();
    Ok(IsotropicThermalExpansionStatePoint {
        resolved,
        linear_coefficient_per_k,
    })
}

/// Total isotropic free linear strain integrated over one temperature path.
#[derive(Clone, Debug, PartialEq)]
pub struct IntegratedIsotropicThermalExpansion {
    reference: IsotropicThermalExpansionStatePoint,
    current: IsotropicThermalExpansionStatePoint,
    selected_claim: ClaimId,
    free_linear_strain: f64,
    identity: ContentHash,
}

impl IntegratedIsotropicThermalExpansion {
    /// Expansion state at the path's reference temperature.
    #[must_use]
    pub const fn reference(&self) -> &IsotropicThermalExpansionStatePoint {
        &self.reference
    }

    /// Expansion state at the path's current temperature.
    #[must_use]
    pub const fn current(&self) -> &IsotropicThermalExpansionStatePoint {
        &self.current
    }

    /// Exact property claim integrated over the path.
    #[must_use]
    pub const fn selected_claim(&self) -> ClaimId {
        self.selected_claim
    }

    /// Signed free linear strain `integral(alpha(T) dT)`.
    #[must_use]
    pub const fn free_linear_strain(&self) -> f64 {
        self.free_linear_strain
    }

    /// Identity binding the card, path endpoints, selected curve, and result.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Integrate one selected isotropic expansion claim exactly over temperature.
///
/// Scalar claims use `alpha * delta_T`; piecewise-linear `alpha(T)` claims are
/// integrated segment by segment without quadrature error. Every non-temperature
/// coordinate must be identical at the two endpoints. The reference selection
/// chooses one claim and the current endpoint is then pinned to that exact claim,
/// preventing a path from silently switching evidence sources.
pub fn integrate_isotropic_thermal_expansion(
    card: &MaterialCard,
    reference_point: &QueryPoint,
    current_point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<IntegratedIsotropicThermalExpansion, MaterialStatePointError> {
    let (reference_temperature_k, current_temperature_k) =
        matching_temperature_path(reference_point, current_point)?;
    let reference =
        resolve_isotropic_thermal_expansion_state_point(card, reference_point, selection)?;
    let selected_claim = reference
        .resolved()
        .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
        .expect("canonical thermal-expansion requirement was resolved")
        .answer()
        .receipt
        .selected;
    let current = resolve_isotropic_thermal_expansion_state_point(
        card,
        current_point,
        MaterialPropertySelection::PinnedByProperty(vec![(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(),
            selected_claim,
        )]),
    )?;
    let claim =
        card.claims()
            .claim(selected_claim)
            .ok_or(MaterialStatePointError::InvalidDerived {
                quantity: "selected thermal-expansion claim",
            })?;
    let free_linear_strain = match &claim.value {
        fs_matdb::PropertyValue::Scalar { value, .. } => {
            value * (current_temperature_k - reference_temperature_k)
        }
        fs_matdb::PropertyValue::Curve {
            abscissa, knots, ..
        } if abscissa == "T" => {
            integrate_piecewise_linear(knots, reference_temperature_k, current_temperature_k)?
        }
        fs_matdb::PropertyValue::Curve { .. } => {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "thermal-expansion curve abscissa must be T",
            });
        }
    };
    if !free_linear_strain.is_finite() {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "integrated free linear strain",
        });
    }
    let mut identity =
        DomainHasher::new("org.frankensim.fs-material.integrated-isotropic-thermal-expansion.v1");
    identity.update(card.content_hash().as_bytes());
    identity.update(selected_claim.0.as_bytes());
    identity.update(reference.resolved().identity().as_bytes());
    identity.update(current.resolved().identity().as_bytes());
    identity.update(&free_linear_strain.to_bits().to_le_bytes());
    Ok(IntegratedIsotropicThermalExpansion {
        reference,
        current,
        selected_claim,
        free_linear_strain,
        identity: identity.finalize(),
    })
}

fn matching_temperature_path(
    reference: &QueryPoint,
    current: &QueryPoint,
) -> Result<(f64, f64), MaterialStatePointError> {
    if reference.axes().len() != current.axes().len() {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "thermal path coordinate set",
        });
    }
    for (axis, reference_value) in reference.axes() {
        let Some(current_value) = current.axes().get(axis) else {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "thermal path coordinate set",
            });
        };
        if axis != "T" && reference_value.to_bits() != current_value.to_bits() {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "non-temperature path coordinate changed",
            });
        }
    }
    let reference_temperature_k =
        reference
            .axes()
            .get("T")
            .copied()
            .ok_or(MaterialStatePointError::InvalidDerived {
                quantity: "reference temperature coordinate",
            })?;
    let current_temperature_k =
        current
            .axes()
            .get("T")
            .copied()
            .ok_or(MaterialStatePointError::InvalidDerived {
                quantity: "current temperature coordinate",
            })?;
    if reference_temperature_k <= 0.0 || current_temperature_k <= 0.0 {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "positive absolute temperature path",
        });
    }
    Ok((reference_temperature_k, current_temperature_k))
}

fn integrate_piecewise_linear(
    knots: &[(f64, f64)],
    start: f64,
    end: f64,
) -> Result<f64, MaterialStatePointError> {
    if start.to_bits() == end.to_bits() {
        return Ok(0.0);
    }
    let (lower, upper, sign) = if start < end {
        (start, end, 1.0)
    } else {
        (end, start, -1.0)
    };
    let mut cuts = Vec::with_capacity(knots.len() + 2);
    cuts.push(lower);
    cuts.extend(
        knots
            .iter()
            .map(|(temperature, _)| *temperature)
            .filter(|temperature| *temperature > lower && *temperature < upper),
    );
    cuts.push(upper);
    let interpolate = |temperature: f64| {
        let upper_index = knots.partition_point(|(knot, _)| *knot < temperature);
        if upper_index == 0 {
            return knots.first().map(|(_, value)| *value);
        }
        if upper_index == knots.len() {
            return knots.last().map(|(_, value)| *value);
        }
        let (x0, y0) = knots[upper_index - 1];
        let (x1, y1) = knots[upper_index];
        Some((temperature - x0).mul_add((y1 - y0) / (x1 - x0), y0))
    };
    let mut integral = 0.0;
    for pair in cuts.windows(2) {
        let y0 = interpolate(pair[0]).ok_or(MaterialStatePointError::InvalidDerived {
            quantity: "empty thermal-expansion curve",
        })?;
        let y1 = interpolate(pair[1]).ok_or(MaterialStatePointError::InvalidDerived {
            quantity: "empty thermal-expansion curve",
        })?;
        integral = (0.5 * (y0 + y1)).mul_add(pair[1] - pair[0], integral);
    }
    Ok(sign * integral)
}

impl IsotropicElasticStatePoint {
    /// Complete evidence-bearing property bundle.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Density at the queried state [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Young's modulus at the queried state [Pa].
    #[must_use]
    pub const fn young_modulus_pa(&self) -> f64 {
        self.young_modulus_pa
    }

    /// Isotropic Poisson ratio.
    #[must_use]
    pub const fn poisson_ratio(&self) -> f64 {
        self.poisson_ratio
    }
}

/// Resolve density and isotropic tangent elasticity at one exact state point.
pub fn resolve_isotropic_elastic_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<IsotropicElasticStatePoint, MaterialStatePointError> {
    let requirements = [
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        )?,
        ScalarPropertyRequirement::try_new(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let value = |name: &str| {
        resolved
            .property(name)
            .expect("canonical isotropic-elastic requirement was resolved")
            .value_si()
    };
    Ok(IsotropicElasticStatePoint {
        density_kg_m3: value(DENSITY_PROPERTY),
        young_modulus_pa: value(YOUNG_MODULUS_PROPERTY),
        poisson_ratio: value(POISSON_RATIO_PROPERTY),
        resolved,
    })
}

/// Evidence-bearing orthotropic tangent elasticity in its material frame.
///
/// Orientation is intentionally not part of this bulk property bundle. A
/// geometry/material-field consumer must supply and identity-bind the mapping
/// from material axes to its spatial frame.
#[derive(Clone, Debug, PartialEq)]
pub struct OrthotropicElasticStatePoint {
    resolved: ResolvedMaterialStatePoint,
    density_kg_m3: f64,
    law: OrthotropicElastic,
}

impl OrthotropicElasticStatePoint {
    /// Complete evidence-bearing property bundle.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Density at the queried state [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Admitted principal-axis orthotropic constitutive law.
    #[must_use]
    pub const fn law(&self) -> &OrthotropicElastic {
        &self.law
    }
}

/// Resolve a complete principal-axis orthotropic tangent at one state point.
///
/// All ten scalars resolve atomically from one immutable card and query point.
/// The derived compliance must be positive definite; no material name selects
/// anisotropy and no isotropic fallback exists.
pub fn resolve_orthotropic_elastic_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    strain_limit: f64,
) -> Result<OrthotropicElasticStatePoint, MaterialStatePointError> {
    if !(strain_limit.is_finite() && strain_limit > 0.0) {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "orthotropic_linear_strain_limit",
        });
    }
    let mut requirements = Vec::with_capacity(10);
    requirements.push(ScalarPropertyRequirement::try_new(
        DENSITY_PROPERTY,
        Density::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    for property in ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            property,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    for property in ORTHOTROPIC_POISSON_RATIO_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            property,
            Dims::NONE,
            ScalarAdmissibility::Finite,
        )?);
    }
    for property in ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            property,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let value = |name: &str| {
        resolved
            .property(name)
            .expect("canonical orthotropic requirement was resolved")
            .value_si()
    };
    let law = OrthotropicElastic::new(
        ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES.map(value),
        ORTHOTROPIC_POISSON_RATIO_PROPERTIES.map(value),
        ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES.map(value),
        strain_limit,
    )
    .map_err(|_| MaterialStatePointError::InvalidDerived {
        quantity: "orthotropic_compliance",
    })?;
    Ok(OrthotropicElasticStatePoint {
        density_kg_m3: value(DENSITY_PROPERTY),
        resolved,
        law,
    })
}

/// The four scalar properties needed by an isotropic elastic contact/body rung.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotropicSolidStatePoint {
    resolved: ResolvedMaterialStatePoint,
    density_kg_m3: f64,
    young_modulus_pa: f64,
    poisson_ratio: f64,
    yield_stress_pa: f64,
}

/// One evidence-bearing absolute complex-index sample at a material state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VisibleComplexIndexSample {
    /// Vacuum wavelength [nm].
    pub wavelength_nm: f64,
    /// Absolute real refractive index.
    pub eta: f64,
    /// Absolute extinction coefficient.
    pub k: f64,
}

/// Visible-band complex refractive index resolved from the same immutable
/// material-card/state-point machinery as mechanical properties.
#[derive(Clone, Debug, PartialEq)]
pub struct VisibleConductorStatePoint {
    resolved: ResolvedMaterialStatePoint,
    samples: [VisibleComplexIndexSample; 9],
}

/// Visible homogeneous-dielectric response resolved at one material state.
///
/// Cauchy coefficients are retained in SI (`B` in m2, `C` in m4). The
/// reference transmittance and distance define the homogeneous Beer-Lambert
/// absorption model consumed by the spectral renderer. Surface roughness is
/// deliberately not a bulk-material property and is bound separately.
#[derive(Clone, Debug, PartialEq)]
pub struct VisibleDielectricStatePoint {
    resolved: ResolvedMaterialStatePoint,
    cauchy_coefficients_si: [f64; 3],
    reference_transmittance_linear_rgb: [f64; 3],
    reference_distance_m: f64,
}

/// Data-selected visible optical constitutive family at one material state.
///
/// Selection depends only on which complete, admissible property schema the
/// immutable material card supplies. Material names and chemistry strings are
/// never inspected. A card that supplies both schemas is ambiguous unless a
/// pinned per-property selection plan makes exactly one family resolvable.
#[derive(Clone, Debug, PartialEq)]
pub enum VisibleOpticalStatePoint {
    /// Opaque complex-index response sampled over the visible band.
    Conductor(VisibleConductorStatePoint),
    /// Homogeneous transmitting response with Cauchy dispersion and
    /// Beer-Lambert absorption.
    Dielectric(VisibleDielectricStatePoint),
}

impl VisibleOpticalStatePoint {
    /// Complete card/state/property-use bundle for the selected family.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Conductor(state) => state.resolved(),
            Self::Dielectric(state) => state.resolved(),
        }
    }
}

impl VisibleDielectricStatePoint {
    /// Complete card/state/property-use bundle.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Cauchy `(A, B_m2, C_m4)` coefficients.
    #[must_use]
    pub const fn cauchy_coefficients_si(&self) -> [f64; 3] {
        self.cauchy_coefficients_si
    }

    /// Linear-RGB transmittance at [`Self::reference_distance_m`].
    #[must_use]
    pub const fn reference_transmittance_linear_rgb(&self) -> [f64; 3] {
        self.reference_transmittance_linear_rgb
    }

    /// Beer-Lambert reference distance [m].
    #[must_use]
    pub const fn reference_distance_m(&self) -> f64 {
        self.reference_distance_m
    }
}

impl VisibleConductorStatePoint {
    /// Complete card/state/property-use bundle.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Canonical visible-band complex-index samples.
    #[must_use]
    pub const fn samples(&self) -> &[VisibleComplexIndexSample; 9] {
        &self.samples
    }
}

/// Resolve a visible complex-index table at one exact material state point.
///
/// No chemistry name selects optical constants. All eighteen dimensionless
/// values must exist in the supplied card and be valid at `point`; the entire
/// table refuses atomically on missing, ambiguous, extrapolated, wrongly
/// dimensioned, or nonphysical data.
pub fn resolve_visible_conductor_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<VisibleConductorStatePoint, MaterialStatePointError> {
    let mut requirements = Vec::with_capacity(18);
    for name in VISIBLE_COMPLEX_IOR_ETA_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            Dims::NONE,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    for name in VISIBLE_COMPLEX_IOR_K_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )?);
    }
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let samples = core::array::from_fn(|index| VisibleComplexIndexSample {
        wavelength_nm: VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM[index],
        eta: resolved
            .property(VISIBLE_COMPLEX_IOR_ETA_PROPERTIES[index])
            .expect("canonical eta requirement was resolved")
            .value_si(),
        k: resolved
            .property(VISIBLE_COMPLEX_IOR_K_PROPERTIES[index])
            .expect("canonical extinction requirement was resolved")
            .value_si(),
    });
    Ok(VisibleConductorStatePoint { resolved, samples })
}

/// Resolve a homogeneous visible dielectric at one exact material state.
///
/// No chemistry name selects an optical model. All seven properties must be
/// present on one immutable material card and valid at the same query point.
/// Transmittance is admitted only in `(0, 1]`; malformed or partial optical
/// data refuses atomically rather than falling back to a visual preset.
pub fn resolve_visible_dielectric_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<VisibleDielectricStatePoint, MaterialStatePointError> {
    let mut requirements = Vec::with_capacity(7);
    requirements.push(ScalarPropertyRequirement::try_new(
        VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY,
        Dims::NONE,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    requirements.push(ScalarPropertyRequirement::try_new(
        VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
        LENGTH_SQUARED_DIMS,
        ScalarAdmissibility::NonNegative,
    )?);
    requirements.push(ScalarPropertyRequirement::try_new(
        VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
        LENGTH_FOURTH_DIMS,
        ScalarAdmissibility::NonNegative,
    )?);
    for name in VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            Dims::NONE,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    requirements.push(ScalarPropertyRequirement::try_new(
        VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
        Dims([1, 0, 0, 0, 0, 0]),
        ScalarAdmissibility::StrictlyPositive,
    )?);
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let value = |name: &str| {
        resolved
            .property(name)
            .expect("canonical visible-dielectric requirement was resolved")
            .value_si()
    };
    let reference_transmittance_linear_rgb = VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES.map(value);
    if reference_transmittance_linear_rgb
        .iter()
        .any(|value| *value > 1.0)
    {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "visible_dielectric_reference_transmittance",
        });
    }
    Ok(VisibleDielectricStatePoint {
        cauchy_coefficients_si: [
            value(VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY),
            value(VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY),
            value(VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY),
        ],
        reference_transmittance_linear_rgb,
        reference_distance_m: value(VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY),
        resolved,
    })
}

/// Resolve exactly one visible optical constitutive family from card data.
///
/// This is the generic material-to-rendering ingress. It admits a family only
/// when every required property resolves at the caller's exact state point.
/// Partial schemas never receive defaults, and two simultaneously complete
/// schemas refuse rather than being selected by material name or call order.
pub fn resolve_visible_optical_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<VisibleOpticalStatePoint, VisibleOpticalStatePointError> {
    let conductor = resolve_visible_conductor_state_point(card, point, selection.clone());
    let dielectric = resolve_visible_dielectric_state_point(card, point, selection);
    match (conductor, dielectric) {
        (Ok(conductor), Err(_)) => Ok(VisibleOpticalStatePoint::Conductor(conductor)),
        (Err(_), Ok(dielectric)) => Ok(VisibleOpticalStatePoint::Dielectric(dielectric)),
        (Ok(_), Ok(_)) => Err(VisibleOpticalStatePointError::AmbiguousFamilies),
        (Err(conductor), Err(dielectric)) => Err(VisibleOpticalStatePointError::NoCompleteFamily {
            conductor: Box::new(conductor),
            dielectric: Box::new(dielectric),
        }),
    }
}

impl IsotropicSolidStatePoint {
    /// Complete evidence-bearing property bundle.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Density at the queried state [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Young's modulus at the queried state [Pa].
    #[must_use]
    pub const fn young_modulus_pa(&self) -> f64 {
        self.young_modulus_pa
    }

    /// Isotropic Poisson ratio at the queried state.
    #[must_use]
    pub const fn poisson_ratio(&self) -> f64 {
        self.poisson_ratio
    }

    /// Yield stress at the queried state [Pa].
    #[must_use]
    pub const fn yield_stress_pa(&self) -> f64 {
        self.yield_stress_pa
    }

    /// Hertz reduced modulus against another isotropic solid [Pa].
    ///
    /// This is the standard two-half-space combination
    /// `1/E* = (1-nu_a^2)/E_a + (1-nu_b^2)/E_b`.  The returned identity binds
    /// both complete resolved material bundles in their ordered surface roles.
    pub fn reduced_modulus_against(
        &self,
        other: &Self,
    ) -> Result<ResolvedReducedModulus, MaterialStatePointError> {
        let compliance = (1.0 - self.poisson_ratio * self.poisson_ratio) / self.young_modulus_pa
            + (1.0 - other.poisson_ratio * other.poisson_ratio) / other.young_modulus_pa;
        let value_pa = compliance.recip();
        if !compliance.is_finite() || compliance <= 0.0 || !value_pa.is_finite() {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "hertz_reduced_modulus_pa",
            });
        }
        let mut hasher = DomainHasher::new("org.frankensim.fs-material.reduced-modulus.v1");
        hasher.update(self.resolved.identity().as_bytes());
        hasher.update(other.resolved.identity().as_bytes());
        hasher.update(&value_pa.to_bits().to_le_bytes());
        Ok(ResolvedReducedModulus {
            value_pa,
            surface_a_state_identity: self.resolved.identity(),
            surface_b_state_identity: other.resolved.identity(),
            identity: hasher.finalize(),
        })
    }
}

/// Resolve the canonical isotropic elastic properties at one state point.
pub fn resolve_isotropic_solid_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<IsotropicSolidStatePoint, MaterialStatePointError> {
    let requirements = [
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        )?,
        ScalarPropertyRequirement::try_new(
            YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let value = |name: &str| {
        resolved
            .property(name)
            .expect("canonical isotropic requirement was resolved")
            .value_si()
    };
    Ok(IsotropicSolidStatePoint {
        density_kg_m3: value(DENSITY_PROPERTY),
        young_modulus_pa: value(YOUNG_MODULUS_PROPERTY),
        poisson_ratio: value(POISSON_RATIO_PROPERTY),
        yield_stress_pa: value(YIELD_STRESS_PROPERTY),
        resolved,
    })
}

/// Ordered two-material Hertz elasticity derived from resolved state points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedReducedModulus {
    /// Reduced modulus [Pa].
    pub value_pa: f64,
    /// Ordered surface-A resolved-state identity.
    pub surface_a_state_identity: ContentHash,
    /// Ordered surface-B resolved-state identity.
    pub surface_b_state_identity: ContentHash,
    /// Identity binding both ordered states and the derived value.
    pub identity: ContentHash,
}

fn resolved_identity(
    domain: &'static str,
    card_identity: ContentHash,
    query_point: &[(String, f64)],
    properties: &[ResolvedScalarProperty],
) -> ContentHash {
    let mut hasher = DomainHasher::new(domain);
    hasher.update(card_identity.as_bytes());
    hasher.update(&(query_point.len() as u64).to_le_bytes());
    for (axis, value) in query_point {
        hasher.update(&(axis.len() as u64).to_le_bytes());
        hasher.update(axis.as_bytes());
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.update(&(properties.len() as u64).to_le_bytes());
    for property in properties {
        hasher.update(&(property.requirement.name.len() as u64).to_le_bytes());
        hasher.update(property.requirement.name.as_bytes());
        for exponent in property.requirement.dims.0 {
            hasher.update(&exponent.to_le_bytes());
        }
        property.requirement.admissibility.encode(&mut hasher);
        hasher.update(&property.value_si().to_bits().to_le_bytes());
        hasher.update(property.answer.receipt.content_hash().as_bytes());
    }
    hasher.finalize()
}

/// Typed, fail-closed material state-point refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialStatePointError {
    /// Requirement count was empty or exceeded the bounded bundle cap.
    RequirementCount {
        /// Offered count.
        observed: usize,
        /// Maximum admitted count.
        maximum: usize,
    },
    /// One requirement was malformed.
    InvalidRequirement {
        /// Offending property name.
        property: String,
        /// Stable explanation.
        reason: &'static str,
    },
    /// One property was required more than once.
    DuplicateRequirement {
        /// Duplicate property key.
        property: String,
    },
    /// A per-property pin plan did not cover the canonical requirement set exactly.
    InvalidSelectionPlan,
    /// The immutable material database refused a query.
    Query {
        /// Property being resolved.
        property: String,
        /// Exact matdb refusal, including out-of-domain and ambiguity cases.
        source: MatDbError,
    },
    /// The selected claim has different dimensions from the consuming operator.
    DimensionMismatch {
        /// Property being resolved.
        property: String,
        /// Consumer-required dimensions.
        expected: Dims,
        /// Claim dimensions.
        found: Dims,
    },
    /// A finite selected value lies outside the consumer's constitutive domain.
    OutsideConsumerDomain {
        /// Property being resolved.
        property: String,
        /// Selected SI value.
        value: f64,
        /// Required numerical domain.
        admissibility: ScalarAdmissibility,
    },
    /// A property receipt was not portable under the current receipt schema.
    Receipt {
        /// Property being resolved.
        property: String,
        /// Exact portable-receipt refusal.
        source: PropertyUsageReceiptError,
    },
    /// A derived physical property was non-finite or outside its mathematical domain.
    InvalidDerived {
        /// Stable derived quantity name.
        quantity: &'static str,
    },
}

impl fmt::Display for MaterialStatePointError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for MaterialStatePointError {}

/// Typed refusal from data-driven visible optical-family selection.
#[derive(Clone, Debug, PartialEq)]
pub enum VisibleOpticalStatePointError {
    /// Neither supported optical property schema resolved completely.
    NoCompleteFamily {
        /// Exact refusal from the complex-index conductor schema.
        conductor: Box<MaterialStatePointError>,
        /// Exact refusal from the Cauchy/Beer-Lambert dielectric schema.
        dielectric: Box<MaterialStatePointError>,
    },
    /// Both schemas resolved, so property presence alone cannot select one
    /// constitutive interpretation without inventing caller intent.
    AmbiguousFamilies,
}

impl fmt::Display for VisibleOpticalStatePointError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for VisibleOpticalStatePointError {}

#[cfg(test)]
mod tests {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
    };

    use super::*;

    fn claim(
        name: &str,
        dims: Dims,
        knots: Vec<(f64, f64)>,
        upper_temperature_k: f64,
    ) -> PropertyClaim {
        PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".to_owned(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots,
                dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, upper_temperature_k),
            uncertainty: UncertaintyModel::RelativeHalfWidth {
                fraction: 0.01,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("synthetic {name} curve"),
                license: "CC0-1.0".to_owned(),
                artifact: None,
            },
        }
    }

    fn solid_card(
        chemistry: &str,
        density: [f64; 2],
        young: [f64; 2],
        poisson: [f64; 2],
        yield_stress: [f64; 2],
    ) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for property in [
            claim(
                DENSITY_PROPERTY,
                Density::DIMS,
                vec![(250.0, density[0]), (600.0, density[1])],
                600.0,
            ),
            claim(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                vec![(250.0, young[0]), (600.0, young[1])],
                600.0,
            ),
            claim(
                POISSON_RATIO_PROPERTY,
                Dims::NONE,
                vec![(250.0, poisson[0]), (600.0, poisson[1])],
                600.0,
            ),
            claim(
                YIELD_STRESS_PROPERTY,
                Pressure::DIMS,
                vec![(250.0, yield_stress[0]), (600.0, yield_stress[1])],
                600.0,
            ),
            claim(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                vec![(250.0, 10.0e-6), (600.0, 24.0e-6)],
                600.0,
            ),
        ] {
            claims
                .insert_claim(property)
                .expect("valid synthetic property");
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: chemistry.to_owned(),
                phase: "solid".to_owned(),
                process: "synthetic-temperature-series".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .expect("material card")
    }

    fn point(temperature_k: f64) -> QueryPoint {
        QueryPoint::new()
            .with("T", temperature_k)
            .expect("finite temperature")
    }

    fn optical_card(upper_temperature_k: f64) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (index, name) in VISIBLE_COMPLEX_IOR_ETA_PROPERTIES.iter().enumerate() {
            claims
                .insert_claim(claim(
                    name,
                    Dims::NONE,
                    vec![
                        (250.0, 1.5 + index as f64 * 0.1),
                        (600.0, 1.4 + index as f64 * 0.1),
                    ],
                    upper_temperature_k,
                ))
                .unwrap();
        }
        for (index, name) in VISIBLE_COMPLEX_IOR_K_PROPERTIES.iter().enumerate() {
            claims
                .insert_claim(claim(
                    name,
                    Dims::NONE,
                    vec![
                        (250.0, 2.0 + index as f64 * 0.1),
                        (600.0, 2.5 + index as f64 * 0.1),
                    ],
                    upper_temperature_k,
                ))
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "test-visible-conductor".to_owned(),
                phase: "solid".to_owned(),
                process: "polished".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    fn dielectric_optical_card(red_transmittance: f64) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for property in [
            claim(
                VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY,
                Dims::NONE,
                vec![(250.0, 1.50), (600.0, 1.48)],
                600.0,
            ),
            claim(
                VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
                LENGTH_SQUARED_DIMS,
                vec![(250.0, 4.2e-15), (600.0, 4.0e-15)],
                600.0,
            ),
            claim(
                VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
                LENGTH_FOURTH_DIMS,
                vec![(250.0, 0.0), (600.0, 0.0)],
                600.0,
            ),
            claim(
                VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES[0],
                Dims::NONE,
                vec![(250.0, red_transmittance), (600.0, red_transmittance)],
                600.0,
            ),
            claim(
                VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES[1],
                Dims::NONE,
                vec![(250.0, 0.99), (600.0, 0.97)],
                600.0,
            ),
            claim(
                VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES[2],
                Dims::NONE,
                vec![(250.0, 0.96), (600.0, 0.92)],
                600.0,
            ),
            claim(
                VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
                Dims([1, 0, 0, 0, 0, 0]),
                vec![(250.0, 0.01), (600.0, 0.01)],
                600.0,
            ),
        ] {
            claims.insert_claim(property).unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "test-visible-dielectric".to_owned(),
                phase: "solid".to_owned(),
                process: "polished".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    fn ambiguous_optical_card() -> MaterialCard {
        let conductor = optical_card(600.0);
        let dielectric = dielectric_optical_card(0.98);
        let mut claims = conductor.claims().clone();
        for (_, claim) in dielectric.claims().claims_ordered() {
            claims.insert_claim(claim.clone()).unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "test-visible-ambiguous".to_owned(),
                phase: "solid".to_owned(),
                process: "polished".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    fn orthotropic_card(poisson: [f64; 3]) -> MaterialCard {
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(claim(
                DENSITY_PROPERTY,
                Density::DIMS,
                vec![(250.0, 800.0), (600.0, 760.0)],
                600.0,
            ))
            .unwrap();
        for (index, property) in ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES.iter().enumerate() {
            claims
                .insert_claim(claim(
                    property,
                    Pressure::DIMS,
                    vec![
                        (250.0, [12.0e9, 3.0e9, 1.5e9][index]),
                        (600.0, [10.0e9, 2.5e9, 1.2e9][index]),
                    ],
                    600.0,
                ))
                .unwrap();
        }
        for (index, property) in ORTHOTROPIC_POISSON_RATIO_PROPERTIES.iter().enumerate() {
            claims
                .insert_claim(claim(
                    property,
                    Dims::NONE,
                    vec![(250.0, poisson[index]), (600.0, poisson[index])],
                    600.0,
                ))
                .unwrap();
        }
        for (index, property) in ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES.iter().enumerate() {
            claims
                .insert_claim(claim(
                    property,
                    Pressure::DIMS,
                    vec![
                        (250.0, [1.2e9, 0.8e9, 0.6e9][index]),
                        (600.0, [1.0e9, 0.7e9, 0.5e9][index]),
                    ],
                    600.0,
                ))
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "test-orthotropic-solid".to_owned(),
                phase: "solid".to_owned(),
                process: "oriented-principal-axis-data".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    fn interface_card(history: &str) -> InterfaceSystemCard {
        let mut claims = ClaimSet::new();
        for property in [
            claim(
                "kinetic-friction-coefficient",
                Dims::NONE,
                vec![(250.0, 0.28), (600.0, 0.20)],
                600.0,
            ),
            claim(
                "adhesion-energy",
                Dims([0, 1, -2, 0, 0, 0]),
                vec![(250.0, 0.05), (600.0, 0.02)],
                600.0,
            ),
        ] {
            claims
                .insert_claim(property)
                .expect("valid synthetic interface property");
        }
        InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: MaterialStateId {
                    chemistry: "copper-c110".to_owned(),
                    phase: "solid".to_owned(),
                    process: "diamond-turned".to_owned(),
                    revision: 0,
                },
                texture_frame: "disc-edge/trace-17".to_owned(),
            },
            SurfaceSpec {
                material: MaterialStateId {
                    chemistry: "soda-lime-glass".to_owned(),
                    phase: "solid".to_owned(),
                    process: "float-polished".to_owned(),
                    revision: 0,
                },
                texture_frame: "base/trace-4".to_owned(),
            },
            SystemContext {
                medium: "dry".to_owned(),
                third_body: None,
                environment: "air".to_owned(),
                history: history.to_owned(),
            },
            claims,
            Vec::new(),
        )
        .expect("interface card")
    }

    #[test]
    fn g0_state_point_interpolates_every_property_and_is_order_invariant() {
        let lead = solid_card(
            "lead-pb99.99",
            [11_360.0, 11_100.0],
            [16.0e9, 8.0e9],
            [0.44, 0.46],
            [18.0e6, 3.0e6],
        );
        let first = resolve_isotropic_solid_state_point(
            &lead,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("solid lead state resolves inside all supplied domains");
        assert_eq!(first.density_kg_m3(), 11_230.0);
        assert_eq!(first.young_modulus_pa(), 12.0e9);
        assert_eq!(first.poisson_ratio(), 0.45);
        assert_eq!(first.resolved().properties().len(), 4);

        let mut requirements = [
            ScalarPropertyRequirement::try_new(
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ];
        let forward = resolve_material_state_point(
            &lead,
            &point(425.0),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        requirements.reverse();
        let reverse = resolve_material_state_point(
            &lead,
            &point(425.0),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(forward.identity(), reverse.identity());
    }

    #[test]
    fn g0_isotropic_thermal_expansion_is_card_resolved_and_domain_bounded() {
        let lead = solid_card(
            "lead-pb99.99",
            [11_360.0, 11_100.0],
            [16.0e9, 8.0e9],
            [0.44, 0.46],
            [18.0e6, 3.0e6],
        );
        let resolved = resolve_isotropic_thermal_expansion_state_point(
            &lead,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("temperature-dependent expansion resolves inside its evidence domain");
        assert_eq!(resolved.linear_coefficient_per_k(), 17.0e-6);
        assert_eq!(resolved.resolved().card_identity(), lead.content_hash());
        assert!(matches!(
            resolve_isotropic_thermal_expansion_state_point(
                &lead,
                &point(700.0),
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::Query {
                source: MatDbError::NoClaimInDomain { .. },
                ..
            })
        ));
    }

    #[test]
    fn g0_piecewise_linear_thermal_expansion_integrates_exactly_and_reverses() {
        let lead = solid_card(
            "lead-pb99.99",
            [11_360.0, 11_100.0],
            [16.0e9, 8.0e9],
            [0.44, 0.46],
            [18.0e6, 3.0e6],
        );
        let forward = integrate_isotropic_thermal_expansion(
            &lead,
            &point(300.0),
            &point(500.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("one selected alpha(T) curve spans the path");
        let reverse = integrate_isotropic_thermal_expansion(
            &lead,
            &point(500.0),
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("the same admitted path reverses");
        assert!((forward.free_linear_strain() - 0.0032).abs() < 1.0e-15);
        assert_eq!(
            reverse.free_linear_strain().to_bits(),
            (-forward.free_linear_strain()).to_bits()
        );
        assert_eq!(forward.selected_claim(), reverse.selected_claim());

        let reference = point(300.0).with("pressure", 1.0e5).unwrap();
        let current = point(500.0).with("pressure", 2.0e5).unwrap();
        assert!(matches!(
            integrate_isotropic_thermal_expansion(
                &lead,
                &reference,
                &current,
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::InvalidDerived {
                quantity: "non-temperature path coordinate changed"
            })
        ));
    }

    #[test]
    fn g0_visible_complex_index_is_state_resolved_and_extrapolation_refuses() {
        let card = optical_card(600.0);
        let resolved = resolve_visible_conductor_state_point(
            &card,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(resolved.samples().len(), 9);
        assert_eq!(resolved.samples()[0].wavelength_nm, 380.0);
        assert_eq!(resolved.samples()[8].wavelength_nm, 780.0);
        assert_eq!(resolved.samples()[0].eta, 1.45);
        assert_eq!(resolved.samples()[0].k, 2.25);
        assert_eq!(resolved.resolved().card_identity(), card.content_hash());

        assert!(matches!(
            resolve_visible_conductor_state_point(
                &card,
                &point(700.0),
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::Query {
                source: MatDbError::NoClaimInDomain { .. },
                ..
            })
        ));
    }

    #[test]
    fn g0_visible_dielectric_is_state_resolved_and_refuses_gain() {
        let card = dielectric_optical_card(0.98);
        let resolved = resolve_visible_dielectric_state_point(
            &card,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("complete homogeneous dielectric state resolves");
        let cauchy = resolved.cauchy_coefficients_si();
        assert!((cauchy[0] - 1.49).abs() <= 1.0e-15);
        assert!((cauchy[1] - 4.1e-15).abs() <= 1.0e-30);
        assert_eq!(cauchy[2], 0.0);
        assert_eq!(
            resolved.reference_transmittance_linear_rgb(),
            [0.98, 0.98, 0.94]
        );
        assert_eq!(resolved.reference_distance_m(), 0.01);
        assert_eq!(resolved.resolved().card_identity(), card.content_hash());

        assert!(matches!(
            resolve_visible_dielectric_state_point(
                &dielectric_optical_card(1.01),
                &point(425.0),
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::InvalidDerived {
                quantity: "visible_dielectric_reference_transmittance"
            })
        ));
    }

    #[test]
    fn g0_visible_optical_family_is_data_selected_and_ambiguous_cards_refuse() {
        let conductor = resolve_visible_optical_state_point(
            &optical_card(600.0),
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("complex-index schema selects conductor optics");
        assert!(matches!(conductor, VisibleOpticalStatePoint::Conductor(_)));

        let dielectric = resolve_visible_optical_state_point(
            &dielectric_optical_card(0.98),
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("Cauchy and transmittance schema selects dielectric optics");
        assert!(matches!(
            dielectric,
            VisibleOpticalStatePoint::Dielectric(_)
        ));

        assert!(matches!(
            resolve_visible_optical_state_point(
                &ambiguous_optical_card(),
                &point(425.0),
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(VisibleOpticalStatePointError::AmbiguousFamilies)
        ));
    }

    #[test]
    fn g0_elastic_resolvers_request_only_their_law_and_admit_orthotropy_atomically() {
        let complete_contact_card = solid_card(
            "test-isotropic-solid",
            [1_000.0, 900.0],
            [10.0e9, 8.0e9],
            [0.25, 0.30],
            [50.0e6, 30.0e6],
        );
        let elastic = resolve_isotropic_elastic_state_point(
            &complete_contact_card,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(elastic.resolved().properties().len(), 3);
        assert_eq!(elastic.density_kg_m3(), 950.0);
        assert_eq!(elastic.young_modulus_pa(), 9.0e9);
        assert_eq!(elastic.poisson_ratio(), 0.275);

        let card = orthotropic_card([0.25, 0.10, 0.20]);
        let orthotropic = resolve_orthotropic_elastic_state_point(
            &card,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly,
            1.0e-3,
        )
        .unwrap();
        assert_eq!(orthotropic.resolved().properties().len(), 10);
        assert_eq!(orthotropic.density_kg_m3(), 780.0);
        assert_eq!(orthotropic.law().e, [11.0e9, 2.75e9, 1.35e9]);
        assert_eq!(orthotropic.law().nu, [0.25, 0.10, 0.20]);
        assert_eq!(orthotropic.law().g, [1.1e9, 0.75e9, 0.55e9]);

        assert!(matches!(
            resolve_orthotropic_elastic_state_point(
                &orthotropic_card([2.0, 2.0, 2.0]),
                &point(425.0),
                MaterialPropertySelection::SingleClaimOnly,
                1.0e-3,
            ),
            Err(MaterialStatePointError::InvalidDerived {
                quantity: "orthotropic_compliance"
            })
        ));
    }

    #[test]
    fn g0_hot_solid_extrapolation_refuses_instead_of_reusing_room_temperature_data() {
        let lead = solid_card(
            "lead-pb99.99",
            [11_360.0, 11_100.0],
            [16.0e9, 8.0e9],
            [0.44, 0.46],
            [18.0e6, 3.0e6],
        );
        assert!(matches!(
            resolve_isotropic_solid_state_point(
                &lead,
                &point(700.0),
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::Query {
                source: MatDbError::NoClaimInDomain { .. },
                ..
            })
        ));
    }

    #[test]
    fn g0_reduced_modulus_binds_ordered_state_identities() {
        let lead = solid_card(
            "lead-pb99.99",
            [11_360.0, 11_100.0],
            [16.0e9, 8.0e9],
            [0.44, 0.46],
            [18.0e6, 3.0e6],
        );
        let glass = solid_card(
            "soda-lime-glass",
            [2500.0, 2480.0],
            [72.0e9, 68.0e9],
            [0.22, 0.23],
            [1.0e9, 0.8e9],
        );
        let lead = resolve_isotropic_solid_state_point(
            &lead,
            &point(293.15),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let glass = resolve_isotropic_solid_state_point(
            &glass,
            &point(293.15),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let lead_glass = lead.reduced_modulus_against(&glass).unwrap();
        let glass_lead = glass.reduced_modulus_against(&lead).unwrap();
        assert_eq!(lead_glass.value_pa.to_bits(), glass_lead.value_pa.to_bits());
        assert_ne!(lead_glass.identity, glass_lead.identity);
        assert_eq!(
            lead_glass.surface_a_state_identity,
            lead.resolved().identity()
        );
        assert_eq!(
            lead_glass.surface_b_state_identity,
            glass.resolved().identity()
        );
    }

    #[test]
    fn g0_interface_state_is_atomic_order_invariant_and_history_bound() {
        let requirements = [
            ScalarPropertyRequirement::try_new(
                "kinetic-friction-coefficient",
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "adhesion-energy",
                Dims([0, 1, -2, 0, 0, 0]),
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
        ];
        let virgin = interface_card("virgin-cleaned");
        let forward = resolve_interface_state_point(
            &virgin,
            &point(425.0),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("complete interface state resolves");
        let reverse = resolve_interface_state_point(
            &virgin,
            &point(425.0),
            &[requirements[1].clone(), requirements[0].clone()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("requirement order is not semantic");
        assert_eq!(forward.identity(), reverse.identity());
        let pins = requirements
            .iter()
            .rev()
            .map(|requirement| {
                (
                    requirement.name().to_owned(),
                    virgin.claims_for(requirement.name())[0].0,
                )
            })
            .collect();
        let pinned = resolve_interface_state_point(
            &virgin,
            &point(425.0),
            &requirements,
            MaterialPropertySelection::PinnedByProperty(pins),
        )
        .expect("exact per-property pins admit independently of offered order");
        assert_ne!(
            forward.identity(),
            pinned.identity(),
            "selection policy is provenance and remains identity-bearing"
        );
        let friction = forward
            .property("kinetic-friction-coefficient")
            .unwrap()
            .value_si();
        assert!((friction - 0.24).abs() < 1.0e-15);
        assert_ne!(
            forward.identity(),
            resolve_interface_state_point(
                &interface_card("run-in-1000-cycles"),
                &point(425.0),
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .expect("second complete state resolves")
            .identity(),
            "interface history is load-bearing even when values happen to match"
        );
    }

    #[test]
    fn g0_interface_state_refuses_partial_or_out_of_domain_resolution() {
        let requirements = [
            ScalarPropertyRequirement::try_new(
                "kinetic-friction-coefficient",
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "missing-wear-coefficient",
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
        ];
        let card = interface_card("virgin-cleaned");
        assert!(matches!(
            resolve_interface_state_point(
                &card,
                &point(425.0),
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::Query {
                source: MatDbError::UnknownProperty { .. },
                ..
            })
        ));
        assert!(matches!(
            resolve_interface_state_point(
                &card,
                &point(700.0),
                &requirements[..1],
                MaterialPropertySelection::SingleClaimOnly,
            ),
            Err(MaterialStatePointError::Query {
                source: MatDbError::NoClaimInDomain { .. },
                ..
            })
        ));
    }
}
