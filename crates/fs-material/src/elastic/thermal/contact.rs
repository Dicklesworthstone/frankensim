//! Source-derived two-temperature Hertz compliance and its thermal derivatives.
//!
//! Each original solid keeps its own temperature, complete fixed external
//! context, immutable card and selected claims. Each selected Young/Poisson,
//! Young/shear or bulk/shear pair is affine on its admitted source interval;
//! the corresponding isotropic moduli retain their exact rational dependence.
//! Yield limits remain constant. The ORIGINAL reduced-modulus owner evaluates
//! E* = [(1-nu_a^2)/E_a + (1-nu_b^2)/E_b]^-1. E* is generally NOT affine.
//!
//! Derivatives are at fixed indentation/geometry. If the native normal law has
//! psi=G(delta)*E*, its surface entropies are -G*dE*/dTi and its internal
//! contribution is psi + sum(Ti*si). The Hessian supplies the reciprocal heat
//! capacity. This module supplies no contact integrator, geometry update,
//! measured calibration, plasticity, phase change or expansion law.

use crate::field_storage::{FieldStorageError, thermal::affine_source_path};
use crate::state_point::{
    DENSITY_PROPERTY, IsotropicElasticParameters, IsotropicSolidStatePoint,
    MaterialStatePointError, POISSON_RATIO_PROPERTY, ResolvedReducedModulus, ScalarAdmissibility,
    YOUNG_MODULUS_PROPERTY, affine_temperature_slope,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

#[derive(Debug)]
pub enum ThermalHertzError {
    Material(MaterialStatePointError),
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for ThermalHertzError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for ThermalHertzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-derived thermal Hertz compliance: {self:?}")
    }
}
impl std::error::Error for ThermalHertzError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Invalid(_) => None,
        }
    }
}
fn finite(value: f64) -> Result<f64, ThermalHertzError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThermalHertzError::Invalid(
            "unrepresentable thermal contact derivative",
        ))
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Surface {
    card: MaterialCard,
    reference: IsotropicSolidStatePoint,
    temperature_range_k: [f64; 2],
    reference_temperature_k: f64,
    parameters: IsotropicElasticParameters,
    source_slopes_per_k: [f64; 2],
    source_path_scales: [f64; 2],
}
impl Surface {
    fn new(
        card: &MaterialCard,
        reference: &IsotropicSolidStatePoint,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, ThermalHertzError> {
        let invalid = ThermalHertzError::Invalid;
        let [low, high] = temperature_range_k;
        if !low.is_finite() || !high.is_finite() || low <= 0.0 || low >= high {
            return Err(invalid(
                "positive increasing source temperature interval required",
            ));
        }
        let original = reference.resolved();
        // This original operation pins the card, source claims and query schema.
        original.at_query_point(card, &original.query())?;
        let reference_temperature_k = *original
            .query()
            .axes()
            .get("T")
            .ok_or(invalid("canonical T coordinate is required"))?;
        if !(low..=high).contains(&reference_temperature_k) {
            return Err(invalid(
                "initial temperature leaves the admitted source interval",
            ));
        }
        let temperature_quantity = *original.axis_quantities().get("T").ok_or(invalid(
            "temperature needs explicit kelvin dimensions or absolute semantics",
        ))?;
        if temperature_quantity != QuantitySpec::dimensional(Temperature::DIMS)
            && temperature_quantity
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                ))
        {
            return Err(invalid("contact temperature must be absolute kelvin"));
        }
        let parameters = match reference.elastic_state().derivation() {
            Some(derivation) => derivation.parameters().clone(),
            None => IsotropicElasticParameters::YoungPoisson {
                young_modulus: original
                    .property(YOUNG_MODULUS_PROPERTY)
                    .ok_or(invalid("original Young source is missing or ambiguous"))?
                    .requirement()
                    .clone(),
                poisson_ratio: original
                    .property(POISSON_RATIO_PROPERTY)
                    .ok_or(invalid("original Poisson source is missing or ambiguous"))?
                    .requirement()
                    .clone(),
            },
        };
        let density_name = reference
            .elastic_state()
            .derivation()
            .map_or(DENSITY_PROPERTY, |d| d.density_requirement().name());
        let yield_name = reference.yield_requirement().name();
        let requirements = parameters.requirements();
        let mut source_slopes_per_k = [0.0; 2];
        let mut source_path_scales = [0.0; 2];
        // Every retained property must cover the admitted interval. Auxiliary
        // properties keep their own consumer domains and exact source receipts;
        // only the selected elastic pair supplies derivatives of this potential.
        for property in original.properties() {
            let name = property.requirement().name();
            let role = requirements.iter().position(|r| r.name() == name);
            let active = role.is_some() || name == density_name || name == yield_name;
            if active
                && !property.has_only_owned_and_fixed_external_conditions("T", temperature_quantity)
            {
                return Err(invalid(
                    "thermal contact requires one owned temperature and fixed external conditions",
                ));
            }
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(invalid("original solid source claim is absent"))?;
            if claim
                .value
                .axes()
                .any(|(axis, dims)| axis == "T" && dims != Temperature::DIMS)
            {
                return Err(invalid(
                    "source temperature axis has incompatible dimensions",
                ));
            }
            let point = |t| {
                original
                    .query()
                    .with("T", t)
                    .map_err(|_| invalid("thermal source query is unrepresentable"))
            };
            let line = card
                .claims()
                .query_line_typed(
                    &claim.key,
                    "T",
                    &point(low)?,
                    &point(high)?,
                    ClaimSelection::Pinned(property.answer().receipt.selected),
                )
                .map_err(|source| MaterialStatePointError::Support {
                    property: name.to_owned(),
                    source: Box::new(source),
                })?;
            let knots: Vec<_> = line
                .knots
                .into_iter()
                .map(|(t, answer)| (t, answer.evidence.value.value))
                .collect();
            let poisson = role == Some(1)
                && matches!(&parameters, IsotropicElasticParameters::YoungPoisson { .. });
            if knots.first().map(|k| k.0) != Some(low)
                || knots.last().map(|k| k.0) != Some(high)
                || knots.iter().any(|(_, value)| {
                    let consumer = match property.requirement().admissibility() {
                        ScalarAdmissibility::Finite => true,
                        ScalarAdmissibility::StrictlyPositive => *value > 0.0,
                        ScalarAdmissibility::NonNegative => *value >= 0.0,
                        ScalarAdmissibility::OpenInterval { lower, upper } => {
                            *value > lower && *value < upper
                        }
                    };
                    !value.is_finite()
                        || !consumer
                        || active
                            && if poisson {
                                *value <= -1.0 || *value >= 0.5
                            } else {
                                *value <= 0.0
                            }
                })
            {
                return Err(invalid(
                    "solid sources must continuously cover their complete admitted thermal interval",
                ));
            }
            if let Some(index) = role {
                let slope = if temperature_quantity
                    == QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )) {
                    // Share the original semantic source-path admission. Legacy
                    // dimension-only T keeps its exact schema and pinned-line
                    // admission above instead of being relabelled as semantic T.
                    affine_source_path(card, original, name, "T", temperature_range_k, !poisson)
                        .map_err(|error| match error {
                            FieldStorageError::Material(source) => {
                                ThermalHertzError::Material(source)
                            }
                            FieldStorageError::Invalid(message) => invalid(message),
                        })?
                        .1
                } else {
                    let (anchor_t, anchor_value) = if poisson {
                        knots[0]
                    } else {
                        (reference_temperature_k, property.value_si())
                    };
                    let slope = affine_temperature_slope(&knots, anchor_t, anchor_value)
                        .map_err(invalid)?;
                    let expected =
                        finite(slope.mul_add(reference_temperature_k - anchor_t, anchor_value))?;
                    let scale = knots.iter().fold(
                        property.value_si().abs().max(expected.abs()),
                        |scale, (_, value)| scale.max(value.abs()),
                    );
                    if (expected - property.value_si()).abs() > 128.0 * f64::EPSILON * scale {
                        return Err(invalid(
                            "signed thermal source reference left its admitted affine path",
                        ));
                    }
                    slope
                };
                if slope == 0.0 && knots.iter().any(|(_, v)| *v != property.value_si()) {
                    return Err(invalid(
                        "nonconstant thermal elastic source derivative underflows",
                    ));
                }
                source_slopes_per_k[index] = slope;
                source_path_scales[index] = knots
                    .iter()
                    .fold(0.0_f64, |scale, (_, value)| scale.max(value.abs()));
            } else if name == yield_name
                && knots
                    .iter()
                    .any(|(_, v)| v.to_bits() != property.value_si().to_bits())
            {
                return Err(invalid("thermal Hertz needs a constant source yield limit"));
            }
        }
        let value = Self {
            card: card.clone(),
            reference: reference.clone(),
            temperature_range_k,
            reference_temperature_k,
            parameters,
            source_slopes_per_k,
            source_path_scales,
        };
        for t in [low, reference_temperature_k, high] {
            value.at_temperature(t)?;
        }
        Ok(value)
    }

    fn at_temperature(&self, t: f64) -> Result<IsotropicSolidStatePoint, ThermalHertzError> {
        if !t.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&t)
        {
            return Err(ThermalHertzError::Invalid(
                "temperature leaves the pinned solid source interval",
            ));
        }
        let point = self
            .reference
            .resolved()
            .query()
            .with("T", t)
            .map_err(|_| ThermalHertzError::Invalid("thermal source query is unrepresentable"))?;
        let result = self.reference.at_query_point(&self.card, &point)?;
        for (index, requirement) in self.parameters.requirements().into_iter().enumerate() {
            let name = requirement.name();
            let reference = self
                .reference
                .resolved()
                .property(name)
                .ok_or(ThermalHertzError::Invalid(
                    "selected elastic source is missing or ambiguous",
                ))?
                .value_si();
            let actual = result
                .resolved()
                .property(name)
                .ok_or(ThermalHertzError::Invalid(
                    "refreshed elastic source is missing or ambiguous",
                ))?
                .value_si();
            let expected = finite(
                self.source_slopes_per_k[index]
                    .mul_add(t - self.reference_temperature_k, reference),
            )?;
            let mut scale = reference.abs().max(expected.abs()).max(actual.abs());
            if index == 1
                && matches!(
                    &self.parameters,
                    IsotropicElasticParameters::YoungPoisson { .. }
                )
            {
                // Equivalent interpolation orders lose relative precision at a
                // signed zero crossing. Admit against the original source path;
                // positive moduli retain their original relative check.
                scale = scale.max(self.source_path_scales[index]);
            }
            if (expected - actual).abs() > 512.0 * f64::EPSILON * scale {
                return Err(ThermalHertzError::Invalid(
                    "refreshed solid differs from its admitted affine elastic sources",
                ));
            }
        }
        if result.yield_stress_pa().to_bits() != self.reference.yield_stress_pa().to_bits() {
            return Err(ThermalHertzError::Invalid(
                "refreshed yield differs from its fixed source",
            ));
        }
        Ok(result)
    }

    fn fixed_poisson(&self) -> bool {
        matches!(
            &self.parameters,
            IsotropicElasticParameters::YoungPoisson { .. }
        ) && self.source_slopes_per_k[1] == 0.0
    }

    fn legacy_identity(&self) -> bool {
        self.reference.elastic_state().derivation().is_none() && self.fixed_poisson()
    }

    // Return c'/c and c''/c for c=(1-nu^2)/E. Ratios avoid forming squared
    // moduli or tiny compliance derivatives before multiplying by E*. All
    // parameter derivatives are those of the original selected affine sources.
    fn compliance_ratios(
        &self,
        state: &IsotropicSolidStatePoint,
    ) -> Result<[f64; 2], ThermalHertzError> {
        let requirements = self.parameters.requirements();
        let mut values = [0.0; 2];
        for i in 0..2 {
            values[i] = state
                .resolved()
                .property(requirements[i].name())
                .ok_or(ThermalHertzError::Invalid(
                    "selected thermal modulus is missing or ambiguous",
                ))?
                .value_si();
        }
        let [a, b] = values;
        let [at, bt] = self.source_slopes_per_k;
        let relative = |slope: f64, value: f64| -> Result<f64, ThermalHertzError> {
            let ratio = finite(slope / value)?;
            if slope != 0.0 && ratio == 0.0 {
                return Err(ThermalHertzError::Invalid(
                    "nonzero source derivative underflows",
                ));
            }
            Ok(ratio)
        };
        let result = match &self.parameters {
            IsotropicElasticParameters::YoungPoisson { .. } => {
                let re = relative(at, a)?;
                let d = finite(1.0 - b * b)?;
                let n = finite(-2.0 * b * bt / d)?;
                [
                    finite(n - re)?,
                    finite(
                        finite(-2.0 * bt * bt / d)? - finite(2.0 * n * re)?
                            + finite(2.0 * re * re)?,
                    )?,
                ]
            }
            IsotropicElasticParameters::YoungShear { .. } => {
                let re = relative(at, a)?;
                let rg = relative(bt, b)?;
                let ratio = finite(a / b)?;
                let denominator = finite(1.0 - 0.25 * ratio)?;
                [
                    finite((-rg + 0.25 * ratio * (2.0 * rg - re)) / denominator)?,
                    finite(
                        (2.0 * rg * rg - 0.5 * ratio * (3.0 * rg * rg - 2.0 * re * rg))
                            / denominator,
                    )?,
                ]
            }
            IsotropicElasticParameters::BulkShear { .. } => {
                let rk = relative(at, a)?;
                let rg = relative(bt, b)?;
                // q=3K+G can overflow even when its reciprocal and the native
                // isotropic response are representable. Normalize first.
                let scale = a.max(b);
                let k = a / scale;
                let g = b / scale;
                let q = finite(3.0 * k + g)?;
                let rq = finite((3.0 * k / q) * rk + (g / q) * rg)?;
                // c = 1/(4G) + 3/(4q). Their normalized compliance weights
                // avoid evaluating either small reciprocal separately.
                let ratio = finite(3.0 * g / q)?;
                let second = ratio / (1.0 + ratio);
                let first = 1.0 - second;
                [
                    finite(-first * rg - second * rq)?,
                    finite(2.0 * (first * rg * rg + second * rq * rq))?,
                ]
            }
        };
        Ok(result)
    }
}

/// Original two-solid compliance on a complete rectangular temperature domain.
/// The source laws are immutable; only actual surface temperatures evolve.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalHertzPair {
    surfaces: std::sync::Arc<[Surface; 2]>,
    identity: ContentHash,
}

/// Native queried states and exact analytical derivatives of the ORIGINAL
/// reduced modulus, in ordered surface-A/surface-B temperature coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalHertz {
    pub surfaces: [IsotropicSolidStatePoint; 2],
    pub reduced_modulus: ResolvedReducedModulus,
    pub gradient_pa_k: [f64; 2],
    pub hessian_pa_k2: [[f64; 2]; 2],
}
impl ThermalHertzPair {
    pub fn new(
        cards: [&MaterialCard; 2],
        references: [&IsotropicSolidStatePoint; 2],
        temperature_ranges_k: [[f64; 2]; 2],
    ) -> Result<Self, ThermalHertzError> {
        let surfaces = [
            Surface::new(cards[0], references[0], temperature_ranges_k[0])?,
            Surface::new(cards[1], references[1], temperature_ranges_k[1])?,
        ];
        let legacy = surfaces.iter().all(Surface::legacy_identity);
        let mut hash = DomainHasher::new(if legacy {
            "org.frankensim.fs-material.thermal-hertz-pair.v1"
        } else {
            "org.frankensim.fs-material.thermal-hertz-pair.v2"
        });
        for surface in &surfaces {
            hash.update(surface.card.content_hash().as_bytes());
            hash.update(surface.reference.identity().as_bytes());
            for value in surface
                .temperature_range_k
                .into_iter()
                .chain([surface.source_slopes_per_k[0]])
            {
                hash.update(&value.to_bits().to_le_bytes());
            }
            if !legacy {
                hash.update(&surface.source_slopes_per_k[1].to_bits().to_le_bytes());
            }
        }
        let result = Self {
            surfaces: std::sync::Arc::new(surfaces),
            identity: hash.finalize(),
        };
        for a in temperature_ranges_k[0] {
            for b in temperature_ranges_k[1] {
                result.at_temperatures([a, b])?;
            }
        }
        Ok(result)
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub fn card_identities(&self) -> [ContentHash; 2] {
        self.surfaces
            .each_ref()
            .map(|surface| surface.card.content_hash())
    }
    pub fn reference_temperatures_k(&self) -> [f64; 2] {
        self.surfaces
            .each_ref()
            .map(|surface| surface.reference_temperature_k)
    }
    pub fn temperature_ranges_k(&self) -> [[f64; 2]; 2] {
        self.surfaces
            .each_ref()
            .map(|surface| surface.temperature_range_k)
    }
    pub fn limiting_yield_stress_pa(&self) -> f64 {
        self.surfaces[0]
            .reference
            .yield_stress_pa()
            .min(self.surfaces[1].reference.yield_stress_pa())
    }

    pub fn at_temperatures(
        &self,
        temperatures_k: [f64; 2],
    ) -> Result<ResolvedThermalHertz, ThermalHertzError> {
        let surfaces = [
            self.surfaces[0].at_temperature(temperatures_k[0])?,
            self.surfaces[1].at_temperature(temperatures_k[1])?,
        ];
        let reduced_modulus = surfaces[0].reduced_modulus_against(&surfaces[1])?;
        let e = reduced_modulus.value_pa;
        let mut gradient_pa_k = [0.0; 2];
        let hessian_pa_k2;
        if self.surfaces.iter().all(Surface::fixed_poisson) {
            // Preserve the original arithmetic in the fixed-Poisson limit.
            let mut weights = [0.0; 2];
            let mut ratios = [0.0; 2];
            for i in 0..2 {
                let bulk = surfaces[i].young_modulus_pa();
                weights[i] = finite(
                    ((1.0 - surfaces[i].poisson_ratio() * surfaces[i].poisson_ratio()) / bulk) * e,
                )?;
                ratios[i] = finite(self.surfaces[i].source_slopes_per_k[0] / bulk)?;
                gradient_pa_k[i] = finite(finite(e * weights[i])? * ratios[i])?;
                if self.surfaces[i].source_slopes_per_k[0] != 0.0
                    && (ratios[i] == 0.0 || gradient_pa_k[i] == 0.0)
                {
                    return Err(ThermalHertzError::Invalid(
                        "nonzero source derivative underflows",
                    ));
                }
            }
            // Concavity of harmonic compliance gives a negative-semidefinite
            // Hessian. The factored form avoids subtracting nearly equal terms.
            let scale = finite(finite(2.0 * e * weights[0])? * weights[1])?;
            let cross = finite(finite(scale * ratios[0])? * ratios[1])?;
            hessian_pa_k2 = [
                [finite(-finite(scale * ratios[0])? * ratios[0])?, cross],
                [cross, finite(-finite(scale * ratios[1])? * ratios[1])?],
            ];
        } else {
            let mut normalized_first = [0.0; 2];
            let mut normalized_second = [0.0; 2];
            for i in 0..2 {
                let response = &surfaces[i];
                let compliance = (1.0 - response.poisson_ratio() * response.poisson_ratio())
                    / response.young_modulus_pa();
                let weight = finite(compliance * e)?;
                let [first, second] = self.surfaces[i].compliance_ratios(response)?;
                normalized_first[i] = finite(weight * first)?;
                normalized_second[i] = finite(weight * second)?;
                gradient_pa_k[i] = finite(-e * normalized_first[i])?;
                if (first != 0.0 && normalized_first[i] == 0.0)
                    || (second != 0.0 && normalized_second[i] == 0.0)
                    || (normalized_first[i] != 0.0 && gradient_pa_k[i] == 0.0)
                {
                    return Err(ThermalHertzError::Invalid(
                        "nonzero contact thermal derivative underflows",
                    ));
                }
            }
            // E*=1/(c_a+c_b): d_i E*=-E*² c_i' and
            // d_ij E*=2 E*³ c_i' c_j' - delta_ij E*² c_i''.
            // A general affine parameter pair need not give a concave E*.
            // The coupled thermal owner still admits the total heat-capacity
            // matrix, including its caloric terms, before publishing a step.
            let cross = finite(finite(2.0 * e * normalized_first[0])? * normalized_first[1])?;
            hessian_pa_k2 = [
                [
                    finite(
                        e * finite(
                            2.0 * normalized_first[0] * normalized_first[0] - normalized_second[0],
                        )?,
                    )?,
                    cross,
                ],
                [
                    cross,
                    finite(
                        e * finite(
                            2.0 * normalized_first[1] * normalized_first[1] - normalized_second[1],
                        )?,
                    )?,
                ],
            ];
        }
        Ok(ResolvedThermalHertz {
            surfaces,
            reduced_modulus,
            gradient_pa_k,
            hessian_pa_k2,
        })
    }
}

#[cfg(test)]
mod tests;
